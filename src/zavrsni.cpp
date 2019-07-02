#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <dirent.h>
#include <set>
#include <map>
#include <unordered_map>
#include <iterator>
#include "spoa/spoa.hpp"
#include "blue_alignment.hpp"
#include "bioparser/bioparser.hpp"

class InputFile {
    public:
        std::string name;
        std::string sequence;
        std::string quality;

     InputFile(
       const char* name_, uint32_t name_length_,
       const char* sequence_, uint32_t sequence_length_,
       const char* quality_, uint32_t quality_length_
    ) :
         name(name_, name_length_),
         sequence(sequence_, sequence_length_),
         quality(quality_, quality_length_)

     { }

     InputFile(
       const char* name_, uint32_t name_length_,
        const char* sequence_, uint32_t sequence_length_
    ) :
         name(name_, name_length_),
         sequence(sequence_, sequence_length_)
     { }
};

typedef struct {
    std::string seq_name;
    std::string seq;
} Sequence;

typedef struct {
	std::vector<Sequence> sequences;
	std::string consensus;
} Cluster;

typedef int (*diffFunction)(std::string& sequence1, std::string& sequence2);
int calculateDifference1(std::string& sequence1, std::string& sequence2)
{
	int diff = 0;
	for (int i = 0, length = sequence1.size(); i < length; i++)
	{
		if (sequence1[i] != sequence2[i]) ++diff;
	}
	return diff;
}

int calculateDifference2(std::string& sequence1, std::string& sequence2)
{
    int size1 = sequence1.size();
    int size2 = sequence2.size();
    std::string cigar;
    unsigned int target_begin;
    int score = blue::pairwise_alignment(sequence1.c_str(), size1, sequence2.c_str(), size2, blue::getType("local"), 1, 0, 0, cigar, target_begin);
    int smaller_seq = size1 < size2 ? size1 : size2;
    return smaller_seq - score;
}

std::vector<Cluster> clustering( std::vector<Sequence>& msas, int difference, diffFunction calculateDifference)
{
    std::vector<Cluster> clusters;
	int clustersSize = 0;
	for (int i = msas.size()-1; i>=0; i--)
	{
 		bool isNewCluster = true;
		for (int j = 0, length = clustersSize; j < length; j++)
		{
			if (length == 0) {
				Cluster first;
				first.sequences.push_back(msas[i]);
				clusters.push_back(first);
				++clustersSize;
				break;
			}

            bool inside = true;
			for(int k = 0; k < clusters[j].sequences.size(); k++)
			{
                int distance = calculateDifference(msas[i].seq, clusters[j].sequences[k].seq);
                if(distance > difference)
                {
                    inside = false;
                    break;
                }
            }

			if (inside)
			{
				clusters[j].sequences.push_back(msas[i]);
				isNewCluster = false;
				break;
			}
		}
		if (isNewCluster)
		{
			Cluster cluster;
			cluster.sequences.push_back(msas[i]);
			clusters.push_back(cluster);
			++clustersSize;
		}
    }
	return clusters;
}

void makeClustering(int numberOfTimes, std::vector<std::string>& consensuses);

int main(int argc, char** argv) {
	if (argc != 3) {
		std::cout << "Niste predali putanju do .msa datoteke ili do direktorija s fastq uzorcima" << std::endl;
		exit(1);
	}

    std::string msa = argv[1];
    std::string file = argv[2];

    std::vector<std::unique_ptr<InputFile>> fastq_objects;

    auto fastq_parser = bioparser::createParser<bioparser::FastqParser, InputFile>(file);
    uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
    while (true) {
        auto status = fastq_parser->parse(fastq_objects, size_in_bytes);
        if (status == false) {
            break;
        }
    }
    std::cout << "Uzorak je parsiran." << std::endl;

	std::ifstream infile(msa);
    std::vector<Sequence> msas;
    std::string line;

    //first line
    std::getline(infile, line);
    while (std::getline(infile, line))
    {
        Sequence seq;
        seq.seq_name = line;
        std::getline(infile, line);
        seq.seq = line;
        msas.push_back(seq);
    }

    std::unordered_map<std::string, std::string> name_seq;
    for(auto& i : fastq_objects)
    {
        name_seq[i->name] = i->sequence;
    }

    std::vector<Cluster> clusters = clustering(msas, 12, &calculateDifference1);
    std::vector<Sequence> consensuses;
    for(int i = 0; i < clusters.size(); i++)
    {
        if(clusters[i].sequences.size() >= 10)
        {
            auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(0), 5, -4, -5, -6);
            auto graph = spoa::createGraph();

            for (const auto& it: clusters[i].sequences) {
                std::string name;
                std::unordered_map<std::string, std::string>::iterator iter = name_seq.find(it.seq_name);
                if(iter != name_seq.end())
                {
                   name = iter->second;
                }

                auto alignment = alignment_engine->align(name, graph);
                graph->add_alignment(alignment, name);
            }

            clusters[i].consensus = graph->generate_consensus();
            consensuses.push_back({"",clusters[i].consensus});

        }
    }

    std::cout << "Pronađeni aleli:" << std::endl << std::endl;

    clusters = clustering(consensuses, 6, &calculateDifference2);
    for(int i = 0; i < clusters.size(); i++)
    {
        auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(0), 5, -4, -5, -6);
        auto graph = spoa::createGraph();

        for (const auto& it : clusters[i].sequences) {
            auto alignment = alignment_engine->align(it.seq, graph);
            graph->add_alignment(alignment, it.seq);
        }

        clusters[i].consensus = graph->generate_consensus();

        fprintf(stderr, "Consensus (%zu)\n",  clusters[i].consensus.size());
        fprintf(stderr, "%s\n\n", clusters[i].consensus.c_str());
    }
	return 0;
}

