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
//to do: dodati help, ovisnost o parametru k, mozda u calcDifference dodati tablicu penaliziranja razl parametara

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

/*std::vector<Cluster> clustering( std::vector<Sequence>& msas, int difference, diffFunction calculateDifference)
{
    std::vector<Cluster> clusters;
	int clustersSize = 0;
	for (int i = 0, size = msas.size(); i < size; i++)
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
}*/


std::vector<Cluster> clustering( std::vector<Sequence>& msas, int difference, diffFunction calculateDifference)
{
    std::vector<Cluster> clusters;
	int clustersSize = 0;
	for (int i = msas.size()-1; i>=0; i--)
	{
		bool isNewCluster = true;
		std::pair<double, int> minDistance = std::make_pair(300.0, 0);

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
            double averageDistance = 0.0;
			for(int k = 0; k < clusters[j].sequences.size(); k++)
			{
                int distance = calculateDifference(msas[i].seq, clusters[j].sequences[k].seq);
                if(distance > difference)
                {
                    inside = false;
                    break;
                } else {
                    averageDistance += distance;
                }
            }

            averageDistance /= (double) clusters[j].sequences.size();

			if (inside)
			{
                if(averageDistance < minDistance.first) {
                    minDistance = std::make_pair(averageDistance, j);
                }
				isNewCluster = false;
			}
		}

		if (isNewCluster)
		{
			Cluster cluster;
			cluster.sequences.push_back(msas[i]);
			clusters.push_back(cluster);
			++clustersSize;
		} else {
            clusters[minDistance.second].sequences.push_back(msas[i]);
		}
	}
	return clusters;
}

inline bool ends_with(std::string const& value, std::string const& ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

int main(int argc, char** argv) {
	if (argc != 3) {
	//	std::cout << "Niste predali putanju do .msa datoteke" << std::endl;
	//	exit(1);
	}

	std::string parentDir = "/home/sanja/Desktop/zavrsniRad/rezultati/Verzije/2.2/spoa_folder/";
	std::string fastq = "/home/sanja/Desktop/zavrsniRad/za_skosier2/some_fastq2/";
	//std::string parentDir = "/home/sanja/Desktop/spoa_folder/";

    DIR *dir;
    struct dirent *ent;
    std::map<std::string, std::string> fastqPaths;
    if ((dir = opendir (fastq.c_str())) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            std::string file = fastq + ent->d_name;
            if(file.size() - fastq.size() < 3) continue;
            std::string name(ent->d_name);

            size_t lastindex = name.find_last_of(".");
            std::string rawname = name.substr(0, lastindex);

            fastqPaths.insert({rawname, file});
        }
    }

    int samples = 0;
    std::vector<int> noOfAllels;
    std::map <std::string, int> allels;
    std::map <std::string, std::vector<std::string>> uzorakAleli;
    std::map <std::string, std::vector<std::string>> alelUzorci;

    if ((dir = opendir (parentDir.c_str())) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            std::string msa = parentDir + ent->d_name;
            if(ends_with(msa, ".consensus")) continue;
            std::string name(ent->d_name);
            size_t lastindex = name.find_last_of(".");
            std::string rawname = name.substr(0, lastindex);

            std::cout << "Trenutno na: " << rawname << std::endl;
            std::map<std::string, std::string>::iterator it = fastqPaths.find(rawname);
            if (it == fastqPaths.end()) {
                std::cout << "Ne postoji pripadajući " << rawname << ".fastq zapis" << std::endl;
                continue;
            }

            ++samples;

            std::string file = fastqPaths.find(rawname)->second;
            std::vector<std::unique_ptr<InputFile>> fastq_objects;
            auto fastq_parser = bioparser::createParser<bioparser::FastqParser, InputFile>(file);
            uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
            while (true) {
                auto status = fastq_parser->parse(fastq_objects, size_in_bytes);
                if (status == false) {
                    break;
                }
            }

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

            clusters = clustering(consensuses, 6, &calculateDifference2);
            noOfAllels.push_back(clusters.size());
            std::vector<std::string> aleliUzorka;

            for(int i = 0; i < clusters.size(); i++)
            {
                std::cout << clusters[i].sequences.size() << std::endl;

                auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(0), 5, -4, -5, -6);
                auto graph = spoa::createGraph();

                for (const auto& it : clusters[i].sequences) {
                    auto alignment = alignment_engine->align(it.seq, graph);
                    graph->add_alignment(alignment, it.seq);
                }

                std::string consensus = graph->generate_consensus();
                allels[consensus]++;
                aleliUzorka.push_back(consensus);

                std::map<std::string, std::vector<std::string>>::iterator iter = alelUzorci.find(consensus);
                if(iter != alelUzorci.end())
                {
                    alelUzorci[consensus].push_back(rawname);
                } else {
                    std::vector<std::string> uzorak;
                    uzorak.push_back(rawname);
                    alelUzorci[consensus] = uzorak;
                }

                std::cout << "Sad sam ispred alel uzorci "<< std::endl;
                std::cout << alelUzorci[consensus].size() << std::endl;

                alelUzorci[consensus];
                std::cout << alelUzorci[consensus].size() << std::endl;

                std::cout << "Consensus " << consensus.size() << std::endl;
                std::cout << consensus << std::endl;
            }
            uzorakAleli[rawname] = aleliUzorka;
        }
    }

    std::cout << "samples " << samples << std::endl;
    for (int i : noOfAllels) std::cout << i << " ";
    std::cout << "" << std::endl;

    int sum = 0;
    std::map<std::string, int>::iterator it = allels.begin();
    std::cout << "\n\nDivokoze: " << std::endl;
	while (it != allels.end())
	{
		std::string consensus = it->first;
        int count = it->second;
        std::cout << count << " --- " << consensus << std::endl;
        it++;
        sum += count;
	}

    std::map<std::string, std::vector<std::string>>::iterator it2 = alelUzorci.begin();
    std::cout << "\n\n alel uzorak: " << std::endl;
	while (it2 != alelUzorci.end())
	{
		std::string consensus = it2->first;
        std::vector<std::string> count = it2->second;
        std::cout << consensus << std::endl;
        for(std::string s : count) std::cout << s << " ";
        std::cout << " " << std::endl;
        it2++;
	}

	std::map<int, int> alelaImaTolkoUzoraka;
	std::map<std::string, std::vector<std::string>>::iterator it3 = uzorakAleli.begin();
    std::cout << "\n\n uzorak alel: " << std::endl;
	while (it3 != uzorakAleli.end())
	{
		std::string consensus = it3->first;
        std::vector<std::string> count = it3->second;
        std::cout << consensus << std::endl;
        alelaImaTolkoUzoraka[count.size()]++;
        for(std::string s : count) std::cout << s << " ";
        std::cout << " " << std::endl;
        it3++;
	}

	std::cout << "broj razl alela: " << noOfAllels.size() << std::endl;
	std::cout << "broj pronadenih alela " << sum << std::endl;

	std::map<int, int>::iterator it4 = alelaImaTolkoUzoraka.begin();
	while (it4 != alelaImaTolkoUzoraka.end())
	{
		std::cout << it4->first << "," << it4->second << std::endl;
		it4++;
	}
	return 0;
}
