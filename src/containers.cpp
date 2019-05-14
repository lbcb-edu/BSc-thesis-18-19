#include <fstream>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <string>
#include "bioparser/bioparser.hpp"

class InputPafFile{
    public:
        std::string seq_name;
        uint32_t seq_length;
        uint32_t seq_start_pos;
        uint32_t seq_end_pos;
        char align_orientation;
        std::string gen_name;
        uint32_t gen_length;
        uint32_t gen_start_pos;
        uint32_t gen_end_pos;
        uint32_t matches;
        uint32_t alignment_length;
        uint32_t mapping_quality;

    InputPafFile (const char* q_name, uint32_t q_name_length,
    uint32_t q_length,
    uint32_t q_begin,
    uint32_t q_end,
    char orientation,
    const char* t_name, std::uint32_t t_name_length,
    uint32_t t_length,
    uint32_t t_begin,
    uint32_t t_end,
    uint32_t matching_bases,
    uint32_t overlap_length,
    uint32_t mapping_quality_) :

        seq_name(q_name, q_name_length),
        seq_length(q_length),
        seq_start_pos(q_begin),
        seq_end_pos(q_end),
        align_orientation(orientation),
        gen_name(t_name, t_name_length),
        gen_length(t_length),
        gen_start_pos(t_begin),
        gen_end_pos(t_end),
        matches(matching_bases),
        alignment_length(overlap_length),
        mapping_quality(mapping_quality_)

    { }
};

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

bool cmp_by_length(std::unique_ptr<InputPafFile>& seq1, std::unique_ptr<InputPafFile>& seq2) {
    if (seq1->seq_name == seq2->seq_name) {
        return seq1->seq_start_pos < seq2->seq_start_pos;
    }
    return seq1->seq_length > seq2->seq_length;

}

int main(int argc, char* argv[]){
    std::string file_name(argv[optind]);
    std::string first(argv[optind+1]);

    std::ifstream infile(file_name);

    //std::vector<my_read> reads;
    /*reads.push_back(my_read());
    int i = 0;*/

    std::vector<std::unique_ptr<InputPafFile>> paf_objects;
    auto paf_parser = bioparser::createParser<bioparser::PafParser, InputPafFile>(file_name);
    auto status = paf_parser->parse_objects(paf_objects, -1);

    std::cout << "Velicina prije sorta " << paf_objects.size() << std::endl;

    std::sort(paf_objects.begin(), paf_objects.end(), cmp_by_length);

    std::cout << "Velicina nakon sorta " << paf_objects.size() << std::endl;

   /* while(infile >> reads[i].seq_name >> reads[i].seq_length >> reads[i].seq_start_pos >> reads[i].seq_end_pos >> reads[i].gen_start_pos >> reads[i].gen_end_pos){
        reads.push_back(my_read());
        i++;
    }

    for(int i = 0; i < reads.size()-1; i++){
        for(int j = i + 1; j < reads.size(); j++){
            if(reads[i].gen_start_pos <= reads[j].gen_start_pos && reads[i].gen_end_pos >= reads[j].gen_end_pos && \
                ((float)(reads[j].seq_end_pos - reads[j].seq_start_pos)/(float)reads[j].seq_length)>0.9) {
                contained.insert(reads[j].seq_name);
            }
        }
    }*/

    /*for (auto &i : paf_objects) {
        std::cout << i->seq_name << "\t" << i->seq_length << "\t" << i->seq_start_pos << "\t" << i->seq_end_pos << std::endl;
    }*/

    for(int i = 0; i < paf_objects.size() - 1; i++) {
        std::unordered_set<std::string> contained;
        for (int j = i + 1; j < paf_objects.size(); j++) {
            /*if ((paf_objects[j])->seq_name == "b5f24896-7153-438e-bebd-ec98d4bd6004_Basecall_2D_000_2d" && (paf_objects[i])->seq_name == "974b7693-9b5e-4610-bbc2-0d7a5ca9d370_Basecall_2D_000_2d") {
                std::cout << "Gledam po j sekvencu " << (paf_objects[j])->seq_name << " duljine mapiranja " << (paf_objects[j])->seq_start_pos << "\t" << (paf_objects[j])->seq_end_pos << std::endl;
            }*/
            if((paf_objects[i])->gen_start_pos <= (paf_objects[j])->gen_start_pos && (paf_objects[i])->gen_end_pos >= (paf_objects[j])->gen_end_pos && \
                ((float)((paf_objects[j])->seq_end_pos - (paf_objects[j])->seq_start_pos))/(float)(paf_objects[j])->seq_length>0.9){
                    /*if ((paf_objects[j])->seq_name == "b5f24896-7153-438e-bebd-ec98d4bd6004_Basecall_2D_000_2d"){
                        std::cout << "Trenutni i " << (paf_objects[i])->seq_name << "\t" << (paf_objects[i])->seq_length << std::endl;
                        std::cout << "Sad sam nasao " << (paf_objects[j])->seq_name << std::endl;
                    }*/
                    //std::cout << "prije " << (paf_objects[j])->seq_name << "\t" << (paf_objects[j])->gen_start_pos << "\t" << (paf_objects[j])->gen_end_pos << std::endl;
                    std::string contained_seq = (paf_objects[j])->seq_name;
                    paf_objects.erase(paf_objects.begin() + j);
                    while ((paf_objects[j])->seq_name == contained_seq) {
                        paf_objects.erase(paf_objects.begin() + j);
                    }
                    contained.insert((paf_objects[j])->seq_name);
                    //std::cout << "poslije " << (paf_objects[j])->seq_name << "\t" << (paf_objects[j])->gen_start_pos << "\t" << (paf_objects[j])->gen_end_pos << std::endl;
                    //std::cout << "Velicina nakon brisanja " << paf_objects.size() << std::endl;
            } else {
                j++;
            }
        }
    }

    for (auto &i : paf_objects) {
        std::cout << i->seq_name << "\t" << i->seq_length << "\t" << i->seq_start_pos << "\t" << i->seq_end_pos << std::endl;
    }

    std::cout << "Velicina paf objectsa nakon " << paf_objects.size() << std::endl;

    infile.close();

    /*std::vector<std::unique_ptr<InputFile>> first_object;
    if (first.find("fastq")>first.length() && first.find("fq")>first.length()){

        auto fasta_parser = bioparser::createParser<bioparser::FastaParser, InputFile>(first);
        fasta_parser->parse_objects(first_object, -1);

        std::cout << "We've parsed file." << std::endl;

    } else {
        auto fastq_parser = bioparser::createParser<bioparser::FastqParser, InputFile>(first);
        uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
            while (true) {
                auto status = fastq_parser->parse_objects(first_object, size_in_bytes);

                if (status == false) {
                    break;
                }
            }
        std::cout << "We've parsed file." << std::endl;
    }

    std::ofstream outfile;
    outfile.open("reduced.fasta");

    for(auto& i : first_object) {
        if(contained.find(i->name) == contained.end()) {
            outfile << ">" << i->name << "\n";
            outfile << i->sequence << "\n";
        }
    }

    outfile.close();*/

    return 0;

}

