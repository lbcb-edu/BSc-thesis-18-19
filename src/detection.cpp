#include "bioparser/bioparser.hpp"
#include <unistd.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <string>
#include <tuple>
#include <algorithm>
#include <stdlib.h>

struct mapping{
    uint32_t seq_len;
    uint32_t seq_start;
    uint32_t seq_end;
    char relative_strand;
    uint32_t gen_start;
    uint32_t gen_end;
    uint32_t gen_length;
    std::string ref_name;

    mapping(uint32_t seq_l, uint32_t seq_s, uint32_t seq_e, char relative_s, uint32_t gen_s, uint32_t gen_e, uint32_t gen_l, std::string ref_n) :
        seq_len(seq_l),
        seq_start(seq_s),
        seq_end(seq_e),
        relative_strand(relative_s),
        gen_start(gen_s),
        gen_end(gen_e),
        gen_length(gen_l),
        ref_name(ref_n)
    { }
};

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

bool sort_function(mapping m1, mapping m2) {
    if (m1.seq_start == m2.seq_start) {
        return m1.seq_end > m2.seq_end;
    }
    return (m1.seq_start < m2.seq_start);
}

bool sort_paf_function(std::unique_ptr<InputPafFile>& m1, std::unique_ptr<InputPafFile>& m2) {
    return (m1->gen_start_pos < m2->gen_start_pos);
}

bool sort_reference_repeatings(std::tuple<uint32_t, uint32_t>& rep1, std::tuple<uint32_t, uint32_t>& rep2) {
    if (std::get<0>(rep1) != std::get<0>(rep2)) {
        return (std::get<0>(rep1) < std::get<0>(rep2));
    } else {
        return (std::get<1>(rep1) > std::get<1>(rep2));
    }
}

bool unique_read_function(mapping& m1, mapping& m2) {
    return (m1.seq_start == m2.seq_start && m1.seq_end == m2.seq_end);
}

bool paf_unique(const std::unique_ptr<InputPafFile>& a, const std::unique_ptr<InputPafFile>& b) {
	return a->seq_name == b->seq_name;
}

bool unique_reference_repeatings(std::tuple<uint32_t, uint32_t>& rep1, std::tuple<uint32_t, uint32_t>& rep2) {
    return std::get<0>(rep1) == std::get<0>(rep2);
}

bool paf_cmp(const std::unique_ptr<InputPafFile>& a, const std::unique_ptr<InputPafFile>& b) {
		if (a->gen_name == b->gen_name) {
	 		if (a->gen_start_pos == b->gen_start_pos) {
	        	return (a->gen_end_pos > b->gen_end_pos);
	        }
	        return (a->gen_start_pos < b->gen_start_pos);
	    }
	    return (a->gen_name > b->gen_name);
}

bool repeats_by_length(mapping& m1, mapping& m2) {
    return ((m1.gen_end - m1.gen_start) > (m2.gen_end - m2.gen_start));
}

// Function finds repeating reference regions by iterating through found repeats and checking if any two repeats map at same region
std::unordered_map<std::string, std::vector<std::tuple<uint32_t, uint32_t>>> find_repeating_regions(std::vector<mapping>& repeats) {
    std::unordered_map<std::string, std::vector<std::tuple<uint32_t, uint32_t>>> return_map;
    for (int i = 0; i < repeats.size() - 1; i++) {
        for (int j = i + 1;  j < repeats.size(); j++) {
            if (repeats[i].ref_name == repeats[j].ref_name) {
                if (repeats[i].gen_start <= repeats[j].gen_start && repeats[i].gen_end >= repeats[j].gen_end && std::find((return_map[repeats[i].ref_name]).begin(),
                (return_map[repeats[i].ref_name]).end(), std::make_tuple(repeats[i].gen_start, repeats[i].gen_end)) == (return_map[repeats[i].ref_name]).end()) {
                    if (return_map.find(repeats[j].ref_name) == return_map.end()){
                        std::vector<std::tuple<uint32_t, uint32_t>> vec;
                        vec.push_back(std::make_tuple(repeats[j].gen_start, repeats[j].gen_end));
                        return_map[repeats[j].ref_name] = vec;
                    } else {
                        (return_map[repeats[j].ref_name]).push_back(std::make_tuple(repeats[j].gen_start, repeats[j].gen_end));
                    }
                    break;
                }
            }
        }
    }

    return return_map;
}


int main(int argc, char* argv[]) {
    std::string paf_file(argv[optind]);
    std::string fasta_file(argv[optind+1]);

    std::vector<std::unique_ptr<InputPafFile>> paf_objects;
    auto paf_parser = bioparser::createParser<bioparser::PafParser, InputPafFile>(paf_file);
    auto status = paf_parser->parse_objects(paf_objects, -1);

    std::unordered_map<std::string, std::vector<mapping>> sequence_mapping_details;

    // Unordered map of read mappings is formed
    for (auto& i : paf_objects){
        std::string key(i->seq_name);
        mapping mapp(i->seq_length, i->seq_start_pos, i->seq_end_pos, i->align_orientation, i->gen_start_pos, i->gen_end_pos, i->gen_length, i->gen_name);
        if (sequence_mapping_details.find(key) == sequence_mapping_details.end()){
            std::vector<mapping> vec;
            vec.push_back(mapp);
            sequence_mapping_details[key] = vec;
        } else {
            std::vector<mapping> current_vector = sequence_mapping_details[key];
            if ((mapp.seq_end - mapp.seq_start) > (current_vector[0].seq_end - current_vector[0].seq_start)) {
                auto it = current_vector.begin();
                current_vector.insert(it, mapp);
            } else {
                current_vector.push_back(mapp);
            }
            sequence_mapping_details[key] = current_vector;
        }
    }



    std::unordered_map<std::string, std::vector<mapping>>::iterator itr;

    std::unordered_map<std::string, std::vector<mapping>> chimeric_reads;
    std::unordered_map<std::string, std::vector<mapping>> repeating_reads;
    std::unordered_set<std::string> chimers;
    std::unordered_set<std::string> repeatings;
    std::unordered_set<std::string> regular;
    std::vector<mapping> all_repeatings;

    /* By checking reads that have multiple mappings and positions of those mappings chimeric and repeating reads
     * are being detected and separated into sets
    */
    for (itr = sequence_mapping_details.begin(); itr != sequence_mapping_details.end(); itr++) {
        if ((itr->second).size() > 1) {
            std::sort((itr->second).begin(), (itr->second).end(), sort_function);
            uint32_t seq_start = (itr->second)[0].seq_start;
            uint32_t seq_end = (itr->second)[0].seq_end;
            uint32_t gen_start = (itr->second)[0].gen_start;
            uint32_t gen_end = (itr->second)[0].gen_end;
            std::string ref_n = (itr->second)[0].ref_name;
            char rel_strand = (itr->second)[0].relative_strand;
            for (int i = 1 ; i < (itr->second).size(); i++) {
                if (!(seq_start <= (itr->second)[i].seq_start && seq_end >= (itr->second)[i].seq_end) || (((itr->second)[i].seq_start - seq_start) > 20/100.0*(seq_end - seq_start) &&
                seq_end < (itr->second)[i].seq_end) ||
                ((seq_end - (itr->second)[i].seq_end) > 0.2*(seq_end - seq_start) && seq_start > (itr->second)[i].seq_start)) {
                    if (((itr->second)[i].seq_end - (itr->second)[i].seq_start) > 0.9*(itr->second)[i].seq_len || (seq_end - seq_start) > 0.9*(itr->second)[i].seq_len) {
                        regular.emplace(itr->first);
                        break;
                    }
                    chimers.emplace(itr->first);
                    int seq_gap = abs((int)((itr->second)[i].seq_start - (itr->second)[i - 1].seq_end));
                    int ref_gap = abs((int)((itr->second)[i].gen_start - (itr->second)[i - 1].gen_end));
                    if((std::min(gen_start, (itr->second)[i].gen_start) > 100 || std::max(gen_end, (itr->second)[i].gen_end) < ((itr->second)[i].gen_length - 100))
                    && !(((float)(std::max(ref_gap, seq_gap) - std::min(ref_gap, seq_gap))/std::max(ref_gap, seq_gap)) < (20.0/100) && (itr->second)[i].ref_name == (itr->second)[i - 1].ref_name && (itr->second)[i].relative_strand == (itr->second)[i].relative_strand)){
                        chimeric_reads[itr->first] = itr->second;
                        break;
                    }
                }
            }
            if (chimers.find(itr->first) == chimers.end() && regular.find(itr->first) == regular.end()) {
                int i = 1;
                if ((itr->second)[0].seq_start == (itr->second)[1].seq_start && (itr->second)[0].seq_end == (itr->second)[1].seq_end) {
                    i = 0;
                }
                all_repeatings.insert(all_repeatings.end(), (itr->second).begin() + i, (itr->second).end());
                repeating_reads[itr->first] = itr->second;
            }
        }
    }

    std::sort(all_repeatings.begin(), all_repeatings.end(), repeats_by_length);

    std::unordered_map<std::string, std::vector<std::tuple<uint32_t, uint32_t>>> reference_repeating_regions = find_repeating_regions(all_repeatings);

    std::unordered_set<std::string> contained_repeats;

    /* After detecting repeating reference regions all found repeatings are being separated in two groups; those who are contained in
     * some reference repeating region and those who are not
    */
    for (auto rep_it = repeating_reads.begin(); rep_it != repeating_reads.end(); rep_it++) {
        bool contained = false;
        int n = 0;
        do {
            mapping longest_mapping = (repeating_reads[rep_it->first])[n];
            for (auto regions_it = reference_repeating_regions.begin(); regions_it != reference_repeating_regions.end(); regions_it++) {
                if ((repeating_reads[rep_it->first])[n].ref_name == regions_it->first) {
                    std::vector<std::tuple<uint32_t, uint32_t>> rep_regions = reference_repeating_regions[regions_it->first];
                    for (int j = 0; j < rep_regions.size(); j++) {
                        if ((repeating_reads[rep_it->first])[n].gen_start >= std::get<0>(rep_regions[j]) && (repeating_reads[rep_it->first])[n].gen_end <= std::get<1>(rep_regions[j])) {
                            contained = true;
                            break;
                        }
                    }
                }
                if (contained) {
                    contained_repeats.emplace(rep_it->first);
                    break;
                }
            }
            n++;
        } while ((repeating_reads[rep_it->first])[n].seq_start == (repeating_reads[rep_it->first])[n-1].seq_start && (repeating_reads[rep_it->first])[n].seq_end == (repeating_reads[rep_it->first])[n-1].seq_end && !contained);
        if (!contained) {
            repeatings.emplace(rep_it->first);
        }
    }


    std::unordered_map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> split_chimers;

    /* In order to use chimeric reads for the assembly chimeric reads are being splitted into multiple reads by their breakpoints;
     * special attention is given to those chimeric reads who have repeating regions as part of the read; those parts are being
     * treated as repeating reads and annotated according to that
    */
    for (itr = chimeric_reads.begin(); itr != chimeric_reads.end(); itr++) {
        for (int i = 0; i < (itr->second).size();) {
            if (i == (itr->second).size() - 1 ){
                if (split_chimers.find(itr->first) == split_chimers.end()) {
                    std::vector<std::pair<uint32_t, uint32_t>> vec;
                    vec.push_back(std::make_pair((itr->second)[i].seq_start, (itr->second)[i].seq_end));
                    split_chimers[itr->first] = vec;
                } else {
                    (split_chimers[itr->first]).push_back(std::make_pair((itr->second)[i].seq_start, (itr->second)[i].seq_end));
                }
                break;
            }
            if (!((itr->second)[i].seq_start <= (itr->second[i + 1]).seq_start && (itr->second[i]).seq_end >= (itr->second[i + 1]).seq_end)) {
                uint32_t start_chimer = (itr->second[i]).seq_start;
                uint32_t end_chimer = (itr->second[i]).seq_end;
                uint32_t gen_endd = (itr->second[i]).gen_end;
                char strand = (itr->second[i]).relative_strand;
                std::string reff = (itr->second[i]).ref_name;
                uint32_t seq_diff = abs((int)((itr->second[i + 1]).seq_start - end_chimer));
                uint32_t gen_diff = abs((int)((itr->second[i + 1]).gen_start - gen_endd));
                while (i < ((itr->second).size() - 1) && (std::max(seq_diff, gen_diff) - std::min(seq_diff, gen_diff))/(float)(std::min(seq_diff, gen_diff)) < 0.2
                && (itr->second[i + 1]).relative_strand == strand && (itr->second[i + 1]).ref_name == reff){
                    end_chimer = (itr->second[i + 1]).seq_end;
                    gen_endd = (itr->second[i + 1]).gen_end;
                    i += 1;
                    seq_diff = abs((int)((itr->second)[i + 1].seq_start - end_chimer));
                    gen_diff = abs((int)((itr->second)[i + 1].gen_start - gen_endd));
                }
                if (split_chimers.find(itr->first) == split_chimers.end()) {
                    std::vector<std::pair<uint32_t, uint32_t>> vec;
                    vec.push_back(std::make_pair(start_chimer, end_chimer));
                    split_chimers[itr->first] = vec;
                } else {
                    (split_chimers[itr->first]).push_back(std::make_pair(start_chimer, end_chimer));
                }
            } else {
                std::vector<std::tuple<uint32_t, uint32_t>> rep_regions = reference_repeating_regions[(itr->second[i]).ref_name];
                bool in_repeat = false;
                for (int j = 0; j < rep_regions.size(); j++) {
                    if (std::get<0>(rep_regions[j]) <= (itr->second)[i].gen_start && std::get<1>(rep_regions[j]) >= (itr->second)[i].gen_end) {
                        in_repeat = true;
                        break;
                    }
                }
                if (!in_repeat){
                    if (split_chimers.find(itr->first) == split_chimers.end()) {
                        std::vector<std::pair<uint32_t, uint32_t>> vec;
                        vec.push_back(std::make_pair((itr->second)[i].seq_start, (itr->second)[i].seq_end));
                        split_chimers[itr->first] = vec;
                    } else {
                        (split_chimers[itr->first]).push_back(std::make_pair((itr->second)[i].seq_start, (itr->second)[i].seq_end));
                    }
                }
                int l = i;
                while ((itr->second)[l].seq_start <= (itr->second[i + 1]).seq_start && (itr->second[l]).seq_end >= (itr->second[i + 1]).seq_end && i < (itr->second).size()) {
                    i++;
                }
            }
            i++;
        }
    }

    std::cout << "Kimernih sekvenci ima " << chimeric_reads.size() << std::endl;
    std::cout << "Repetitivnih sekvenci ima " << repeatings.size() + contained_repeats.size() << std::endl;

    // Input file with reads is being parsed by bioparser library
    std::vector<std::unique_ptr<InputFile>> first_object;
    if (fasta_file.find("fastq")>fasta_file.length() && fasta_file.find("fq")>fasta_file.length()){

        auto fasta_parser = bioparser::createParser<bioparser::FastaParser, InputFile>(fasta_file);
        fasta_parser->parse_objects(first_object, -1);
    } else {
        auto fastq_parser = bioparser::createParser<bioparser::FastqParser, InputFile>(fasta_file);
        uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
            while (true) {
                auto status = fastq_parser->parse_objects(first_object, size_in_bytes);

                if (status == false) {
                    break;
                }
            }
    }

    std::ofstream cleaned_sequences;
    cleaned_sequences.open("cleaned_sequences.fasta");

    std::ofstream chimeric;
    chimeric.open("chimeric_reads.fasta");

    std::ofstream repeating;
    repeating.open("repeating_reads.fasta");

     /* Iterating through vector of reads, reads are being separated into files according to their clasiffication as chimeric or
      * repeating read or none of the above
     */
     for(auto& i : first_object) {
        if (split_chimers.find(i->name) != split_chimers.end()) {
            std::vector<std::pair<uint32_t, uint32_t>> positions = split_chimers[i->name];
            for(int j = 0; j < positions.size(); j++) {
                uint32_t len = std::get<1>(positions[j]) - std::get<0>(positions[j]);
                std::string cut_piece = (i->sequence).substr(std::get<0>(positions[j]), len);
                std::string seq_name = (std::to_string(std::get<0>(positions[j]))).append(i->name);
                seq_name.append(std::to_string(std::get<1>(positions[j])));
                cleaned_sequences << ">" << seq_name << "\n";
                cleaned_sequences << cut_piece << "\n";
            }
        }
        if(chimeric_reads.find(i->name) != chimeric_reads.end()){
            std::vector<mapping> info = chimeric_reads[i->name];
            chimeric << ">" << i->name << "\t";
            std::sort(info.begin(), info.end(), sort_function);
            std::vector<mapping>::iterator vec_it;
            vec_it = std::unique(info.begin(), info.end(), unique_read_function);
            std::set<std::tuple<uint32_t, uint32_t>> written;
            for (auto& i : info) {
                if (written.find(std::make_tuple(i.seq_start, i.seq_end)) == written.end()) {
                    chimeric << i.seq_start << " " << i.seq_end << "\t";
                    written.emplace(std::make_tuple(i.seq_start, i.seq_end));
                }
            }
            chimeric << "\n";
            chimeric << i->sequence << "\n";
        }
        else if (repeatings.find(i->name) != repeatings.end()) {
            repeating << ">" << i->name << "\t";
            std::vector<mapping> info = repeating_reads[i->name];
            std::sort(info.begin(), info.end(), sort_function);
            int j = 1;
            if (info[0].seq_start == info[1].seq_start && info[0].seq_end == info[1].seq_end) {
                j = 0;
            }
            std::set<std::tuple<uint32_t, uint32_t>> written;
            while (j < info.size()) {
                std::tuple<uint32_t, uint32_t> tup = std::make_tuple(info[j].seq_start, info[j].seq_end);
                if (written.find(tup) == written.end()) {
                    repeating << std::get<0>(tup) << " " << std::get<1>(tup) << "\t";
                    written.emplace(tup);
                }
                j++;
            }
            repeating << "\n";
            repeating << i->sequence << "\n";

            std::string name = (i->name);
            name.append(" ");
            j = 1;
            if (info[0].seq_start == info[1].seq_start && info[0].seq_end == info[1].seq_end) {
                j = 0;
            }
            std::set<std::tuple<uint32_t, uint32_t>> written1;
            while(j < info.size()){
                std::string start_region = "XB:i:";
                std::string end_region = "XE:i:";
                std::tuple<uint32_t, uint32_t> tup = std::make_tuple(info[j].seq_start, info[j].seq_end);
                if(written1.find(tup) == written1.end()){
                    start_region.append(std::to_string(info[j].seq_start));
                    end_region.append(std::to_string(info[j].seq_end));
                    written1.emplace(tup);
                    name.append(start_region);
                    name.append(" ");
                    name.append(end_region);
                    name.append(" ");
                }
                j++;
            }
            cleaned_sequences << ">" << name << "\n";
            cleaned_sequences << i->sequence << "\n";

        }
        else {
            cleaned_sequences << ">" << i->name << "\n";
            cleaned_sequences << i->sequence << "\n";
        }
	if (contained_repeats.find(i->name) != contained_repeats.end()){
	    repeating << ">" << i->name << "\t";
            std::vector<mapping> info = repeating_reads[i->name];
            std::sort(info.begin(), info.end(), sort_function);
            int j = 1;
            if (info[0].seq_start == info[1].seq_start && info[0].seq_end == info[1].seq_end) {
                j = 0;
            }
	    std::set<std::tuple<uint32_t, uint32_t>> written;
            while (j < info.size()) {
                std::tuple<uint32_t, uint32_t> tup = std::make_tuple(info[j].seq_start, info[j].seq_end);
                if (written.find(tup) == written.end()) {
                    repeating << std::get<0>(tup) << " " << std::get<1>(tup) << "\t";
                    written.emplace(tup);
                }
                j++;
            }
            repeating << "\n";
            repeating << i->sequence << "\n";
	}
    }

    return 0;

}
