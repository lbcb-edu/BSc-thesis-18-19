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
    uint32_t seq_start;
    uint32_t seq_end;
    uint32_t gen_start;
    uint32_t gen_end;
    uint32_t gen_length;
    std::string ref_name;

    mapping(uint32_t seq_s, uint32_t seq_e, uint32_t gen_s, uint32_t gen_e, uint32_t gen_l, std::string ref_n) :
        seq_start(seq_s),
        seq_end(seq_e),
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

std::unordered_map<std::string, std::vector<std::tuple<uint32_t, uint32_t>>> find_repeating_regions(std::vector<mapping>& repeats) {
    std::unordered_map<std::string, std::vector<std::tuple<uint32_t, uint32_t>>> return_map;
    for (int i = 0; i < repeats.size() - 1; i++) {
        //std::cout << repeats[i].gen_start << "\t" << repeats[i].gen_end << "\t" << repeats[i].seq_start << "\t" << repeats[i].seq_end << std::endl;
        for (int j = i + 1;  j < repeats.size(); j++) {
            if (repeats[i].ref_name == repeats[j].ref_name) {
                if (repeats[i].gen_start <= repeats[j].gen_start && repeats[i].gen_end >= repeats[j].gen_end && std::find((return_map[repeats[i].ref_name]).begin(),
                (return_map[repeats[i].ref_name]).end(), std::make_tuple(repeats[i].gen_start, repeats[i].gen_end)) == (return_map[repeats[i].ref_name]).end()) {
                    //std::cout << "Na i je " << repeats[i].gen_start << "\t" << repeats[i].gen_end << "\t" << repeats[i].seq_start << "\t" << repeats[i].seq_end << std::endl;
                    //std::cout << "Na j je " << repeats[j].gen_start << "\t" << repeats[j].gen_end << "\t" << repeats[j].seq_start << "\t" << repeats[j].seq_end << "\n" << std::endl;
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

void clear_contained_reads(std::vector<std::unique_ptr<InputPafFile>> &paf_objects) {
	std::vector<std::unique_ptr<InputPafFile>>::iterator it = paf_objects.begin();
	/*it = std::unique (paf_objects.begin(), paf_objects.end(), paf_unique);
	paf_objects.resize(std::distance(paf_objects.begin(), it));
	paf_objects.erase(std::remove_if(paf_objects.begin(), paf_objects.end(), [](std::unique_ptr<InputPafFile> &p){return p->seq_end_pos - p->seq_start_pos < 0.9 * p->seq_length;}), paf_objects.end());
	/*auto paf_cmp = [](const std::unique_ptr<InputPafFile>& a, const std::unique_ptr<InputPafFile>& b) {
		if (a->gen_name == b->gen_name) {
	 		if (a->gen_start_pos == b->gen_end_pos) {
	        	return (a->gen_end_pos > b->gen_end_pos);
	        }
	        return (a->gen_start_pos < b->gen_end_pos);
	    }
	    return (a->gen_name > b->gen_name);
	};*/
	std::sort(paf_objects.begin(), paf_objects.end(), paf_cmp);
	it = paf_objects.begin();
	int s, e;
	std::string n;
	while (it != --paf_objects.end()) {
		s = (*it)->gen_start_pos;
		e = (*it)->gen_end_pos;
		n = (*it)->gen_name;
		paf_objects.erase(std::remove_if(it+1, paf_objects.end(), [&s, &e, &n](std::unique_ptr<InputPafFile> &p){return p->gen_start_pos >= s && p->gen_end_pos <= e && p->gen_name == n
		&& (p->seq_end_pos - p->seq_start_pos > 0.9 * p->seq_length);}), paf_objects.end());
		if (it != --paf_objects.end()) {
			it++;
		}
	}
}


int main(int argc, char* argv[]) {
    std::string paf_file(argv[optind]);
    std::string fasta_file(argv[optind+1]);

    std::vector<std::unique_ptr<InputPafFile>> paf_objects;
    auto paf_parser = bioparser::createParser<bioparser::PafParser, InputPafFile>(paf_file);
    auto status = paf_parser->parse_objects(paf_objects, -1);

    std::unordered_map<std::string, std::vector<mapping>> sequence_mapping_details;


    for (auto& i : paf_objects){
        //std::cout << i->seq_name << "\t" << i->seq_start_pos << "\t" << i->seq_end_pos << std::endl;
        std::string key(i->seq_name);
        mapping mapp(i->seq_start_pos, i->seq_end_pos, i->gen_start_pos, i->gen_end_pos, i->gen_length, i->gen_name);
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

    /*for (auto it = sequence_mapping_details.begin(); it != sequence_mapping_details.end();) {
        if ((it->second).size() == 1) {
            //std::cout << "tu sam " << (it->second).size()  << std::endl;
            it = sequence_mapping_details.erase(it);
        } else{
            it++;
        }
    }*/

    std::unordered_map<std::string, std::vector<mapping>> chimeric_reads;
    std::unordered_map<std::string, std::vector<mapping>> repeating_reads;
    std::unordered_set<std::string> chimers;
    std::unordered_set<std::string> repeatings;
    std::vector<mapping> all_repeatings;

    for (itr = sequence_mapping_details.begin(); itr != sequence_mapping_details.end(); itr++) {
        if ((itr->second).size() > 1) {
            std::sort((itr->second).begin(), (itr->second).end(), sort_function);
            uint32_t seq_start = (itr->second)[0].seq_start;
            uint32_t seq_end = (itr->second)[0].seq_end;
            uint32_t gen_start = (itr->second)[0].gen_start;
            uint32_t gen_end = (itr->second)[0].gen_end;
            for (int i = 1 ; i < (itr->second).size(); i++) {
                if (!(seq_start <= (itr->second)[i].seq_start && seq_end >= (itr->second)[i].seq_end)) {
                    chimers.emplace(itr->first);
                    int seq_gap = abs((int)(std::max(seq_start, (itr->second)[i].seq_start) - std::min(seq_end, (itr->second)[i].seq_end)));
                    int ref_gap = abs((int)(std::max(gen_start, (itr->second)[i].gen_start) - std::min(seq_end, (itr->second)[i].seq_end)));
                    if(std::min(gen_start, (itr->second)[i].gen_start) > 100 || std::max(gen_end, (itr->second)[i].gen_end) < ((itr->second)[i].gen_length - 100)
                    && (std::min(ref_gap, seq_gap)/std::max(ref_gap, seq_gap)) > 0.12){
                        chimeric_reads[itr->first] = itr->second;
                        break;
                    }
                }
            }
            if (chimers.find(itr->first) == chimers.end()) {
                int i = 1;
                if ((itr->second)[0].seq_start == (itr->second)[1].seq_start && (itr->second)[0].seq_end == (itr->second)[1].seq_end) {
                    i = 0;
                }
                all_repeatings.insert(all_repeatings.end(), (itr->second).begin() + i, (itr->second).end());
                repeating_reads[itr->first] = itr->second;
            }
        }
    }


    /*for (auto &i : paf_objects) {
        std::cout << i->seq_name << "\t" << i->gen_name << "\t" << i->gen_start_pos << "\t" << i->gen_end_pos << std::endl;
    }*/

    std::sort(all_repeatings.begin(), all_repeatings.end(), repeats_by_length);

    std::unordered_map<std::string, std::vector<std::tuple<uint32_t, uint32_t>>> reference_repeating_regions = find_repeating_regions(all_repeatings);


    for (auto rep_it = repeating_reads.begin(); rep_it != repeating_reads.end(); rep_it++) {
        bool contained = false;
        int n = 0;
        //std::cout << "Sekvenca " << rep_it->first << " s mapiranjem " << longest_mapping.gen_start << "\t" << longest_mapping.gen_end << std::endl;
        do {
            mapping longest_mapping = (repeating_reads[rep_it->first])[n];
            for (auto regions_it = reference_repeating_regions.begin(); regions_it != reference_repeating_regions.end(); regions_it++) {
                if (longest_mapping.ref_name == regions_it->first) {
                    std::vector<std::tuple<uint32_t, uint32_t>> rep_regions = reference_repeating_regions[regions_it->first];
                    for (int j = 0; j < rep_regions.size(); j++) {
                        if (longest_mapping.gen_start >= std::get<0>(rep_regions[j]) && longest_mapping.gen_end <= std::get<1>(rep_regions[j])) {
                            std::cout << "Sekvenca " << rep_it->first << " s mapiranjem " << longest_mapping.gen_start << "\t" << longest_mapping.gen_end << std::endl;
                            contained = true;
                            break;
                        }
                    }
                }
                if (contained) {
                    break;
                }
            }
            if (!contained) {
                repeatings.emplace(rep_it->first);
            }
            n++;
        } while ((repeating_reads[rep_it->first])[n].seq_start == (repeating_reads[rep_it->first])[n-1].seq_start && (repeating_reads[rep_it->first])[n].seq_end == (repeating_reads[rep_it->first])[n-1].seq_end && !contained);
    }


    /*for(int i = 0; i < all_repeatings.size(); i++) {
        for (int j = 0; j < all_repeatings.size(); j++) {
            if (all_repeatings[i].ref_name == all_repeatings[j].ref_name) {
                if (abs((int)(all_repeatings[i].gen_start - all_repeatings[j].gen_start)) + abs((int)(all_repeatings[i].gen_end - all_repeatings[j].gen_end)) < 500) {
                    if (reference_repeating_regions.find(all_repeatings[i].ref_name) == reference_repeating_regions.end()){
                        std::vector<std::tuple<uint32_t, uint32_t>> vec;
                        vec.push_back(std::make_tuple(all_repeatings[i].gen_start, all_repeatings[i].gen_end));
                        reference_repeating_regions[all_repeatings[i].ref_name] = vec;
                    } else {
                        (reference_repeating_regions[all_repeatings[i].ref_name]).push_back(std::make_tuple(all_repeatings[i].gen_start, all_repeatings[i].gen_end));
                    }
                }
            }
        }
    }*/

    /*for (auto it = reference_repeating_regions.begin(); it != reference_repeating_regions.end(); it++) {
        std::sort((it->second).begin(), (it->second).end(), sort_reference_repeatings);
        std::vector<std::tuple<uint32_t, uint32_t>>::iterator vec_it;
        vec_it = std::unique((it->second).begin(), (it->second).end(), unique_reference_repeatings);
        (it->second).resize(std::distance((it->second).begin(), vec_it));
        /*std::cout << "Repetativne regije na referenci " << it->first << std::endl;
        for (int i = 0; i < (it->second).size(); i++) {
            std::cout << std::get<0>((it->second)[i]) << "\t" << std::get<1>((it->second)[i]) << std::endl;
        }*/
    //}

    /*for (itr = sequence_mapping_details.begin(); itr != sequence_mapping_details.end(); itr++) {
        std::cout << "Sekvenca " << itr->first << " ima ova mapiranja:" << std::endl;
        for (int i = 0; i < (itr->second).size(); i++) {
            std::cout << (itr->second)[i].seq_start << "\t" << (itr->second)[i].seq_end << "\t" << (itr->second)[i].gen_start << "\t" << (itr->second)[i].gen_end << std::endl;
        }
    }*/

    /*std::cout << "----------------------------------------------------------" << std::endl;

    for (itr = chimeric_reads.begin(); itr != chimeric_reads.end(); itr++) {
        std::cout << "Sekvenca " << itr->first << " ima ova mapiranja:" << std::endl;
        for (int i = 0; i < (itr->second).size(); i++) {
            std::cout << (itr->second)[i].seq_start << "\t" << (itr->second)[i].seq_end << "\t" << (itr->second)[i].gen_start << "\t" << (itr->second)[i].gen_end << "\t" << (itr->second)[i].ref_name << std::endl;
        }
    }

    std::cout << "_________________________________________________________" << std::endl;

    for (itr = repeating_reads.begin(); itr != repeating_reads.end(); itr++) {
        std::cout << "Sekvenca " << itr->first << " ima ova mapiranja:" << std::endl;
        for (int i = 0; i < (itr->second).size(); i++) {
            std::cout << (itr->second)[i].seq_start << "\t" << (itr->second)[i].seq_end << "\t" << (itr->second)[i].gen_start << "\t" << (itr->second)[i].gen_end << "\t" << (itr->second)[i].ref_name << std::endl;
        }
    }*/

    std::cout << "Kimernih sekvenci ima " << chimeric_reads.size() << std::endl;
    std::cout << "Repetativnih sekvenci ima " << repeating_reads.size() << std::endl;

    std::vector<uint64_t> chimeric_annotations;

    for(itr = chimeric_reads.begin(); itr != chimeric_reads.end(); itr++) {
        for(int i = 0; i < (itr->second).size() - 1; i++) {
            //std::cout << "na prvom " << (itr->second)[i].seq_start << " i " << (itr->second)[i].seq_end << ", na drugom " << (itr->second)[i+1].seq_start  << " i " << (itr->second)[i+1].seq_end << std::endl;
            if((itr->second)[i].seq_start != (itr->second)[i+1].seq_start && (itr->second)[i].seq_end != (itr->second)[i+1].seq_end){
                uint64_t positions = (itr->second)[i].seq_end;
                positions = positions << 32;
                positions = positions + (itr->second)[i+1].seq_start;
                chimeric_annotations.push_back(positions);

            }
        }
    }

    /*for(int i = 0; i < chimeric_annotations.size(); i++) {
        std::cout << chimeric_annotations[i] << std::endl;
    }*/

    std::vector<std::unique_ptr<InputFile>> first_object;
    if (fasta_file.find("fastq")>fasta_file.length() && fasta_file.find("fq")>fasta_file.length()){

        auto fasta_parser = bioparser::createParser<bioparser::FastaParser, InputFile>(fasta_file);
        fasta_parser->parse_objects(first_object, -1);

        std::cout << "We've parsed file." << std::endl;

    } else {
        auto fastq_parser = bioparser::createParser<bioparser::FastqParser, InputFile>(fasta_file);
        uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
            while (true) {
                auto status = fastq_parser->parse_objects(first_object, size_in_bytes);

                if (status == false) {
                    break;
                }
            }
        std::cout << "We've parsed file." << std::endl;
    }

    /*std::ofstream repeating_regions;
    repeating_regions.open("reference_repeating_regions.rpt");

    for (auto it = reference_repeating_regions.begin(); it != reference_repeating_regions.end(); it++) {
        std::sort((it->second).begin(), (it->second).end(), sort_reference_repeatings);
        std::vector<std::tuple<uint32_t, uint32_t>>::iterator vec_it;
        vec_it = std::unique((it->second).begin(), (it->second).end());
        (it->second).resize(std::distance((it->second).begin(), vec_it));
        for (int i = 0; i < (it->second).size(); i++) {
            repeating_regions << ">" << (it->first) << " :" << std::get<0>((it->second)[i]) << "-" << std::get<1>((it->second)[i]) << "\n";
        }
    }*/

    std::ofstream cleaned_sequences;
    cleaned_sequences.open("cleaned_sequences.fasta");

    std::ofstream chimeric;
    chimeric.open("chimeric_reads.fasta");

    std::ofstream repeating;
    repeating.open("repeating_reads.fasta");

     for(auto& i : first_object) {
        if(chimers.find(i->name) != chimers.end()){
            chimeric << ">" << i->name << "\t";
            std::vector<mapping> info = chimeric_reads[i->name];
            std::sort(info.begin(), info.end(), sort_function);
            std::vector<mapping>::iterator vec_it;
            vec_it = std::unique(info.begin(), info.end(), unique_read_function);
            //info.resize(std::distance(info.begin(), vec_it));
            std::set<std::tuple<uint32_t, uint32_t>> written;
            for (auto& i : info) {
                if (written.find(std::make_tuple(i.seq_start, i.seq_end)) == written.end()) {
                    chimeric << i.seq_start << " " << i.seq_end << "\t";
                    written.emplace(std::make_tuple(i.seq_start, i.seq_end));
                }
            }
            chimeric << "\n";
            chimeric << i->sequence << "\n";
        } else if (repeatings.find(i->name) != repeatings.end()) {
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
        } else {
            cleaned_sequences << ">" << i->name << "\n";
            cleaned_sequences << i->sequence << "\n";
        }
    }

    return 0;

}
