#include <cstdint>
#include <memory>
#include <vector>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#include "map.hpp"
#include "radixsort.hpp"
#include "revcomp.hpp"
#include "region.hpp"
#include "brown_minimizers.hpp"
#include "sam.hpp"

// Clip end of sequence that contains multiple Ns
// Args: seq - sequence to be clipped
//       k   - k-mer length
//       w   - window length
// Return: size of clipping; if negative, sequence should be rejected
int32_t clip(const std::string& seq, const uint32_t k, const uint32_t w) {
  std::size_t found = seq.find("NNN", 10);
  if (found == std::string::npos) return 0;
  if (found < (w + k - 1)) return -1;
  return seq.size() - found;
}

// Process (as) single-end sequencing read
// Args: sam          - SAM format output string
//       ref_index    - map of minimizer values to positions in target minimizers list
//       t_minimizers - list of target minimizers
//       reference    - FastAQ representation of reference
//       read         - FastAQ representation of read
//       parameters   - mapping parameters
// Return: none
void process_single(std::string& sam, const std::unordered_map<uint64_t, index_pos_t>& ref_index, 
                    const std::vector<minimizer_t>& t_minimizers, const std::unique_ptr<fastaq::FastAQ>& reference, 
                    const std::unique_ptr<fastaq::FastAQ>& read, const mapping_params& parameters) {
  int32_t clipped = clip(read->sequence, parameters.k, parameters.w);
  if (clipped < 0) {
    sam += unmapped_sam(read->name, read->sequence, read->quality, 0, 0, 0);
    return;
  }
  std::vector<minimizer_t> q_minimizers = brown::minimizers(read->sequence.c_str(), read->sequence.size() - clipped,
                                                            parameters.k, parameters.w);
  split_hits_t hits;
  find_minimizer_hits(hits.first, hits.second, ref_index, t_minimizers, q_minimizers);
  radixsort(hits.first);
  radixsort(hits.second);
  std::pair<bin_t, bin_t> candidates(extract_candidates(hits.first, parameters.threshold, read->sequence.size()),
                                     extract_candidates(hits.second, parameters.threshold, read->sequence.size()));
  if (candidates.first.size() == 0 && candidates.second.size() == 0) {
    sam += unmapped_sam(read->name, read->sequence, read->quality, 0, 0, 0);
    return;
  }
  std::unordered_set<uint32_t> processed;
  std::vector<std::pair<uint32_t, std::string>> mappings;
  for (auto& bin : candidates.first) {
    region_t reg = find_region(bin.second);
    expand_region(reg, read->sequence.size(), parameters.k, reference->sequence.size() - 1);
    if (processed.find(std::get<1>(reg.first)) != processed.end()) continue;
    processed.insert(std::get<1>(reg.first));
    mappings.emplace_back(sam_format_single(read->name, read->sequence, read->quality,
                                            reference->name, reference->sequence, reg, 
                                            parameters, clipped));
  }
  for (auto& bin : candidates.second) {
    region_t reg = find_region(bin.second);
    expand_region(reg, read->sequence.size(), parameters.k, reference->sequence.size() - 1);
    if (processed.find(std::get<1>(reg.first)) != processed.end()) continue;
    processed.insert(std::get<1>(reg.first));
    std::string rc = reverse_complement(read->sequence, 0, read->sequence.size());
    std::string rq = std::string(read->quality.rbegin(), read->quality.rend());
    mappings.emplace_back(sam_format_single(read->name, rc, rq,
                                            reference->name, reference->sequence, reg, 
                                            parameters, clipped));
  }
  std::stable_sort(mappings.begin(), mappings.end(), 
                   [] (const std::pair<uint32_t, std::string>& a, const std::pair<uint32_t, std::string>& b) {
                     return a.first > b.first;
                   }
  );
  uint32_t stop = parameters.all ? mappings.size() : 1;
  for (uint32_t i = 0; i < stop; ++i) {
    if (i) {
      std::size_t f_start = mappings[i].second.find('\t', 0) + 1;
      std::size_t f_end = mappings[i].second.find('\t', f_start);
      mappings[i].second.replace(f_start, f_end - f_start, 
          std::to_string(std::stoi(mappings[i].second.substr(f_start, f_end - f_start)) | 0x100));
      std::size_t pos = 0;
      for (uint32_t j = 0; j < 9; ++j) pos = mappings[i].second.find('\t', pos) + 1;
      std::size_t end = mappings[i].second.find('\t', pos);
      mappings[i].second.replace(pos, end - pos, "*");
      pos += 2;
      end = mappings[i].second.find('\t', pos);
      mappings[i].second.replace(pos, end - pos, "*");
    }
    sam += mappings[i].second;
  }
}

// Process paired-end sequencing reads
// Args: mappings   - list of SAM format output string with priority values
//       checked    - pair of corresponding lists containing lists of hits that were checked for insert size
//       reference  - FastAQ representation of reference
//       first      - FastAQ representation of first read from pair
//       second     - FastAQ representation of second read from pair
//       parameters - mapping parameters
// Return: none
void process_pairs(std::vector<std::pair<uint32_t, std::string>>& mappings, paired_checked_t& checked,
                   const std::unique_ptr<fastaq::FastAQ>& reference, 
                   const std::unique_ptr<fastaq::FastAQ>& first, const std::unique_ptr<fastaq::FastAQ>& second,
                   const mapping_params& parameters, const int32_t clipped1, const int32_t clipped2) {
  std::unordered_set<uint32_t> processed;
  for (uint32_t j = 0; j < checked.first.size(); ++j) {
    std::pair<region_t, region_t> region_pair(find_region(checked.first[j]),
                                              find_region(checked.second[j]));
    expand_region(region_pair.first, first->sequence.size(), 
                  parameters.k, reference->sequence.size() - 1);
    expand_region(region_pair.second, second->sequence.size(), 
                  parameters.k, reference->sequence.size() - 1);
    if (processed.find(std::get<1>(region_pair.first.first)) != processed.end()) continue;
    processed.insert(std::get<1>(region_pair.first.first));
    std::string rc;
    std::string rq;
    std::string* query1 = std::get<2>(region_pair.first.first)
                          ? &(rc = reverse_complement(first->sequence, 0, first->sequence.size()))
                          : &(first->sequence);
    std::string* quality1 = std::get<2>(region_pair.first.first)
                            ? &(rq = std::string(first->quality.rbegin(), first->quality.rend()))
                            : &(first->quality);
    std::string* query2 = std::get<2>(region_pair.second.first)
                          ? &(rc = reverse_complement(second->sequence, 0, second->sequence.size()))
                          : &(second->sequence);
    std::string* quality2 = std::get<2>(region_pair.second.first)
                            ? &(rq = std::string(second->quality.rbegin(), second->quality.rend()))
                            : &(second->quality);
    std::pair<uint32_t, std::pair<std::string, std::string>> sam_pair = sam_format_pair(
        first->name,
        *query1, *quality1, 
        *query2, *quality2, 
        reference->name, reference->sequence, region_pair, parameters, clipped1, clipped2);

    mappings.emplace_back(sam_pair.first, sam_pair.second.first);
    mappings.emplace_back(sam_pair.first, sam_pair.second.second);
  }
}

// Map single reads
// Args: ref_index    - map of minimizer values to positions in target minimizers list
//       t_minimizers - list of target minimizers
//       reference    - FastAQ representation of reference
//       reads        - list of FastAQ representations of reads
//       parameters   - mapping parameters
//       t_start      - starting index of reads to be mapped       // used for parallelization
//       t_end        - last index of reads to be mapped           // [t_start, t_end)
// Return: mapping result in SAM format
std::string map_single(const std::unordered_map<uint64_t, index_pos_t>& ref_index, 
                       const std::vector<minimizer_t>& t_minimizers,
                       const std::unique_ptr<fastaq::FastAQ>& reference, 
                       const std::vector<std::unique_ptr<fastaq::FastAQ>>& reads, const mapping_params& parameters,
                       uint32_t t_start, uint32_t t_end) {
  std::string sam;
  for (uint32_t i = t_start; i < t_end; ++i) {
    process_single(sam, ref_index, t_minimizers, reference, reads[i], parameters);
  }
  return sam;
}

// Map paired-end reads as single reads
// Args: ref_index    - map of minimizer values to positions in target minimizers list
//       t_minimizers - list of target minimizers
//       reference    - FastAQ representation of reference
//       paired_reads - pair of lists of FastAQ representations of paired reads
//       parameters   - mapping parameters
//       t_start      - starting index of reads to be mapped       // used for parallelization as:
//       t_end        - last index of reads to be mapped           // [t_start, t_end)
// Return: mapping result in SAM format
std::string map_as_single(const std::unordered_map<uint64_t, index_pos_t>& ref_index, 
                          const std::vector<minimizer_t>& t_minimizers,
                          const std::unique_ptr<fastaq::FastAQ>& reference, const paired_reads_t& paired_reads,
                          const mapping_params& parameters, uint32_t t_start, uint32_t t_end) {
  std::string sam;
  for (uint32_t i = t_start; i < t_end; ++i) {
    process_single(sam, ref_index, t_minimizers, reference, paired_reads.first[i], parameters);
    process_single(sam, ref_index, t_minimizers, reference, paired_reads.second[i], parameters);
  }
  return sam;
}

// Map paired-end reads as single reads
// Args: ref_index    - map of minimizer values to positions in target minimizers list
//       t_minimizers - list of target minimizers
//       reference    - FastAQ representation of reference
//       paired_reads - pair of lists of FastAQ representations of paired reads
//       parameters   - mapping parameters
//       t_start      - starting index of reads to be mapped       // used for parallelization as:
//       t_end        - last index of reads to be mapped           // [t_start, t_end)
// Return: mapping result in SAM format
std::string map_paired(const std::unordered_map<uint64_t, index_pos_t>& ref_index, 
                       const std::vector<minimizer_t>& t_minimizers,
                       const std::unique_ptr<fastaq::FastAQ>& reference, const paired_reads_t& paired_reads,
                       const mapping_params& parameters, uint32_t t_start, uint32_t t_end) {
  std::string sam;
  for (uint32_t i = t_start; i < t_end; ++i) {
    int32_t clipped1 = clip(paired_reads.first[i]->sequence, parameters.k, parameters.w);
    int32_t clipped2 = clip(paired_reads.second[i]->sequence, parameters.k, parameters.w);
    if (clipped1 < 0 || clipped2 < 0) {
      sam += unmapped_sam(paired_reads.first[i]->name, paired_reads.first[i]->sequence,
                          paired_reads.first[i]->quality, 1, 1, 0)
             + unmapped_sam(paired_reads.second[i]->name, paired_reads.second[i]->sequence, 
                            paired_reads.second[i]->quality, 1, 0, 1);
      continue;
    }
    paired_minimizers_t p_minimizers(brown::minimizers(paired_reads.first[i]->sequence.c_str(),
                                                       paired_reads.first[i]->sequence.size() - clipped1,
                                                       parameters.k, parameters.w),
                                     brown::minimizers(paired_reads.second[i]->sequence.c_str(),
                                                       paired_reads.second[i]->sequence.size() - clipped2,
                                                       parameters.k, parameters.w));
    split_hits_t hits1;
    split_hits_t hits2;
    find_minimizer_hits(hits1.first, hits1.second, ref_index, t_minimizers, p_minimizers.first);
    find_minimizer_hits(hits2.first, hits2.second, ref_index, t_minimizers, p_minimizers.second);

    radixsort(hits1.first);
    radixsort(hits1.second);
    radixsort(hits2.first);
    radixsort(hits2.second);

    std::pair<bin_t, bin_t> candidates1(
        extract_candidates(hits1.first, parameters.threshold, paired_reads.first[i]->sequence.size()),
        extract_candidates(hits2.second, parameters.threshold, paired_reads.second[i]->sequence.size()));
    std::pair<bin_t, bin_t> candidates2(
        extract_candidates(hits1.second, parameters.threshold, paired_reads.first[i]->sequence.size()),
        extract_candidates(hits2.first, parameters.threshold, paired_reads.second[i]->sequence.size()));
    paired_checked_t checked1 = check_pairing(candidates1, parameters.insert_size,
                                              paired_reads.second[i]->sequence.size());
    paired_checked_t checked2 = check_pairing(candidates2, parameters.insert_size,
                                              paired_reads.first[i]->sequence.size());

    std::vector<std::pair<uint32_t, std::string>> mappings;
    process_pairs(mappings, checked1, reference, paired_reads.first[i], paired_reads.second[i], parameters, clipped1, clipped2);
    process_pairs(mappings, checked2, reference, paired_reads.first[i], paired_reads.second[i], parameters, clipped1, clipped2);
    if (mappings.size() == 0) {
      sam += unmapped_sam(paired_reads.first[i]->name, paired_reads.first[i]->sequence,
                          paired_reads.first[i]->quality, 1, 1, 0)
             + unmapped_sam(paired_reads.second[i]->name, paired_reads.second[i]->sequence, 
                            paired_reads.second[i]->quality, 1, 0, 1);
      continue;
    }
    std::stable_sort(mappings.begin(), mappings.end(), 
                     [] (const std::pair<uint32_t, std::string>& a, const std::pair<uint32_t, std::string>& b) {
                       return a.first > b.first;
                     }
    );
    std::string sam1;
    std::string sam2;
    uint32_t stop = parameters.all ? mappings.size() - 1 : 1;
    for (uint32_t j = 0; j < stop; j+=2) {
      if (j) {
        std::size_t f_start = mappings[j].second.find('\t', 0) + 1;
        std::size_t f_end = mappings[j].second.find('\t', f_start);
        mappings[j].second.replace(f_start, f_end - f_start, 
            std::to_string(std::stoi(mappings[j].second.substr(f_start, f_end - f_start)) | 0x100));
        f_start = mappings[j + 1].second.find('\t', 0) + 1;
        f_end = mappings[j + 1].second.find('\t', f_start);
        mappings[j + 1].second.replace(f_start, f_end - f_start, 
            std::to_string(std::stoi(mappings[j + 1].second.substr(f_start, f_end - f_start)) | 0x100));
        std::size_t pos = 0;
        for (uint32_t k = 0; k < 9; ++k) pos = mappings[j].second.find('\t', pos) + 1;
        std::size_t end = mappings[j].second.find('\t', pos);
        mappings[j].second.replace(pos, end - pos, "*");
        pos += 2;
        end = mappings[j].second.find('\t', pos);
        mappings[j].second.replace(pos, end - pos, "*");
        pos = 0;
        for (uint32_t k = 0; k < 9; ++k) pos = mappings[j + 1].second.find('\t', pos) + 1;
        end = mappings[j + 1].second.find('\t', pos);
        mappings[j + 1].second.replace(pos, end - pos, "*");
        pos += 2;
        end = mappings[j + 1].second.find('\t', pos);
        mappings[j + 1].second.replace(pos, end - pos, "*");
      }
      sam1 += mappings[j].second;
      sam2 += mappings[j + 1].second;
    }
    sam += sam1 + sam2;
  }
  return sam;
}