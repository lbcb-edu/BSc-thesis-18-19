#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <cctype>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <set>
#include <queue>
#include <deque>
#include <unordered_map>
#include <tuple>
#include <utility>
#include <algorithm>

#include "srmapper.hpp"
#include "fastaq.hpp"
#include "brown_minimizers.hpp"
#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"

// minimizer: value, position, origin
typedef std::tuple<uint64_t, uint32_t, bool> minimizer_t;
// index: position, range
typedef std::pair<uint32_t, uint32_t> index_pos_t;
// hit: query position, reference position, relative strand
typedef std::tuple<uint32_t, uint32_t, bool> minimizer_hit_t;
// region = two "hits" with swapped positions for reverse complement
typedef std::pair<minimizer_hit_t, minimizer_hit_t> region_t;
// paired read
typedef std::pair<std::vector<std::unique_ptr<fastaq::FastAQ>>, std::vector<std::unique_ptr<fastaq::FastAQ>>> paired_reads_t;
// paired read minimizers
typedef std::pair<std::vector<minimizer_t>, std::vector<minimizer_t>> paired_minimizers_t;
// paired read hits
typedef std::pair<std::vector<minimizer_hit_t>, std::vector<minimizer_hit_t>> split_hits_t;
// bin
typedef std::unordered_map<uint32_t, std::vector<minimizer_hit_t>> bin_t;
// paired checked and grouped hits (candidates)
typedef std::pair<std::vector<std::vector<minimizer_hit_t>>, std::vector<std::vector<minimizer_hit_t>>> paired_checked_t;

const std::set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
const std::set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};

static struct option long_options[] = {
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},
  {"paired", no_argument, NULL, 'p'},
  {"window_length", required_argument, NULL, 'w'},
  {"kmers", required_argument, NULL, 'k'},
  {"frequency", required_argument, NULL, 'f'},
  {"insert_size", required_argument, NULL, 'i'},
  {"region_size", required_argument, NULL, 'r'},
  {"threads", required_argument, NULL, 't'},
  {"threshold", required_argument, NULL, 'T'},
  {NULL, no_argument, NULL, 0}
};

typedef struct {
  uint32_t k;
  uint32_t w;
  float f;
  uint32_t insert_size;
  uint32_t region_size;
  uint32_t threshold;
} mapping_params;

void help(void) {
  printf("srmapper - tool for mapping short reads to reference genome.\n\n"

         "Usage: srmapper [OPTIONS] reference [reads]\n"
         "  reference - FASTA file containing reference genome\n"
         "  reads     - one or two FASTA/FASTQ file containing a set of fragments\n"
         "              one = single-read sequencing reads\n"
         "              two = paired-end read sequencing reads\n\n"

         "Supported file extensions: .fasta\n"
         "                           .fa\n"
         "                           .fastq\n"
         "                           .fq\n"
         "                           .fasta.gz\n"
         "                           .fa.gz\n"
         "                           .fastq.gz\n"
         "                           .fq.gz\n\n"
         
         "OPTIONS:\n"
         "  -h  or  --help           print help (displayed now) and exit\n"
         "  -v  or  --version        print version info and exit\n"
         "  -p  or  --paired         use pairing information (needs insert size information)"
         "  -k  or  --kmers          <uint>\n"
         "                             default: 15\n"
         "                             constraints: largest supported is 16\n"
         "                             number of letters in substrings\n"
         "  -w  or  --window_length  <uint>\n"
         "                             default: 5\n"
         "                             length of window\n"
         "  -f  or  --frequency      <float>\n"
         "                             default: 0.001\n"
         "                             constraints: must be from [0, 1]\n"
         "                             number of most frequent minimizers that\n"
         "                             are not taken into account\n"
         "  -i  or  --insert_size    <uint>\n"
         "                             default: 215\n"
         "                             fragment insert size\n"
         "  -r  or  --region_size    <uint>\n"
         "                             default: 100\n"
         "                             region size to count hits from\n"
         "  -t  or  --threads        <uint>\n"
         "                             default: 3\n"
         "                             number of threads\n"
         "  -T  or  --threshold      <uint>\n"
         "                             default: 1\n"
         "                             number of hits needed in order to consider\n"
         "                             a region a candidate for mapping\n"
  );
}

void version(void) {
  printf("srmapper %d.%d\n",
    srmapper_VERSION_MAJOR,
    srmapper_VERSION_MINOR
  );
}

bool check_extension(const std::string& filename, const std::set<std::string>& extensions) {
  for (const auto& it : extensions) {
    if (filename.size() > it.size()) {
      if (filename.compare(filename.size()-it.size(), std::string::npos, it) == 0) {
        return true;
      }
    }
  }
  return false;
}

bool hit_ordering(const minimizer_hit_t& a, const minimizer_hit_t& b) {
  if (std::get<2>(a) == std::get<2>(b)) {
    if(std::get<1>(a) == std::get<1>(b)) {
      return std::get<0>(a) < std::get<0>(b);
    }
    return std::get<1>(a) < std::get<1>(b);
  }
  return std::get<2>(a) < std::get<2>(b);
}

void radixsort(std::vector<minimizer_hit_t>& hits) {
  if (hits.size() <= 1) return;
  uint32_t max_val = std::get<1>(hits[0]);
  const uint32_t size = hits.size();
  for (uint32_t i = 1; i < size; ++i) {
    if (std::get<1>(hits[i]) > max_val) {
      max_val = std::get<1>(hits[i]);
    }
  }
  std::vector<minimizer_hit_t> out(size);
  for (uint32_t e = 1; max_val / e > 0; e *= 10) {
    uint32_t count[10];
    for (uint32_t i = 0; i < 10; ++i) {
      count[i] = 0;
    }
    for (uint32_t i = 0; i < size; ++i) {
      count[(std::get<1>(hits[i]) / e) % 10]++;
    }
    for (uint32_t i = 1; i < 10; ++i) {
      count[i] += count[i - 1];
    }
    for (uint32_t i = size - 1; i != (uint32_t)(-1); --i) {
      out[count[(std::get<1>(hits[i]) / e) % 10] - 1] = hits[i];
      count[(std::get<1>(hits[i]) / e) % 10]--;
    }
    hits.assign(out.begin(), out.end());
  }
}

void prep_ref(std::vector<minimizer_t>& t_minimizers, const float f) {
  std::unordered_map<uint64_t, uint32_t> ref_min_frequency;
  for (const auto& minimizer : t_minimizers) {
    ref_min_frequency[std::get<0>(minimizer)]++;
  }

  std::vector<uint32_t> occurences;
  occurences.reserve(ref_min_frequency.size());
  for (const auto& entry : ref_min_frequency) {
    occurences.push_back(entry.second);
  }

  uint32_t position = (uint32_t)((1.0f - f) * (occurences.size() - 1.0f));
  std::sort(occurences.begin(), occurences.end());

  uint32_t cutoff_freq = occurences[position] == 1 ? 2 : occurences[position];

  std::vector<minimizer_t> temp;
  temp.reserve(t_minimizers.size());
  for (const auto& minimizer : t_minimizers) {
    if (ref_min_frequency[std::get<0>(minimizer)] < cutoff_freq) {
      temp.push_back(minimizer);
    }
  }
  std::swap(t_minimizers, temp);
  std::sort(t_minimizers.begin(), t_minimizers.end(),
      [] (const minimizer_t& a, const minimizer_t& b) {
        return (std::get<0>(a) < std::get<0>(b));
      });
  t_minimizers.shrink_to_fit();
}

std::unordered_map<uint64_t, index_pos_t> index_ref(const std::vector<minimizer_t>& t_minimizers) {
  std::unordered_map<uint64_t, index_pos_t> ref_index;
  uint32_t pos = 0;
  uint32_t num = 0;
  uint32_t prev_min = std::get<0>(t_minimizers[0]);
  for (const auto& minimizer : t_minimizers) {
    if (prev_min != std::get<0>(minimizer)) {
      ref_index[prev_min] = std::make_pair(pos, num);
      pos += num;
      num = 1;
      prev_min = std::get<0>(minimizer); 
    } else {
      num++;
    }
  }
  ref_index[prev_min] = std::make_pair(pos, num);
  return ref_index;
}

std::vector<uint32_t> find_ref_region_hits(
    const std::unordered_map<uint64_t, index_pos_t>& ref_index,
    const std::vector<minimizer_t>& t_minimizers,
    const std::vector<minimizer_t>& q_minimizers,
    const uint32_t region_size) {

  std::vector<uint32_t> hits;

  for (const auto& minimizer : q_minimizers) {
    auto found = ref_index.find(std::get<0>(minimizer));
    if (found != ref_index.end()) {
      for (uint32_t i = 0; i < found->second.second; i++) {
        hits.push_back(std::get<1>(t_minimizers[found->second.first + i]) / region_size);  
      }
    }
  }

  return hits;
}

void find_minimizer_hits(
    std::vector<minimizer_hit_t>& fwd,
    std::vector<minimizer_hit_t>& rev,
    const std::unordered_map<uint64_t, index_pos_t>& ref_index,
    const std::vector<minimizer_t>& t_minimizers,
    const std::vector<minimizer_t>& q_minimizers) {
  for (const auto& minimizer : q_minimizers) {
    auto found = ref_index.find(std::get<0>(minimizer));
    if (found != ref_index.end()) {
      for (uint32_t i = 0; i < found->second.second; i++) {
        if (std::get<2>(minimizer) == std::get<2>(t_minimizers[found->second.first + i])) {
          fwd.emplace_back(
              std::get<1>(minimizer),
              std::get<1>(t_minimizers[found->second.first + i]),
              0);  
        } else {
          rev.emplace_back(
              std::get<1>(minimizer),
              std::get<1>(t_minimizers[found->second.first + i]),
              1);
        }
      }
    }
  }
}

bin_t extract_candidates(
    const std::vector<minimizer_hit_t>& hits, 
    const uint32_t threshold, const uint32_t region_size) {
  bin_t candidates;
  std::vector<minimizer_hit_t> temp_c;
  std::vector<minimizer_hit_t> temp_n;
  uint32_t current = 0;
  uint32_t count = 0;
  uint32_t next = 0;
  uint32_t next_count = 0;
  for (uint32_t i = 0; i < hits.size(); ++i) {
    if (std::get<1>(hits[i]) / region_size == current) {
      count++;
      temp_c.push_back(hits[i]);
      continue;
    }
    if (std::get<1>(hits[i]) / region_size == current + 1) {
      count++;
      next_count++;
      next = current + 1;
      temp_c.push_back(hits[i]);
      temp_n.push_back(hits[i]);
      continue;
    }
    if (count >= threshold) {
      // NEEDS HELP
      candidates[current] = temp_c;
    }
    if (next == current + 1) {
      current = next;
      count = next_count;
      temp_c.assign(temp_n.begin(), temp_n.end());
      temp_n.clear();
    } else {
      current = std::get<1>(hits[i]) / region_size;
      count = 1;
    }
  }
  if (count >= threshold) {
    candidates[current] = temp_c;
  }
  return candidates;
}

paired_checked_t check_pairing(std::pair<bin_t, bin_t>& candidates,
                               const uint32_t insert_size,
                               const uint32_t region_size) {
  paired_checked_t checked;
  for (const auto& bin : candidates.first) {
    auto found = candidates.second.find(bin.first + ((insert_size - 100) / region_size));
    if (found != candidates.second.end()) {
      checked.first.emplace_back(bin.second.begin(), bin.second.end());
      checked.second.emplace_back(found->second.begin(), found->second.end());
      continue;
    }
    found = candidates.second.find(bin.first + ((insert_size - 100) / region_size) - 1);
    if (found != candidates.second.end()) {
      checked.first.emplace_back(bin.second.begin(), bin.second.end());
      checked.second.emplace_back(found->second.begin(), found->second.end());
      continue;
    }
    found = candidates.second.find(bin.first + ((insert_size - 100) / region_size) + 1);
    if (found != candidates.second.end()) {
      checked.first.emplace_back(bin.second.begin(), bin.second.end());
      checked.second.emplace_back(found->second.begin(), found->second.end());
      continue;
    }
    found = candidates.second.find(bin.first - ((insert_size - 100) / region_size));
    if (found != candidates.second.end()) {
      checked.first.emplace_back(bin.second.begin(), bin.second.end());
      checked.second.emplace_back(found->second.begin(), found->second.end());
      continue;
    }
    found = candidates.second.find(bin.first - ((insert_size - 100) / region_size) - 1);
    if (found != candidates.second.end()) {
      checked.first.emplace_back(bin.second.begin(), bin.second.end());
      checked.second.emplace_back(found->second.begin(), found->second.end());
      continue;
    }
    found = candidates.second.find(bin.first - ((insert_size - 100) / region_size) + 1);
    if (found != candidates.second.end()) {
      checked.first.emplace_back(bin.second.begin(), bin.second.end());
      checked.second.emplace_back(found->second.begin(), found->second.end());
      continue;
    }
  }
  return checked;
}

std::pair<region_t, unsigned int> find_region(std::vector<minimizer_hit_t>& hits, unsigned int k) {
  std::pair<region_t, unsigned int> region;
  region = std::make_pair(std::make_pair(std::make_tuple(0, 0, 0), std::make_tuple(0, 0, 0)), 0);
  if (hits.empty()) {
    return region;
  }

  std::sort(hits.begin(), hits.end(), 
            [] (const minimizer_hit_t& a, const minimizer_hit_t& b) {
              if (std::get<2>(a) == 0) {
                if (std::get<0>(a) == std::get<0>(b)) {
                  return std::get<1>(a) < std::get<1>(b);
                }
                return std::get<0>(a) < std::get<0>(b);
              } else {
                if (std::get<0>(a) == std::get<0>(b)) {
                  return std::get<1>(a) < std::get<1>(b);
                }
                return std::get<0>(a) > std::get<0>(b);
              }
            }
  );

  std::vector<region_t> lis(hits.size(),
      std::make_pair(std::make_tuple(0, 0, 0), std::make_tuple(0, 0, 0)));

  std::vector<unsigned int> matches(hits.size(), 0);

  unsigned int len = 1;
  lis[0] = std::make_pair(hits[0], hits[0]);
  matches[0] = k;

  for (unsigned int i = 1; i < hits.size(); ++i) {
    if (std::get<1>(hits[i]) > std::get<1>(lis[len - 1].second)) {
      unsigned int diff = std::get<0>(hits[i]) - std::get<0>(lis[len - 1].second);
      diff = std::min(diff, std::get<1>(hits[i]) - std::get<1>(lis[len - 1].second));
      diff = std::min(diff, k);
      matches[len] = matches[len - 1] + diff;
      lis[len] = lis[len - 1];
      lis[len].second = hits[i];
      len++;
    } else {
      auto pair_it = std::upper_bound(lis.begin(), lis.begin() + len, hits[i],
          [] (const minimizer_hit_t& a,
              const region_t& b) {
            return (std::get<1>(a) < std::get<1>(b.second));
          });
      if (pair_it == lis.begin()) {
        *pair_it = std::make_pair(hits[i], hits[i]);
      } else {
        *pair_it = std::make_pair((pair_it - 1)->first, hits[i]);
        unsigned int diff = std::get<0>(hits[i]) - std::get<0>((pair_it - 1)->second);
        diff = std::min(diff, std::get<1>(hits[i]) - std::get<1>((pair_it - 1)->second));
        diff = std::min(diff, k);
        matches[pair_it - lis.begin()] = matches[pair_it - 1 - lis.begin()] + diff;
      }
    }
  }
  region_t max_reg = lis[len - 1];
  if (std::get<2>(max_reg.first) == 1) {
    unsigned int temp = std::get<0>(max_reg.first);
    std::get<0>(max_reg.first) = std::get<0>(max_reg.second);
    std::get<0>(max_reg.second) = temp;
  }
  region = std::make_pair(max_reg, matches[len - 1]);
  return region; 
}

void expand_regions(std::vector<region_t>& regions, uint32_t read_size, uint32_t k) {
  for (auto& reg : regions) {
    if (std::get<2>(reg.first) == 0) {
      std::get<1>(reg.first) -= std::get<0>(reg.first);
      std::get<1>(reg.second) += read_size - std::get<0>(reg.second);
    } else {
      std::get<1>(reg.first) -= read_size - std::get<0>(reg.second) - k;
      std::get<1>(reg.second) += std::get<0>(reg.first) + k;
    }
    std::get<0>(reg.first) = 0;
    std::get<0>(reg.second) = read_size - 1;
  }
}

void map_single(const std::unordered_map<uint64_t, index_pos_t>& ref_index, const std::vector<minimizer_t>& t_minimizers,
                const std::vector<std::unique_ptr<fastaq::FastAQ>>& reference, 
                const std::vector<std::unique_ptr<fastaq::FastAQ>>& reads, const mapping_params& parameters) {
  for (uint32_t i = 0; i < reads.size(); ++i) {
    std::vector<minimizer_t> q_minimizers = brown::minimizers(reads[i]->sequence.c_str(), reads[i]->sequence.size(),
                                                             parameters.k, parameters.w);
    split_hits_t hits;
    find_minimizer_hits(hits.first, hits.second, ref_index, t_minimizers, q_minimizers);
    radixsort(hits.first);
    radixsort(hits.second);
    std::pair<bin_t, bin_t> candidates(extract_candidates(hits.first, parameters.threshold, parameters.region_size),
                                       extract_candidates(hits.second, parameters.threshold, parameters.region_size));
    std::vector<region_t> regions;
    for (auto& bin : candidates.first) {
      auto reg = find_region(bin.second, parameters.k);
      if (reg.second > 1 * parameters.k) {
        regions.push_back(reg.first);
      }
    }
    for (auto& bin : candidates.second) {
      auto reg = find_region(bin.second, parameters.k);
      if (reg.second > 1 * parameters.k) {
        regions.push_back(reg.first);
      }
    }
    expand_regions(regions, reads[i]->sequence.size(), parameters.k);
    for (const auto& region : regions) {
      std::cout << reads[i]->name << "\t" << std::get<2>(region.first) 
                << "\t" << std::get<0>(region.first) + 1 << "\t" << std::get<1>(region.first) + 1 
                << "\t" << std::get<0>(region.second) + 1 << "\t" << std::get<1>(region.second) + 1 
                << std::endl;
    }
  }
}

void map_paired(const std::unordered_map<uint64_t, index_pos_t>& ref_index, const std::vector<minimizer_t>& t_minimizers,
                const std::vector<std::unique_ptr<fastaq::FastAQ>>& reference, const paired_reads_t& paired_reads,
                const mapping_params& parameters) {
  for (uint32_t i = 0; i < paired_reads.first.size(); ++i) {
    paired_minimizers_t p_minimizers(brown::minimizers(paired_reads.first[i]->sequence.c_str(),
                                                       paired_reads.first[i]->sequence.size(),
                                                       parameters.k, parameters.w),
                                     brown::minimizers(paired_reads.second[i]->sequence.c_str(),
                                                       paired_reads.second[i]->sequence.size(),
                                                       parameters.k, parameters.w));
    split_hits_t hits1; // first 0, second 1
    split_hits_t hits2;
    find_minimizer_hits(hits1.first, hits1.second, ref_index, t_minimizers, p_minimizers.first);
    find_minimizer_hits(hits2.first, hits2.second, ref_index, t_minimizers, p_minimizers.second);

    radixsort(hits1.first);
    radixsort(hits1.second);
    radixsort(hits2.first);
    radixsort(hits2.second);

    std::pair<bin_t, bin_t> candidates1(
        extract_candidates(hits1.first, parameters.threshold, parameters.region_size),
        extract_candidates(hits2.second, parameters.threshold, parameters.region_size));
    std::pair<bin_t, bin_t> candidates2(
        extract_candidates(hits2.first, parameters.threshold, parameters.region_size),
        extract_candidates(hits1.second, parameters.threshold, parameters.region_size));
    paired_checked_t checked1 = check_pairing(candidates1, parameters.insert_size, 
                                              parameters.region_size);
    paired_checked_t checked2 = check_pairing(candidates2, parameters.insert_size, 
                                              parameters.region_size);
    std::vector<region_t> regions;
    for (auto& hits : checked1.first) {
      std::pair<region_t, unsigned int> reg = find_region(hits, parameters.k);
      if (reg.second > 1 * parameters.k) {
        regions.push_back(reg.first);
      }
    }
    for (auto& hits : checked2.second) {
      std::pair<region_t, unsigned int> reg = find_region(hits, parameters.k);
      if (reg.second > 1 * parameters.k) {
        regions.push_back(reg.first);
      }
    }
    for (auto& hits : checked1.second) {
      std::pair<region_t, unsigned int> reg = find_region(hits, parameters.k);
      if (reg.second > 1 * parameters.k) {
        regions.push_back(reg.first);
      }
    }
    for (auto& hits : checked2.first) {
      std::pair<region_t, unsigned int> reg = find_region(hits, parameters.k);
      if (reg.second > 1 * parameters.k) {
        regions.push_back(reg.first);
      }
    }
    expand_regions(regions, paired_reads.first[i]->sequence.size(), parameters.k);
    for (const auto& region : regions) {
      std::cout << paired_reads.first[i]->name << "\t" << std::get<2>(region.first) 
                << "\t" << std::get<0>(region.first) + 1 << "\t" << std::get<1>(region.first) + 1 
                << "\t" << std::get<0>(region.second) + 1 << "\t" << std::get<1>(region.second) + 1 
                << std::endl;
    }
  }
}

int main(int argc, char **argv) {
  int optchr;
  mapping_params parameters;
  parameters.k = 18;
  parameters.w = 3;
  parameters.f = 0.001f;
  parameters.insert_size = 215;
  parameters.region_size = 100;
  parameters.threshold = 1;
  bool paired = false;
  uint32_t threads = 3;

  while ((optchr = getopt_long(argc, argv, "hvk:w:f:i:r:t:T:p", long_options, NULL)) != -1) {
    switch (optchr) {
      case 'h': {
        help();
        exit(0);
      }
      case 'v': {
        version();
        exit(0);
      }
      case 'p': {
        paired = true;
        break;
      }
      case 'k': {
        parameters.k = atoi(optarg);
        break;
      }
      case 'w': {
        parameters.w = atoi(optarg);
        break;
      }
      case 'f': {
        parameters.f = atof(optarg);
        if (parameters.f < 0.0f || parameters.f > 1.0f) {
          fprintf(stderr, "[srmapper] error: f must be from [0, 1].\n"); 
          exit(1); 
        }
        break;
      }
      case 'i': {
        parameters.insert_size = atoi(optarg);
        break;
      }
      case 'r': {
        parameters.region_size = atoi(optarg);
        break;
      }
      case 't': {
        threads = atoi(optarg);
        break;
      }
      case 'T': {
        parameters.threshold = atoi(optarg);
        break;
      }
      default: {
        fprintf(stderr, "[srmapper] error: Unknown option. Type %s --help for usage.\n", argv[0]);
        exit(1);
      }
    }
  }

  if (argc - optind > 3 || argc - optind < 2) {
    fprintf(stderr, "[srmapper] error: Expected read(s) and reference. Use --help for usage.\n");
    exit(1);
  }
  if (argc - optind == 2 && paired) {
    fprintf(stderr, "[srmapper] error: Expected paired reads in order to use pairing information (option -p).\n");
    exit(1);
  }

  // fprintf(stderr, "\nMapping process started with parameters utilizing %u threads:\n"
  //                 "  k             = %u\n"
  //                 "  window length = %u\n"
  //                 "  top f freq    = %g\n"
  //                 "  insert size   = %u\n"
  //                 "  region size   = %u\n"
  //                 "  threshold     = %u\n",
  //                 threads,
  //                 parameters.k, parameters.w, parameters.f,
  //                 parameters.insert_size, parameters.region_size, parameters.threshold);

  fprintf(stderr, "\nLoading reference... ");

  std::string reference_file(argv[optind]);
  if (!check_extension(reference_file, fasta_formats)) {
      fprintf(stderr, "[srmapper] error: Unsupported reference file format. Check --help for supported file formats.\n");
      exit(1);
  }
  std::vector<std::unique_ptr<fastaq::FastAQ>> reference;
  fastaq::FastAQ::parse(reference, reference_file, check_extension(reference_file, fasta_formats));

  fprintf(stderr, "\rLoaded reference.        \n"
                  "\nIndexing reference... ");

  std::vector<minimizer_t> t_minimizers = brown::minimizers(reference[0]->sequence.c_str(), 
                                                            reference[0]->sequence.size(), 
                                                            parameters.k, parameters.w);
  prep_ref(t_minimizers, parameters.f);
  std::unordered_map<uint64_t, index_pos_t> ref_index = index_ref(t_minimizers);

  fprintf(stderr, "\rIndexed reference.        \n\n");
  
  fastaq::FastAQ::print_statistics(reference, reference_file);
  
  if (argc - optind == 3) {
    fprintf(stderr, "\nLoading paired-end reads... ");

    std::string reads_file1(argv[optind + 1]);
    std::string reads_file2(argv[optind + 2]);

    if (!(check_extension(reads_file1, fasta_formats) || check_extension(reads_file1, fastq_formats))
        || !(check_extension(reads_file2, fasta_formats) || check_extension(reads_file2, fastq_formats))) {
      fprintf(stderr, "[srmapper] error: Unsupported paired-end reads formats. Check --help for supported file formats.\n");
      exit(1);
    }
    paired_reads_t paired_reads;
    fastaq::FastAQ::parse(paired_reads.first, reads_file1, check_extension(reads_file1, fasta_formats));
    fastaq::FastAQ::parse(paired_reads.second, reads_file2, check_extension(reads_file2, fasta_formats));

    fprintf(stderr, "\rLoaded paired-end reads.        \n\n");

    fastaq::stats pr1_stats = fastaq::FastAQ::print_statistics(paired_reads.first, reads_file1);
    fastaq::stats pr2_stats = fastaq::FastAQ::print_statistics(paired_reads.second, reads_file2);
    fprintf(stderr, "\n");
    if (paired) {
      fprintf(stderr, "Using insert size information.\n");
      if (pr1_stats.num != pr2_stats.num) {
        fprintf(stderr, "[srmapper] error: Paired-end read files must have equal number of reads (pairs).\n");
        exit(1);
      }
      if (pr1_stats.max - pr1_stats.min > 0 || pr2_stats.max - pr2_stats.min > 0) {
        fprintf(stderr, "[srmapper] warning: Reads are not of fixed size.\n");
      }
      map_paired(ref_index, t_minimizers, reference, paired_reads, parameters);
    } else {

    }
  } else {
    fprintf(stderr, "\nLoading reads... ");
    std::string reads_file(argv[optind + 1]);

    if (!(check_extension(reads_file, fasta_formats) || check_extension(reads_file, fastq_formats))) {
      fprintf(stderr, "[srmapper] error: Unsupported format. Check --help for supported file formats.\n");
      exit(1);
    }

    std::vector<std::unique_ptr<fastaq::FastAQ>> reads;
    fastaq::FastAQ::parse(reads, reads_file, check_extension(reads_file, fasta_formats));

    fprintf(stderr, "\rLoaded reads.        \n");
    
    fastaq::FastAQ::print_statistics(reads, reads_file);

    map_single(ref_index, t_minimizers, reference, reads, parameters);
  }

  return 0;
}