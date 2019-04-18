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
// paired read
typedef std::pair<std::vector<std::unique_ptr<fastaq::FastAQ>>, std::vector<std::unique_ptr<fastaq::FastAQ>>> paired_reads_t;
// paired read minimizers
typedef std::pair<std::vector<minimizer_t>, std::vector<minimizer_t>> paired_minimizers_t;
// paired read hits
typedef std::pair<std::vector<minimizer_hit_t>, std::vector<minimizer_hit_t>> paired_hits_t;

const std::set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
const std::set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};

static struct option long_options[] = {
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},
  {"window_length", required_argument, NULL, 'w'},
  {"kmers", required_argument, NULL, 'k'},
  {"frequency", required_argument, NULL, 'f'},
  {"insert_size", required_argument, NULL, 'i'},
  {"region_size", required_argument, NULL, 'r'},
  {"threshold", required_argument, NULL, 't'},
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
         "                             default: 215\n"
         "                             region size to count hits from\n"
         "  -t  or  --threshold      <uint>\n"
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

void radixsort(std::vector<uint32_t>& array) {
  if (array.size() <= 1) return;
  uint32_t max_val = array[0];
  const uint32_t size = array.size();
  for (uint32_t i = 1; i < size; ++i) {
    if (array[i] > max_val) {
      max_val = array[i];
    }
  }
  std::vector<uint32_t> out(size);
  for (uint32_t e = 1; max_val / e > 0; e *= 10) {
    uint32_t count[10];
    for (uint32_t i = 0; i < 10; ++i) {
      count[i] = 0;
    }
    for (uint32_t i = 0; i < size; ++i) {
      count[(array[i] / e) % 10]++;
    }
    for (uint32_t i = 1; i < 10; ++i) {
      count[i] += count[i - 1];
    }
    for (uint32_t i = size - 1; i != (uint32_t)(-1); --i) {
      out[count[(array[i] / e) % 10] - 1] = array[i];
      count[(array[i] / e) % 10]--;
    }
    for (uint32_t i = 0; i < size; ++i) {
      array[i] = out[i];
    }
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

std::vector<minimizer_hit_t> find_minimizer_hits(
    const std::unordered_map<uint64_t, index_pos_t>& ref_index,
    const std::vector<minimizer_t>& t_minimizers,
    const std::vector<minimizer_t>& q_minimizers) {

  std::vector<minimizer_hit_t> hits;

  for (const auto& minimizer : q_minimizers) {
    auto found = ref_index.find(std::get<0>(minimizer));
    if (found != ref_index.end()) {
      for (uint32_t i = 0; i < found->second.second; i++) {
        if (std::get<2>(minimizer) == std::get<2>(t_minimizers[found->second.first + i])) {
          hits.emplace_back(
              std::get<1>(minimizer),
              std::get<1>(t_minimizers[found->second.first + i]),
              0);  
        } else {
          hits.emplace_back(
              std::get<1>(minimizer),
              std::get<1>(t_minimizers[found->second.first + i]),
              1);
        }
      }
    }
  }

  return hits;
}

std::vector<uint32_t> extract_candidates(const std::vector<uint32_t>& reg_hits, const uint32_t threshold) {
  std::vector<uint32_t> candidates;
  uint32_t current = 0;
  uint32_t count = 0;
  uint32_t next = 0;
  uint32_t next_count = 0;
  for (uint32_t i = 0; i < reg_hits.size(); ++i) {
    // std::cerr << "checking " << reg_hits[i] << " | current " << current << ", count " << count << ", next " << next << ", next_count " << next_count << std::endl;
    if (reg_hits[i] == current) {
      // std::cerr << "  increase count" << std::endl;
      count++;
      continue;
    }
    if (reg_hits[i] == current + 1) {
      // std::cerr << "  increase both counts, set next" << std::endl;
      count++;
      next_count++;
      next = current + 1;
      continue;
    }
    if (count >= threshold) {
      candidates.push_back(current);
    }
    if (next == current + 1) {
      current = next;
      count = next_count;
    } else {
      current = reg_hits[i];
      count = 1;
    }
  }
  if (count >= threshold) {
    candidates.push_back(current);
  }
  return candidates;
}

void map_paired(const std::unordered_map<uint64_t, index_pos_t>& ref_index, const std::vector<minimizer_t>& t_minimizers,
         const std::vector<std::unique_ptr<fastaq::FastAQ>>& reference, const paired_reads_t& paired_reads,
         const mapping_params& parameters) {
  uint32_t n_reads_with_candidates = 0;
  for (uint32_t i = 0; i < paired_reads.first.size(); ++i) {
    paired_minimizers_t pminimizers(brown::minimizers(paired_reads.first[i]->sequence.c_str(),
                                                      paired_reads.first[i]->sequence.size(),
                                                      parameters.k, parameters.w),
                                    brown::minimizers(paired_reads.second[i]->sequence.c_str(),
                                                      paired_reads.second[i]->sequence.size(),
                                                      parameters.k, parameters.w));
    std::pair<std::vector<uint32_t>, std::vector<uint32_t>> phits(
        find_ref_region_hits(ref_index, t_minimizers, pminimizers.first, parameters.region_size),
        find_ref_region_hits(ref_index, t_minimizers, pminimizers.second, parameters.region_size));
    // std::cerr << "phits first size " << phits.first.size() << std::endl;
    // for (const auto& phf : phits.first) {
    //   std::cerr << phf << std::endl;
    // }
    // std::cerr << "phits second size " << phits.second.size() << std::endl;
    // for (const auto& phs : phits.second) {
    //   std::cerr << phs << std::endl;
    // }
    // std::cerr << std::endl;

    radixsort(phits.first);
    radixsort(phits.second);

    // for (const auto& phf : phits.first) {
    //   std::cerr << phf << std::endl;
    // }
    // std::cerr << std::endl;
    // for (const auto& phs : phits.second) {
    //   std::cerr << phs << std::endl;
    // }
    // std::cerr << std::endl;

    std::pair<std::vector<uint32_t>, std::vector<uint32_t>> candidates(
        extract_candidates(phits.first, parameters.threshold),
        extract_candidates(phits.second, parameters.threshold));
    if (candidates.first.size() && candidates.second.size()) n_reads_with_candidates++;
  }
  std::cerr << n_reads_with_candidates << std::endl;
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

  while ((optchr = getopt_long(argc, argv, "hvk:w:f:i:r:t:", long_options, NULL)) != -1) {
    switch (optchr) {
      case 'h': {
        help();
        exit(0);
      }
      case 'v': {
        version();
        exit(0);
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
    fprintf(stderr, "[srmapper] error: Expected read(s) and reference! Use --help for usage.\n");
    exit(1);
  }

  fprintf(stderr, "\nMapping process started with parameters:\n"
                  "  k             = %u\n"
                  "  window length = %u\n"
                  "  top f freq    = %g\n"
                  "  insert size   = %u\n"
                  "  region size   = %u\n"
                  "  threshold     = %u\n",
                  parameters.k, parameters.w, parameters.f,
                  parameters.insert_size, parameters.region_size, parameters.threshold);

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

    if (pr1_stats.num != pr2_stats.num) {
      fprintf(stderr, "[srmapper] error: Paired-end read files must have equal number of reads (pairs).\n");
      exit(1);
    }
    if (pr1_stats.max - pr1_stats.min > 0.5 * pr1_stats.avg || pr2_stats.max - pr2_stats.min > 0.5 * pr2_stats.avg) {
      fprintf(stderr, "[srmapper] warning: Reads are not of fixed size.\n");
    }

    map_paired(ref_index, t_minimizers, reference, paired_reads, parameters);
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

    // map_single(ref_index, t_minimizers, reference, reads, k, w);
  }

  return 0;
}