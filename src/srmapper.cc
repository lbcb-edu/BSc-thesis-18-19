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
#include <unordered_set>
#include <tuple>
#include <utility>
#include <algorithm>
#include <chrono>
#include <cstring>

#include "srmapper.hpp"
#include "fastaq.hpp"
#include "brown_minimizers.hpp"
#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "ksw2.h"

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

std::unordered_map<uint8_t, uint64_t> c = {{'C', 0}, {'A', 1}, {'T', 2}, {'U', 2}, {'G', 3}};
std::unordered_map<char, char> complement_map = {{'C', 'G'}, {'A', 'T'}, {'T', 'A'}, {'U', 'A'}, {'G', 'C'}};

static struct option long_options[] = {
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},
  {"paired", no_argument, NULL, 'p'},
  {"match", required_argument, NULL, 'm'},
  {"mismatch", required_argument, NULL, 'M'},
  {"gap-open", required_argument, NULL, 'o'},
  {"gap-extend", required_argument, NULL, 'e'},
  {"band", required_argument, NULL, 'b'},
  {"kmers", required_argument, NULL, 'k'},
  {"window_length", required_argument, NULL, 'w'},
  {"frequency", required_argument, NULL, 'f'},
  {"insert_size", required_argument, NULL, 'i'},
  {"region_size", required_argument, NULL, 'r'},
  {"threads", required_argument, NULL, 't'},
  {"threshold", required_argument, NULL, 'T'},
  {NULL, no_argument, NULL, 0}
};

typedef struct {
  int32_t mch;
  int32_t mis;
  int32_t gapo;
  int32_t gape;
  uint32_t band;
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
         "  reads     - one or two FASTA/FASTQ file containing a set of fragments\n\n"

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
         "  -p  or  --paired         use pairing information (needs insert size information)\n"
         "  -m  or  --match          <int>\n"
         "                             default: 1\n"
         "                             match value\n"
         "  -M  or  --mismatch       <int>\n"
         "                             default: -2\n"
         "                             mismatch value\n"
         "  -o  or  --gap-open       <int>\n"
         "                             default: 2\n"
         "                             gap open value\n"
         "  -e  or  --gap-extend     <int>\n"
         "                             default: 1\n"
         "                             gap extend value\n"
         "  -b  or  --band           <uint>\n"
         "                             default: -1\n"
         "                             ksw2 alignment band, band < 0 => disabled\n"
         "  -k  or  --kmers          <uint>\n"
         "                             default: 18\n"
         "                             constraints: largest supported is 32\n"
         "                             number of letters in substrings\n"
         "  -w  or  --window_length  <uint>\n"
         "                             default: 3\n"
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
         "                             default: 2\n"
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

std::string reverse_complement(const std::string& original, unsigned int pos, unsigned int length) {
  std::string rc(original.begin() + pos, original.begin() + pos + length);
  unsigned int j = pos + length - 1;
  for (unsigned int i = 0; i < length; ++i) {
    rc[i] = complement_map[original[j--]];
  }
  return rc;
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

paired_checked_t check_pairing(std::pair<bin_t, bin_t>& candidates, const uint32_t insert_size,
                               const uint32_t read_size, const uint32_t region_size) {
  paired_checked_t checked;
  for (const auto& bin : candidates.first) {
    if (!std::get<2>(bin.second[0])) {
      auto found = candidates.second.find(bin.first + ((insert_size - read_size) / region_size));
      if (found != candidates.second.end()) {
        checked.first.emplace_back(bin.second.begin(), bin.second.end());
        checked.second.emplace_back(found->second.begin(), found->second.end());
        continue;
      }
      found = candidates.second.find(bin.first + ((insert_size - read_size) / region_size) - 1);
      if (found != candidates.second.end()) {
        checked.first.emplace_back(bin.second.begin(), bin.second.end());
        checked.second.emplace_back(found->second.begin(), found->second.end());
        continue;
      }
      found = candidates.second.find(bin.first + ((insert_size - read_size) / region_size) + 1);
      if (found != candidates.second.end()) {
        checked.first.emplace_back(bin.second.begin(), bin.second.end());
        checked.second.emplace_back(found->second.begin(), found->second.end());
        continue;
      }
    } else {
      auto found = candidates.second.find(bin.first - ((insert_size - read_size) / region_size));
      if (found != candidates.second.end()) {
        checked.first.emplace_back(bin.second.begin(), bin.second.end());
        checked.second.emplace_back(found->second.begin(), found->second.end());
        continue;
      }
      found = candidates.second.find(bin.first - ((insert_size - read_size) / region_size) - 1);
      if (found != candidates.second.end()) {
        checked.first.emplace_back(bin.second.begin(), bin.second.end());
        checked.second.emplace_back(found->second.begin(), found->second.end());
        continue;
      }
      found = candidates.second.find(bin.first - ((insert_size - read_size) / region_size) + 1);
      if (found != candidates.second.end()) {
        checked.first.emplace_back(bin.second.begin(), bin.second.end());
        checked.second.emplace_back(found->second.begin(), found->second.end());
        continue;
      }
    }
  }
  return checked;
}

region_t find_region(std::vector<minimizer_hit_t>& hits) {
  if (hits.empty()) {
    return std::make_pair(std::make_tuple(0, 0, 0), std::make_tuple(0, 0, 0));
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

  unsigned int len = 1;
  lis[0] = std::make_pair(hits[0], hits[0]);

  for (unsigned int i = 1; i < hits.size(); ++i) {
    if (std::get<1>(hits[i]) > std::get<1>(lis[len - 1].second)) {
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
      }
    }
  }
  region_t region = lis[len - 1];
  if (std::get<2>(region.first) == 1) {
    unsigned int temp = std::get<0>(region.first);
    std::get<0>(region.first) = std::get<0>(region.second);
    std::get<0>(region.second) = temp;
  }
  return region; 
}

void expand_region(region_t& reg, uint32_t read_size, uint32_t k, uint32_t max_size) {
  if (std::get<2>(reg.first) == 0) {
    std::get<1>(reg.first) -= std::get<1>(reg.first) < std::get<0>(reg.first) ? std::get<1>(reg.first) : std::get<0>(reg.first);
    std::get<1>(reg.second) += max_size < std::get<1>(reg.second) + read_size - std::get<0>(reg.second) ? max_size - std::get<1>(reg.second) : read_size - std::get<0>(reg.second);
  } else {
    std::get<1>(reg.first) -= std::get<1>(reg.first) < read_size - std::get<0>(reg.second) - k ? std::get<1>(reg.first) : read_size - std::get<0>(reg.second) - k;
    std::get<1>(reg.second) += max_size < std::get<1>(reg.second) + std::get<0>(reg.first) + k ? max_size : std::get<0>(reg.first) + k;
  }
  std::get<0>(reg.first) = 0;
  std::get<0>(reg.second) = read_size - 1;
}

std::tuple<uint32_t, int32_t, std::string> ksw2(const std::string& target, const std::string& query, const region_t& region, 
                                                const mapping_params& parameters) {
  int a = parameters.mch;
  int b = parameters.mis;
  int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
  int len = query.size();
  uint8_t *ts, *qs;
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  ts = (uint8_t*)malloc(len);
  qs = (uint8_t*)malloc(len);
  for (int i = 0; i < len; ++i) {
    ts[i] = c[(uint8_t)target.c_str()[i + std::get<1>(region.first)]]; // encode to 0/1/2/3
    qs[i] = c[(uint8_t)query.c_str()[i]];
  }
  ksw_extz2_sse(0, len, qs, len, ts, 5, mat, parameters.gapo, parameters.gape, parameters.band, -1, 0, 0, &ez);
  std::string cigar;
  uint32_t matches = 0;
  uint32_t t_pos = std::get<1>(region.first);
  uint32_t q_pos = std::get<0>(region.first);
  for (int i = 0; i < ez.n_cigar; ++i) {
    if ("MIDNSHP=X"[ez.cigar[i]&0xf] == 'M') {
      for (uint32_t j = 0; j < (ez.cigar[i]>>4); ++j) {
        if (target[t_pos + j] == query[q_pos + j]) matches++;
      }
      t_pos += ez.cigar[i]>>4;
      q_pos += ez.cigar[i]>>4;
    } else if ("MIDNSHP=X"[ez.cigar[i]&0xf] == 'I') {
      q_pos += ez.cigar[i]>>4;
    } else if ("MIDNSHP=X"[ez.cigar[i]&0xf] == 'D') {
      t_pos += ez.cigar[i]>>4;
    }
    cigar += std::to_string(ez.cigar[i]>>4) + "MIDNSHP=X"[ez.cigar[i]&0xf];
  }
  free(ez.cigar); free(ts); free(qs);
  return std::make_tuple(matches, ez.score, cigar);
}

std::string sam_format_single(const std::string& qname, const std::string& query, const std::string& qual,
                              const std::string& rname, const std::string& ref, const region_t& region, 
                              const mapping_params& parameters) {
  std::string sam_name = qname.substr(0, qname.find('/', 0));
  std::tuple<uint32_t, int32_t, std::string> cigar = ksw2(ref, query, region, parameters);
  int flag = std::get<2>(region.first) ? 0x10 : 0x0;
  uint32_t mapq = std::min((int)((double)std::get<0>(cigar) / query.size() * 60), 60);
  std::string sam = sam_name + "\t" +
                    std::to_string(flag) + "\t" +
                    rname + "\t" +
                    std::to_string(std::get<1>(region.first) + 1) + "\t" +
                    std::to_string(mapq) + "\t" +
                    std::get<2>(cigar) + "\t" +
                    "*" + "\t" +
                    "0" + "\t" +
                    "0" + "\t" +
                    query + "\t" +
                    qual + "\t" +
                    "NM:i:" + std::to_string(query.size() - std::get<0>(cigar)) + "\t" +
                    "AS:i:" + std::to_string(std::get<1>(cigar)) + "\n";
  return sam;
}

std::pair<std::string, std::string> sam_format_pair(const std::string& qname,
    const std::string& query1, const std::string& qual1, const std::string& query2, const std::string& qual2,
    const std::string& rname, const std::string& ref,
    const std::pair<region_t, region_t>& region_pair, const mapping_params& parameters) {
  std::string sam_name = qname.substr(0, qname.find('/', 0));
  int32_t insert_size = std::get<1>(region_pair.first.first) < std::get<1>(region_pair.second.first)
                        ? std::get<1>(region_pair.second.second) - std::get<1>(region_pair.first.first)
                        : std::get<1>(region_pair.second.first) - std::get<1>(region_pair.first.second);
  std::tuple<uint32_t, int32_t, std::string> cigar1 = ksw2(ref, query1, region_pair.first, parameters);
  std::tuple<uint32_t, int32_t, std::string> cigar2 = ksw2(ref, query2, region_pair.second, parameters);
  int prop_aligned = std::get<0>(cigar1) + std::get<0>(cigar2) > 0.5 * insert_size ? 0x2 : 0x0;
  int flag1 = 0x1 | prop_aligned | (std::get<2>(region_pair.first.first) ? 0x10 : 0x0) | (std::get<2>(region_pair.second.first) ? 0x20 : 0x0) | 0x40;
  int flag2 = 0x1 | prop_aligned | (std::get<2>(region_pair.second.first) ? 0x10 : 0x0) | (std::get<2>(region_pair.first.first) ? 0x20 : 0x0) | 0x80;
  uint32_t mapq = std::min((int)((double)(std::get<0>(cigar1) + std::get<0>(cigar2)) / abs(insert_size) * 60), 60);
  std::string sam1 = sam_name + "\t" +
                     std::to_string(flag1) + "\t" +
                     rname + "\t" +
                     std::to_string(std::get<1>(region_pair.first.first) + 1) + "\t" +
                     std::to_string(mapq) + "\t" +
                     std::get<2>(cigar1) + "\t" +
                     "=" + "\t" +
                     std::to_string(std::get<1>(region_pair.second.first) + 1) + "\t" +
                     std::to_string(insert_size) + "\t" +
                     query1 + "\t" +
                     qual1 + "\t" +
                     "NM:i:" + std::to_string(query1.size() - std::get<0>(cigar1)) + "\t" +
                     "AS:i:" + std::to_string(std::get<1>(cigar1)) + "\n";
  std::string sam2 = sam_name + "\t" +
                     std::to_string(flag2) + "\t" +
                     rname + "\t" +
                     std::to_string(std::get<1>(region_pair.second.first) + 1) + "\t" +
                     std::to_string(mapq) + "\t" +
                     std::get<2>(cigar2) + "\t" +
                     "=" + "\t" +
                     std::to_string(std::get<1>(region_pair.first.first) + 1) + "\t" +
                     std::to_string(-insert_size) + "\t" +
                     query2 + "\t" +
                     qual2 + "\t" +
                     "NM:i:" + std::to_string(query2.size() - std::get<0>(cigar2)) + "\t" +
                     "AS:i:" + std::to_string(std::get<1>(cigar2)) + "\n";
  return std::make_pair(sam1, sam2);
}

std::string unmapped_sam(const std::string& qname, const std::string& query, const std::string& qual, bool pair, bool first, bool last) {
  std::string sam_name = qname.substr(0, qname.find('/', 0));
  int flag = 0x4;
  if (pair) {
    flag |= 0x1 | 0x8;
    if (first) flag |= 0x40;
    if (last) flag |= 0x80;
  }
  std::string sam = sam_name + "\t" +
                    std::to_string(flag) + "\t" +
                    "*" + "\t" +
                    "0" + "\t" +
                    "255" + "\t" +
                    "*" + "\t" +
                    "*" + "\t" +
                    "0" + "\t" +
                    "0" + "\t" +
                    query + "\t" +
                    qual + "\n";
  return sam;
}

void process_pairs(std::vector<std::pair<uint32_t, std::string>>& mappings, paired_checked_t& checked,
                   const std::vector<std::unique_ptr<fastaq::FastAQ>>& reference, 
                   const paired_reads_t& paired_reads, uint32_t i, const mapping_params& parameters) {
  std::unordered_set<uint32_t> processed;
  for (uint32_t j = 0; j < checked.first.size(); ++j) {
    std::pair<region_t, region_t> region_pair(find_region(checked.first[j]),
                                              find_region(checked.second[j]));
    expand_region(region_pair.first, paired_reads.first[i]->sequence.size(), 
                  parameters.k, reference[0]->sequence.size() - 1);
    expand_region(region_pair.second, paired_reads.second[i]->sequence.size(), 
                  parameters.k, reference[0]->sequence.size() - 1);
    if (processed.find(std::get<1>(region_pair.first.first)) != processed.end()) continue;
    processed.insert(std::get<1>(region_pair.first.first));
    std::string rc;
    std::string rq;
    std::string* query1 = std::get<2>(region_pair.first.first)
                          ? &(rc = reverse_complement(paired_reads.first[i]->sequence, 0, paired_reads.first[i]->sequence.size()))
                          : &(paired_reads.first[i]->sequence);
    std::string* quality1 = std::get<2>(region_pair.first.first)
                            ? &(rq = std::string(paired_reads.first[i]->quality.rbegin(), paired_reads.first[i]->quality.rend()))
                            : &(paired_reads.first[i]->quality);
    std::string* query2 = std::get<2>(region_pair.second.first)
                          ? &(rc = reverse_complement(paired_reads.second[i]->sequence, 0, paired_reads.second[i]->sequence.size()))
                          : &(paired_reads.second[i]->sequence);
    std::string* quality2 = std::get<2>(region_pair.second.first)
                            ? &(rq = std::string(paired_reads.second[i]->quality.rbegin(), paired_reads.second[i]->quality.rend()))
                            : &(paired_reads.second[i]->quality);
    std::pair<std::string, std::string> sam_pair = sam_format_pair(
        paired_reads.first[i]->name,
        *query1, *quality1, 
        *query2, *quality2, 
        reference[0]->name, reference[0]->sequence, region_pair, parameters);

    mappings.emplace_back(0, sam_pair.first);
    mappings.emplace_back(0, sam_pair.second);
  }
}

void map_single(const std::unordered_map<uint64_t, index_pos_t>& ref_index, const std::vector<minimizer_t>& t_minimizers,
                const std::vector<std::unique_ptr<fastaq::FastAQ>>& reference, 
                const std::vector<std::unique_ptr<fastaq::FastAQ>>& reads, const mapping_params& parameters) {
  std::string sam;
  for (uint32_t i = 0; i < reads.size(); ++i) {
    std::vector<minimizer_t> q_minimizers = brown::minimizers(reads[i]->sequence.c_str(), reads[i]->sequence.size(),
                                                             parameters.k, parameters.w);
    split_hits_t hits;
    find_minimizer_hits(hits.first, hits.second, ref_index, t_minimizers, q_minimizers);
    radixsort(hits.first);
    radixsort(hits.second);
    std::pair<bin_t, bin_t> candidates(extract_candidates(hits.first, parameters.threshold, parameters.region_size),
                                       extract_candidates(hits.second, parameters.threshold, parameters.region_size));
    if (candidates.first.size() == 0 && candidates.second.size() == 0) {
      sam += unmapped_sam(reads[i]->name, reads[i]->sequence, reads[i]->quality, 0, 0, 0);
      continue;
    }
    std::unordered_set<uint32_t> processed;
    std::vector<std::pair<uint32_t, std::string>> mappings;
    for (auto& bin : candidates.first) {
      auto reg = find_region(bin.second);
      expand_region(reg, reads[i]->sequence.size(), parameters.k, reference[0]->sequence.size() - 1);
      if (processed.find(std::get<1>(reg.first)) != processed.end()) continue;
      processed.insert(std::get<1>(reg.first));
      mappings.emplace_back(0, sam_format_single(reads[i]->name, reads[i]->sequence, reads[i]->quality,
                            reference[0]->name, reference[0]->sequence, reg, 
                            parameters));
    }
    for (auto& bin : candidates.second) {
      auto reg = find_region(bin.second);
      expand_region(reg, reads[i]->sequence.size(), parameters.k, reference[0]->sequence.size() - 1);
      if (processed.find(std::get<1>(reg.first)) != processed.end()) continue;
      processed.insert(std::get<1>(reg.first));
      std::string rc = reverse_complement(reads[i]->sequence, 0, reads[i]->sequence.size());
      std::string rq = std::string(reads[i]->quality.rbegin(), reads[i]->quality.rend());
      mappings.emplace_back(0, sam_format_single(reads[i]->name, rc, rq,
                                                 reference[0]->name, reference[0]->sequence, reg, 
                                                 parameters));
    }
    for (const auto& m : mappings) {
      sam += m.second;
    }
  }
  std::cout << sam;
}

void map_paired(const std::unordered_map<uint64_t, index_pos_t>& ref_index, const std::vector<minimizer_t>& t_minimizers,
                const std::vector<std::unique_ptr<fastaq::FastAQ>>& reference, const paired_reads_t& paired_reads,
                const mapping_params& parameters) {
  std::string sam;
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
        extract_candidates(hits1.second, parameters.threshold, parameters.region_size),
        extract_candidates(hits2.first, parameters.threshold, parameters.region_size));
    paired_checked_t checked1 = check_pairing(candidates1, parameters.insert_size, parameters.region_size,
                                              paired_reads.second[i]->sequence.size());
    paired_checked_t checked2 = check_pairing(candidates2, parameters.insert_size, parameters.region_size,
                                              paired_reads.first[i]->sequence.size());

    std::vector<std::pair<uint32_t, std::string>> mappings;
    process_pairs(mappings, checked1, reference, paired_reads, i, parameters);
    process_pairs(mappings, checked2, reference, paired_reads, i, parameters);
    if (mappings.size() == 0) {
      sam += unmapped_sam(paired_reads.first[i]->name, paired_reads.first[i]->sequence, paired_reads.first[i]->quality, 1, 1, 0)
             + unmapped_sam(paired_reads.second[i]->name, paired_reads.second[i]->sequence, paired_reads.second[i]->quality, 1, 0, 1);
      continue;
    }
    for (const auto& m : mappings) {
      sam += m.second;
    }
  }
  std::cout << sam;
}

int main(int argc, char **argv) {
  int optchr;
  mapping_params parameters;
  parameters.mch = 1;
  parameters.mis = -2;
  parameters.gapo = 2;
  parameters.gape = 1;
  parameters.band = -1;
  parameters.k = 18;
  parameters.w = 3;
  parameters.f = 0.001f;
  parameters.insert_size = 215;
  parameters.region_size = 100;
  parameters.threshold = 2;
  bool paired = false;
  uint32_t threads = 3;

  while ((optchr = getopt_long(argc, argv, "hvpm:M:o:e:b:k:w:f:i:r:t:T:", long_options, NULL)) != -1) {
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
      case 'm': {
        parameters.mch = atoi(optarg);
        break;
      }
      case 'M': {
        parameters.mis = atoi(optarg);
        break;
      }
      case 'o': {
        parameters.gapo = atoi(optarg);
        break;
      }
      case 'e': {
        parameters.gape = atoi(optarg);
        break;
      }
      case 'b': {
        parameters.band = atoi(optarg);
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

  auto i_start = std::chrono::steady_clock::now();

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

  auto i_end = std::chrono::steady_clock::now();

  auto i_interval = std::chrono::duration_cast<std::chrono::duration<double>>(i_end - i_start);

  std::cerr << "\nIndexing time: " << i_interval.count() << " sec" << std::endl;
  
  auto m_start = std::chrono::steady_clock::now();

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
      map_single(ref_index, t_minimizers, reference, paired_reads.first, parameters);
      map_single(ref_index, t_minimizers, reference, paired_reads.second, parameters);
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

  auto m_end = std::chrono::steady_clock::now();

  auto m_interval = std::chrono::duration_cast<std::chrono::duration<double>>(m_end - m_start);

  std::cerr << "\nMapping time: " << m_interval.count() << " sec" << std::endl;

  return 0;
}