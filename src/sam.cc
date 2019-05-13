#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <algorithm>

#include "sam.hpp"
#include "mapping_params.hpp"
#include "ksw2.h"

std::unordered_map<uint8_t, uint64_t> c = {{'C', 0}, {'A', 1}, {'T', 2}, {'U', 2}, {'G', 3}};

// Perform KSW2 algorithm on found region, return cigar
// Args: target     - target sequence
//       query      - query sequence
//       region     - region on query and target
//       parameters - wrapper of significant mapping parameters
// Return: number of exact matches, alignment score, cigar string
std::tuple<uint32_t, int32_t, std::string> ksw2(const std::string& target, const std::string& query, 
                                                const region_t& region, const mapping_params& parameters) {
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

// Generates mapping in SAM format
// Args: qname      - query name from read file
//       query      - query sequence
//       qual       - query sequence per base qualities
//       rname      - reference name from reference file
//       ref        - reference sequence
//       region     - mapping region
//       parameters - mapping parameters
// Return: mapping in SAM format
std::pair<uint32_t, std::string> sam_format_single(const std::string& qname, const std::string& query, 
                                                   const std::string& qual, const std::string& rname, 
                                                   const std::string& ref, const region_t& region, 
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
  return std::make_pair(mapq, sam);
}

// Generates paired mapping in SAM format
// Args: qname       - query name from read file
//       query1      - query sequence from first read
//       qual1       - query sequence per base qualities from first read
//       query2      - query sequence from second read
//       qual2       - query sequence per base qualities from second read
//       rname       - reference name from reference file
//       ref         - reference sequence
//       region_pair - paired mapping region
//       parameters  - mapping parameters
// Return: paired mapping in SAM format
std::pair<uint32_t, std::pair<std::string, std::string>> sam_format_pair(const std::string& qname,
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
  int flag1 = 0x1 | prop_aligned | (std::get<2>(region_pair.first.first) ? 0x10 : 0x0) 
                  | (std::get<2>(region_pair.second.first) ? 0x20 : 0x0) | 0x40;
  int flag2 = 0x1 | prop_aligned | (std::get<2>(region_pair.second.first) ? 0x10 : 0x0) 
                  | (std::get<2>(region_pair.first.first) ? 0x20 : 0x0) | 0x80;
  uint32_t mapq = std::min((int)((double)(std::get<0>(cigar1) + std::get<0>(cigar2)) / (query1.size() + query2.size()) * 60), 60);
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
  return std::make_pair(mapq, std::make_pair(sam1, sam2));
}

// Generates SAM format output for unmapped read
// Args: qname - query name from read file
//       query - query sequence
//       qual  - query sequence per base qualities
//       pair  - paired reads flag
//       first - first in segment flag
//       last  - last in segment flag
// Return: unmapped read SAM format output
std::string unmapped_sam(const std::string& qname, const std::string& query, const std::string& qual, 
                         bool pair, bool first, bool last) {
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