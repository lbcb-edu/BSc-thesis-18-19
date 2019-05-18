#pragma once

#include "mapping_params.hpp"

// Hit: query position, reference position, relative strand
typedef std::tuple<uint32_t, uint32_t, bool> minimizer_hit_t;
// Region: two "hits" with swapped positions for reverse complement that "map" the query to the target
typedef std::pair<minimizer_hit_t, minimizer_hit_t> region_t;

std::pair<uint32_t, std::string> sam_format_single(const std::string& qname, const std::string& query, 
                                                   const std::string& qual, const std::string& rname, 
                                                   const std::string& ref, const region_t& region, 
                                                   const mapping_params& parameters, const int32_t clipped);

std::pair<uint32_t, std::pair<std::string, std::string>> sam_format_pair(const std::string& qname,
    const std::string& query1, const std::string& qual1, const std::string& query2, const std::string& qual2,
    const std::string& rname, const std::string& ref,
    const std::pair<region_t, region_t>& region_pair, const mapping_params& parameters,
    const int32_t clipped1, const int32_t clipped2);

std::string unmapped_sam(const std::string& qname, const std::string& query, const std::string& qual, 
                         bool pair, bool first, bool last);