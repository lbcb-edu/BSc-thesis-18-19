#pragma once

// Minimizer: value, position, origin
typedef std::tuple<unsigned int, unsigned int, bool> minimizer;
// Index: position, amount
typedef std::pair<unsigned int, unsigned int> minimizer_index_t;
// Minimizer hit: query position, reference position
typedef std::pair<unsigned int, unsigned int> minimizer_hit_t;
// Region hit: region number, hits number
typedef std::pair<unsigned int, int> region_hits;


void find_minimizer_hits(
    const std::unordered_map<unsigned int, minimizer_index_t>& ref_index,
    const std::vector<minimizer>& t_minimizers,
    const std::vector<minimizer>& q_minimizers,
    std::vector<minimizer_hit_t>& same,
    std::vector<minimizer_hit_t>& rev);

std::unordered_map<unsigned int, int> parse_hits_map(std::vector<minimizer_hit_t> & hits, int k_value);

std::vector<region_hits> find_top_3(std::unordered_map<unsigned int, int>& map, int threshold);

std::vector<region_hits> find_top_all(std::unordered_map<unsigned int, int>& map);

region_hits find_top(std::unordered_map<unsigned int, int>& map, int threshold);

void find_regions(
	std::vector<region_hits>& start_hits_top,
	std::vector<region_hits>& start_hits_top_rev,
  	std::vector<region_hits>& end_hits_top,
  	std::vector<region_hits>& end_hits_top_rev,
  	int& hits_number,
  	unsigned int& starting_region,
  	unsigned int& ending_region,
  	bool& rev,
  	int region_size);

void find_positions(
	std::vector<minimizer_hit_t>& start_hits_same,
	std::vector<minimizer_hit_t>& start_hits_rev,
	std::vector<minimizer_hit_t>& end_hits_same,
	std::vector<minimizer_hit_t>& end_hits_rev,
	bool rev, 
	unsigned int starting_region,
	unsigned int ending_region, 
	unsigned int& start,
	unsigned int& end,
	unsigned int& ref_start,
	unsigned int& ref_end,
	int k_value);

void find_what_exists(
	std::vector<region_hits>& start_hits_top, 
	std::vector<region_hits>& start_hits_top_rev,
  	std::vector<region_hits>& end_hits_top, 
  	std::vector<region_hits>& end_hits_top_rev, 
  	bool& bool_start, bool& bool_end);