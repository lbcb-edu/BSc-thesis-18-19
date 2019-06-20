#pragma once

// Minimizer: value, position, origin
typedef std::tuple<unsigned int, unsigned int, bool> minimizer;
// Index: position, amount
typedef std::pair<unsigned int, unsigned int> minimizer_index_t;
// Minimizer hit: query position, reference position
typedef std::pair<unsigned int, unsigned int> minimizer_hit_t;
// Region hit: region number, hits number
typedef std::pair<unsigned int, int> region_hits;

/**
 * Method finds minimizer hits, checks strand and stores that match in last two arguments
 * 
 * @param ref_index	map of minimizers
 * @param t_minimizers target minimizers
 * @param q_minimizers query minimizers
 * @param same vector where match is stored if strands are same
 * @param rev vector where match is stored if strands are diffrent
 */
void find_minimizer_hits(
    const std::unordered_map<unsigned int, minimizer_index_t>& ref_index,
    const std::vector<minimizer>& t_minimizers,
    const std::vector<minimizer>& q_minimizers,
    std::vector<minimizer_hit_t>& same,
    std::vector<minimizer_hit_t>& rev);


/**
 * Method returns map where key is region number and value is number of matches in that particular region
 */
std::unordered_map<unsigned int, int> parse_hits_map(
	std::vector<minimizer_hit_t> & hits,
	int k_value);

/**
 * Method returns maximum 15 regions which have highest match number. Returnd value is pair where first
 * value is region number and second is number of matches
 */
std::vector<region_hits> find_top_15(
	std::unordered_map<unsigned int, int>& map,
	int threshold);

/**
 * Method returns vector containing pair where first value is region number and second is number of matches
 */
std::vector<region_hits> find_top_all(
	std::unordered_map<unsigned int,
	int>& map);

/**
 * Method returns region with most matches
 */
region_hits find_top(
	std::unordered_map<unsigned int, int>& map,
	int threshold);

/**
 * Method finds and stores starting and ending region in given arguments
 */
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

/**
 * Method find and stores starting and ending position on query
 * and starting and ending position on reference
 */
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

/**
 * Method checks which parts where found while checing matches.
 * Result is stored in given boolean arguments
 */
void find_what_exists(
	std::vector<region_hits>& start_hits_top, 
	std::vector<region_hits>& start_hits_top_rev,
  	std::vector<region_hits>& end_hits_top, 
  	std::vector<region_hits>& end_hits_top_rev, 
  	bool& bool_start, bool& bool_end);