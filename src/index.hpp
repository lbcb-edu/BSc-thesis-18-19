#pragma once

// Minimizer: value, position, origin
typedef std::tuple<unsigned int, unsigned int, bool> minimizer;
// Index: position, amount
typedef std::pair<unsigned int, unsigned int> minimizer_index_t;

/**
 * Method removes top f frequent minimizers
 * 
 * @param t_minimizers	vector of minimizers
 * @param f 	percentage
 */
void prep_ref(std::vector<minimizer>& t_minimizers, const float f);

/**
 * Method returns map where key is minimizer value and value is position and amout of minimizers from given vector.
 * 
 * @param t_minimizers	vector of minimizers
 */
std::unordered_map<unsigned int, minimizer_index_t> index_ref(const std::vector<minimizer>& t_minimizers);