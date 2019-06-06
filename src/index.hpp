#pragma once

// Minimizer: value, position, origin
typedef std::tuple<unsigned int, unsigned int, bool> minimizer;
// Index: position, amount
typedef std::pair<unsigned int, unsigned int> minimizer_index_t;

void prep_ref(std::vector<minimizer>& t_minimizers, const float f);

std::unordered_map<unsigned int, minimizer_index_t> index_ref(const std::vector<minimizer>& t_minimizers);