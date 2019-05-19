#pragma once

// value, position, origin
typedef std::tuple<unsigned int, unsigned int, bool> triplet_t;
// position, range
typedef std::pair<unsigned int, unsigned int> minimizer_index_t;

void prep_ref(std::vector<triplet_t>& t_minimizers, const float f);

std::unordered_map<unsigned int, minimizer_index_t> index_ref(const std::vector<triplet_t>& t_minimizers);