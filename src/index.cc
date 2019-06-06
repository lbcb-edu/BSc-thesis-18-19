#include <cstdint>
#include <utility>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "index.hpp"

void prep_ref(std::vector<minimizer>& t_minimizers, const float f) {
  std::unordered_map<unsigned int, unsigned int> ref_min_frequency;
  for (const auto& minimizer : t_minimizers) {
    ref_min_frequency[std::get<0>(minimizer)]++;
  }

  std::vector<unsigned int> occurences;
  occurences.reserve(ref_min_frequency.size());

  for (const auto& entry : ref_min_frequency) {
    occurences.push_back(entry.second);
  }

  unsigned int position = (unsigned int)((1.0f - f) * (occurences.size() - 1.0f));
  std::sort(occurences.begin(), occurences.end());

  unsigned int cutoff_freq = occurences[position] == 1 ? 2 : occurences[position];

  std::vector<minimizer> temp;
  temp.reserve(t_minimizers.size());

  for (const auto& minimizer : t_minimizers) {
    if (ref_min_frequency[std::get<0>(minimizer)] < cutoff_freq) {
      temp.push_back(minimizer);
    }
  }

  std::swap(t_minimizers, temp);

  std::sort(t_minimizers.begin(), t_minimizers.end(),
      [] (const minimizer& a, const minimizer& b) {
        return (std::get<0>(a) < std::get<0>(b));
      });

  t_minimizers.shrink_to_fit();
}

std::unordered_map<unsigned int, minimizer_index_t> index_ref(const std::vector<minimizer>& t_minimizers) {
  std::unordered_map<unsigned int, minimizer_index_t> ref_index;
  unsigned int pos = 0;
  unsigned int num = 0;
  unsigned int prev_min = std::get<0>(t_minimizers[0]);
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