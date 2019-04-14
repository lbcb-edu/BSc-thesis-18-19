#include <iostream>
#include <vector>
#include <unordered_set>
#include <tuple>
#include <queue>
#include <limits>
#include <algorithm>
#include <functional>
#include <map>

#include "brown_minimizers.hpp"

namespace brown {

std::map<char, int> char_to_val = {{'C', 0}, {'A', 1}, {'T', 2}, {'U', 2}, {'G', 3}};

typedef std::tuple<unsigned int, unsigned int, bool> triplet_t;

struct triplet_hash {
  std::size_t operator() (const triplet_t& k) const noexcept {
    return std::hash<unsigned int>{}(std::get<0>(k)) ^ std::hash<unsigned int>{}(std::get<1>(k)) ^ std::hash<bool>{}(std::get<2>(k));
  }
};

struct triplet_equal {
  bool operator() (const triplet_t& v0, const triplet_t& v1) const noexcept {
    return (
            std::get<0>(v0) == std::get<0>(v1) &&
            std::get<1>(v0) == std::get<1>(v1) &&
            std::get<2>(v0) == std::get<2>(v1)
           );
  }
};

struct triplet_ordering {
  bool operator() (const triplet_t& v0, const triplet_t& v1) const noexcept {
    return (std::get<1>(v0) < std::get<1>(v1));
  }
};

inline unsigned int value(const char* sequence, unsigned int pos, unsigned int k) {
  unsigned int val = 0;
  for (unsigned int i = 0; i < k; ++i) {
    val = (val << 2) | char_to_val[sequence[pos + i]];
  }
  return val;
}

inline unsigned int value_reverse_complement(const char* sequence,
                                      unsigned int pos,
                                      unsigned int k) {

  unsigned int val = 0;
  for (unsigned int i = k - 1; i < (unsigned int)(-1); --i) {
    val = (val << 2) | ((~char_to_val[sequence[pos + i]] & 3));
  }
  return val;
}

void interior_minimizers_fill(
    std::unordered_set<triplet_t, triplet_hash, triplet_equal>& minimizers_set,
    const char* sequence,
    unsigned int sequence_length,
    unsigned int k,
    unsigned int window_length) {

  unsigned int mask = (1 << (2 * k)) - 1;
  unsigned int original;
  unsigned int rev_com;
  unsigned int foriginal;
  unsigned int frev_com;

  bool cache_empty = true;
  triplet_t cache;
  for (unsigned int i = 0; i <= sequence_length - window_length - k + 1; ++i) {
    unsigned int m = std::numeric_limits<unsigned int>::max();

    if (cache_empty) {
      foriginal = value(sequence, i, k) >> 2;
      frev_com = (value_reverse_complement(sequence, i, k) << 2) & mask;
    } else {
      foriginal = original;
      frev_com = rev_com;
    }
    original = foriginal;
    rev_com = frev_com;

    unsigned int start;
    start = cache_empty ? 0 : window_length - 1;

    for (unsigned int j = start; j < window_length; ++j) {
      original = ((original << 2) | char_to_val[sequence[k - 1 + i + j]]) & mask;
      rev_com = ((rev_com >> 2) | ((~char_to_val[sequence[k - 1 + i + j]] & 3) << (2 * (k - 1)))) & mask;
      if (original != rev_com) {
        m = std::min(m, std::min(original, rev_com));
      }
    }

    if (!cache_empty) {
      if (std::get<1>(cache) >= i && std::get<0>(cache) <= m) {
        m = std::get<0>(cache);
      } else {
        cache_empty = true;
      }
    }

    original = foriginal;
    rev_com = frev_com;

    for (unsigned int j = start; j < window_length; ++j) {
      original = ((original << 2) | char_to_val[sequence[k - 1 + i + j]]) & mask;
      rev_com = ((rev_com >> 2) | ((~char_to_val[sequence[k - 1 + i + j]] & 3) << (2 * (k - 1)))) & mask;
      if (original < rev_com && original == m) {
        minimizers_set.emplace(m, i + j, 0);
        cache = std::make_tuple(m, i + j, 0);
        cache_empty = false;
      } else if (rev_com < original && rev_com == m) {
        minimizers_set.emplace(m, i + j, 1);
        cache = std::make_tuple(m, i + j, 1);
        cache_empty = false;
      }
    }
  }
}

void end_minimizers_fill(
    std::unordered_set<triplet_t, triplet_hash, triplet_equal>& minimizers_set,
    unsigned int k,
    unsigned int window_length,
    std::function<unsigned int(unsigned int)> position,
    std::function<unsigned int(unsigned int, unsigned int)> start_o,
    std::function<unsigned int(unsigned int, unsigned int)> start_r,
    std::function<unsigned int(unsigned int, int, unsigned int)> insert_o,
    std::function<unsigned int(unsigned int, int, unsigned int)> insert_r) {

  unsigned int mask = (1 << (2 * k)) - 1;
  unsigned int original;
  unsigned int rev_com;
  unsigned int m = std::numeric_limits<unsigned int>::max();

  for (unsigned int i = 0; i < window_length - 1; ++i) {
    if (i == 0) {
      original = start_o(position(i), mask);
      rev_com = start_r(position(i), mask);
    }

    original = insert_o(original, position(i), mask);
    rev_com = insert_r(rev_com, position(i), mask);

    if (original != rev_com) {
      m = std::min(m, std::min(original, rev_com));
    }

    if (original < rev_com && original == m) {
      minimizers_set.emplace(m, position(i), 0);
    } else if (rev_com < original && rev_com == m) {
      minimizers_set.emplace(m, position(i), 1);
    }
  }
}

std::vector<triplet_t> minimizers(
    const char* sequence,
    unsigned int sequence_length,
    unsigned int k,
    unsigned int window_length) {

  if (k > 16) {
    fprintf(stderr, "[brown::minimizers] error: Largest supported value for k is 16!\n"
                    "  k = %d.\n", k);
    exit(1);
  }
  if (sequence_length < window_length + k - 1) {
    fprintf(stderr, "[brown::minimizers] error: Sequence length too short for given parameters!\n"
                    "  Length = %d, k = %d, window length = %d.\n",
        sequence_length, k, window_length);
    exit(1);
  }

  std::unordered_set<triplet_t, triplet_hash, triplet_equal> minimizers_set;
  
  interior_minimizers_fill(minimizers_set, sequence, sequence_length, k, window_length);

  end_minimizers_fill(minimizers_set, k, window_length,
                      [] (unsigned int i) {return i;},
                      [&k, &sequence] (unsigned int pos, unsigned int mask) {
                        return (value(sequence, pos, k) >> 2) & mask;
                      },
                      [&k, &sequence] (unsigned int pos, unsigned int mask) {
                        return (value_reverse_complement(sequence, pos, k) << 2) & mask;
                      },
                      [&k, &sequence] (unsigned int original, unsigned int pos, unsigned int mask) {
                        return ((original << 2) | char_to_val[sequence[k - 1 + pos]]) & mask;
                      },
                      [&k, &sequence] (unsigned int rev_com, unsigned int pos, unsigned int mask) {
                        return ((rev_com >> 2) | ((~char_to_val[sequence[k - 1 + pos]] & 3) << (2 * (k - 1)))) & mask;
                      }
                     );

  end_minimizers_fill(minimizers_set, k, window_length,
                      [&sequence_length, &k] (unsigned int i) {return sequence_length - k - i;},
                      [&k, &sequence] (unsigned int pos, unsigned int mask) {
                        return (value(sequence, pos, k) << 2) & mask;
                      },
                      [&k, &sequence] (unsigned int pos, unsigned int mask) {
                        return (value_reverse_complement(sequence, pos, k) >> 2) & mask;
                      },
                      [&k, &sequence] (unsigned int original, unsigned int pos, unsigned int mask) {
                        return ((original >> 2) | (char_to_val[sequence[pos]] << (2 * (k - 1)))) & mask;
                      },
                      [&sequence] (unsigned int rev_com, unsigned int pos, unsigned int mask) {
                        return ((rev_com << 2) | (~char_to_val[sequence[pos]] & 3)) & mask;
                      }
                     );

  std::vector<triplet_t> minimizers_vector(
      minimizers_set.begin(), minimizers_set.end());

  return minimizers_vector;
}

}