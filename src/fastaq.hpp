#pragma once

#include <string>
#include <vector>
#include <limits>
#include <cstdint>

#include "bioparser/bioparser.hpp"

/**
 * Class that represents FASTA and FASTAQ format
 */
class FastAQ {
public:
  std::string name;
  std::string sequence;
  std::string quality;

  FastAQ(
    const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length) : FastAQ(name, name_length, sequence, sequence_length, "", 0){}

  FastAQ(
    const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length,
    const char* quality, uint32_t quality_length) {
      this->name = {name, name_length};
      this->sequence = {sequence, sequence_length};
      this->quality = {quality, quality_length};
  }

  /**
   * Method print statistics about file
   */
  static void print_statistics(const std::vector<std::unique_ptr<FastAQ>> &fastaq_objects, const std::string file) {
      uint32_t num = fastaq_objects.size();
      double average = 0;
      uint32_t max = 0;
      uint32_t min = std::numeric_limits<uint32_t>::max();
      for (uint32_t i = 0; i < num; i++) {
        average += fastaq_objects[i]->sequence.size();
        if (fastaq_objects[i]->sequence.size() > max) {
          max = fastaq_objects[i]->sequence.size();
        }
        if (fastaq_objects[i]->sequence.size() < min) {
          min = fastaq_objects[i]->sequence.size();
        }
      }
      average /= num;
      fprintf(stderr, "Stats for: %s\n"
                      "  Number of sequences: %u\n"
                      "  Average length:      %g\n"
                      "  Maximum length:      %u\n"
                      "  Minimum length:      %u\n",
                      file.c_str(), num, average, max, min);
  }
};