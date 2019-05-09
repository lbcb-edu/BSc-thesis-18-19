#pragma once

#include <cstdint>

// Mapping parameters wrapper
struct mapping_params {
  int32_t mch;
  int32_t mis;
  int32_t gapo;
  int32_t gape;
  uint32_t band;
  uint32_t k;
  uint32_t w;
  float f;
  uint32_t insert_size;
  uint32_t region_size;
  uint32_t threshold;
};