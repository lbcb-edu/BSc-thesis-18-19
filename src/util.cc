#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <limits>
#include <utility>
#include <unordered_map>
#include <algorithm>
#include <cstring>

#include "util.hpp"


void find_minimizer_hits(
    const std::unordered_map<unsigned int, minimizer_index_t>& ref_index,
    const std::vector<minimizer>& t_minimizers,
    const std::vector<minimizer>& q_minimizers,
    std::vector<minimizer_hit_t>& same,
    std::vector<minimizer_hit_t>& rev) {

  for (const auto& minimizer : q_minimizers) {
    auto found = ref_index.find(std::get<0>(minimizer));
    if (found != ref_index.end()) {
      for (unsigned int i = 0; i < found->second.second; i++) {
        if (std::get<2>(minimizer) == std::get<2>(t_minimizers[found->second.first + i])) {
          same.emplace_back(
            std::get<1>(minimizer),
            std::get<1>(t_minimizers[found->second.first + i]));  
        } else {
          rev.emplace_back(
            std::get<1>(minimizer),
            std::get<1>(t_minimizers[found->second.first + i]));
        }
      }
    }
  }
}

std::unordered_map<unsigned int, int> parse_hits_map(std::vector<minimizer_hit_t> & hits, int k_value){
  std::unordered_map<unsigned int, int> result;
  for(auto& hit : hits){
    unsigned int region = std::get<1>(hit) / k_value;
    result[region]++;
    if(region != 0) result[region - 1]++;
  }
  return result;
}

std::vector<region_hits> find_top_3(std::unordered_map<unsigned int, int>& map, int threshold){
  std::vector<region_hits> result;
    for (auto& it: map) {
      if(it.second < threshold) continue;
      if(result.size() < 3){
        result.push_back(it);
      }else{
        if(it.second > std::get<1>(result[2])){
          result.pop_back();
          result.push_back(it);
          std::sort(result.begin(), result.end(),
            [] (const region_hits& a, const region_hits& b) {
              return (std::get<1>(a) > std::get<1>(b));
          });
        }
      }
    }
    std::sort(result.begin(), result.end(),
            [] (const region_hits& a, const region_hits& b) {
              return (std::get<1>(a) > std::get<1>(b));
    });
    return result;
}

std::vector<region_hits> find_top_all(std::unordered_map<unsigned int, int>& map){
  std::vector<region_hits> result;
    for (auto& it: map) {
      result.push_back(it);
    }
    std::sort(result.begin(), result.end(),
            [] (const region_hits& a, const region_hits& b) {
              return (std::get<1>(a) > std::get<1>(b));
    });
    return result;
}

void find_regions(
	std::vector<region_hits>& start_hits_top,
	std::vector<region_hits>& start_hits_top_rev,
  	std::vector<region_hits>& end_hits_top,
  	std::vector<region_hits>& end_hits_top_rev,
  	int& hits_number,
  	unsigned int& starting_region,
  	unsigned int& ending_region,
  	bool& rev,
  	int region_size){

  for(auto& start : start_hits_top){
    for(auto& end : end_hits_top){
      if(end.first == (start.first + region_size) || end.first == (start.first + region_size - 1)){
        int sum = start.second + end.second;
        if(sum > hits_number){
          hits_number = sum;
          starting_region = start.first;
          ending_region = end.first;
        }
      }
    }
  }

  for(auto start : start_hits_top_rev){
    for(auto end : end_hits_top_rev){
      if(start.first == (end.first + region_size) || start.first == (end.first + region_size - 1)){
        int sum = start.second + end.second;
        if(sum > hits_number){
          hits_number = sum;
          starting_region = start.first;
          ending_region = end.first;
          rev = true;
        }
      }
    }
  }
}

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
	int k_value){

  if(!rev){
    for(auto& minimizer : start_hits_same){
      if((minimizer.second / k_value) == starting_region || (minimizer.second / k_value) - 1 == starting_region) {
        if(minimizer.first < start){
          start = minimizer.first;
          ref_start = minimizer.second;
        } 
      }
    }
    for(auto& minimizer : end_hits_same){
      if((minimizer.second / k_value) == ending_region || (minimizer.second / k_value) - 1 == ending_region) {
        if(minimizer.first > end){
          end = minimizer.first;
          ref_end = minimizer.second;
        } 
      }
    }
  }else{
    for(auto& minimizer : start_hits_rev){
      if((minimizer.second / k_value) == starting_region || (minimizer.second / k_value) - 1 == starting_region){
        if(minimizer.first < start){
          start = minimizer.first;
          ref_end = minimizer.second;
        }
      }
    }
    for(auto& minimizer : end_hits_rev){
      if((minimizer.second / k_value) == ending_region || (minimizer.second / k_value) - 1 == ending_region){
        if(minimizer.first > end){
          end = minimizer.first;
          ref_start = minimizer.second;
        }
      }
    }
  }
}


void find_what_exists(
	std::vector<region_hits>& start_hits_top, 
	std::vector<region_hits>& start_hits_top_rev,
  	std::vector<region_hits>& end_hits_top, 
  	std::vector<region_hits>& end_hits_top_rev, 
  	bool& bool_start, bool& bool_end){

  bool bool_start_same = start_hits_top.size() > 0 ? true : false;
  bool bool_end_same = end_hits_top.size() > 0 ? true : false;
  bool bool_start_rev = start_hits_top_rev.size() > 0 ? true : false;
  bool bool_end_rev = end_hits_top_rev.size() > 0 ? true : false;

  bool_start = bool_start_same || bool_start_rev;
  bool_end = bool_end_same || bool_end_rev;

  if(bool_start && bool_end){
    int max_start;
    int max_end;

    if(bool_start_same && bool_start_rev){
      max_start = std::max(start_hits_top[0].second, start_hits_top_rev[0].second);
    }else if(bool_start_same){
      max_start = start_hits_top[0].second;
    }else{
      max_start = start_hits_top_rev[0].second;
    }

    if(bool_end_same && bool_end_rev){
      max_end = std::max(end_hits_top[0].second, end_hits_top_rev[0].second);
    }else if(bool_end_same){
      max_end = end_hits_top[0].second;
    }else{
      max_end = end_hits_top_rev[0].second;
    }

    if(max_start > max_end){
      bool_end = false;
    }else{
      bool_start = false;
    }
  }   
}