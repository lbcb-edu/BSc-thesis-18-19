#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <cctype>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <limits>
#include <ctime>
#include <utility>
#include <unordered_map>
#include <queue>
#include <deque>
#include <map>
#include <algorithm>

#include "long_reads_mapper.hpp"
#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "brown_minimizers.hpp"


const std::set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
const std::set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};

typedef std::tuple<unsigned int, unsigned int, bool> minimizer;
typedef std::tuple<unsigned int, unsigned int, bool> triplet_t;

// position, amount
typedef std::pair<unsigned int, unsigned int> minimizer_index_t;
// query position, reference position, relative strand
//typedef std::tuple<unsigned int, unsigned int, bool> minimizer_hit_t;
typedef std::pair<unsigned int, unsigned int> minimizer_hit_t;

//number of hits, position
typedef std::tuple<unsigned int, unsigned int> stat;

//number of LIS, starting position, relative strand
typedef std::tuple<unsigned int, unsigned int, bool> best_match;

//region number, number of hits
typedef std::pair<unsigned int, int> region_hits;


static struct option long_options[] = {
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},
  {"match", required_argument, NULL, 'M'},
  {"mismatch", required_argument, NULL, 'm'},
  {"gap", required_argument, NULL, 'g'},
  {"window_length", required_argument, NULL, 'w'},
  {"kmers", required_argument, NULL, 'k'},
  {"threads", required_argument, NULL, 't'},
  {NULL, no_argument, NULL, 0}
};


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

  static void parse(
    std::vector<std::unique_ptr<FastAQ>> &fastaq_objects,
    const std::string file, const int file_format) {
      if (file_format == 1) {
        auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FastAQ>(file);
        fasta_parser->parse_objects(fastaq_objects, -1);
      } else {
        auto fastq_parser = bioparser::createParser<bioparser::FastqParser, FastAQ>(file);
        fastq_parser->parse_objects(fastaq_objects, -1);
      }
    }

  static void print_statistics(
    const std::vector<std::unique_ptr<FastAQ>> &fastaq_objects,
    const std::string file) {
      int num = fastaq_objects.size();
      double average = 0;
      uint32_t max = 0;
      uint32_t min = std::numeric_limits<int>::max();
      for (int i = 0; i < num; i++) {
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
                      "  Number of sequences: %d\n"
                      "  Average length:      %g\n"
                      "  Maximum length:      %d\n"
                      "  Minimum length:      %d\n",
                      file.c_str(), num, average, max, min);
  }
};


void help(void) {
  printf("long_reads_mapper - for mapping long erroneous fragments reference  genome.\n\n"

         "Usage: long_reads_mapper [OPTIONS] [file1 file2]\n"
         "file1 - FASTA/FASTQ file containing a set of fragments\n"
         "file2 - FASTA file containing reference genome\n\n"

         "Supported file extensions: .fasta\n"
         "                           .fa\n"
         "                           .fastq\n"
         "                           .fq\n"
         "                           .fasta.gz\n"
         "                           .fa.gz\n"
         "                           .fastq.gz\n"
         "                           .fq.gz\n\n"
         "OPTIONS:\n"
         "  -h  or  --help           print help (displayed now) and exit\n"
         "  -v  or  --version        print version info and exit\n"
         "  -M  or  --match          <int>\n"
         "                             default: 4\n"
         "                             match number\n"
         "  -m  or  --mismatch       <int>\n"
         "                             default: -1\n"
         "                             mismatch number\n"
         "  -g  or  --gap            <int>\n"
         "                             default: -2\n"
         "                             gap number\n"
         "  -w  or  --window_length  <int>\n"
         "                             default: 5\n"
         "                             length of window\n"
         "  -k  or  --kmers          <int>\n"
         "                             default: 15\n"
         "                             constraints: largest supported is 16\n"
         "                             number of letters in substrings\n"
         "  -t  or  --threads        <int>\n"
         "                             default: 3\n"
         "                             number of threads\n"
  );
}


void version(void) {
  printf("brown_mapper %d.%d\n",
    long_reads_mapper_VERSION_MAJOR,
    long_reads_mapper_VERSION_MINOR
  );
}


bool contains_extension(const std::string file, const std::set<std::string> &extensions) {
  for (const auto& it: extensions) {
    if (file.size() > it.size()) {
      if (file.compare(file.size()-it.size(), std::string::npos, it) == 0) {
        return true;
      }
    }
  }
  return false;
}
























void prep_ref(std::vector<triplet_t>& t_minimizers, const float f) {
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

  std::vector<triplet_t> temp;
  temp.reserve(t_minimizers.size());

  for (const auto& minimizer : t_minimizers) {
    if (ref_min_frequency[std::get<0>(minimizer)] < cutoff_freq) {
      temp.push_back(minimizer);
    }
  }

  std::swap(t_minimizers, temp);

  std::sort(t_minimizers.begin(), t_minimizers.end(),
      [] (const triplet_t& a, const triplet_t& b) {
        return (std::get<0>(a) < std::get<0>(b));
      });

  t_minimizers.shrink_to_fit();
}


std::unordered_map<unsigned int, minimizer_index_t> index_ref(
    const std::vector<triplet_t>& t_minimizers) {

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

void find_minimizer_hits(
    const std::unordered_map<unsigned int, minimizer_index_t>& ref_index,
    const std::vector<triplet_t>& t_minimizers,
    const std::vector<triplet_t>& q_minimizers,
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


bool hit_ordering_top_3(const region_hits& a, const region_hits& b) {
  return std::get<1>(a) > std::get<1>(b);
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
          std::sort(result.begin(), result.end(), hit_ordering_top_3);
        }
      }
    }
    std::sort(result.begin(), result.end(), hit_ordering_top_3);
    return result;
}

std::string map_paf(const std::vector<triplet_t>& t_minimizers,
    const std::unordered_map<unsigned int, minimizer_index_t>& ref_index,
    const std::vector<std::unique_ptr<FastAQ>>& fastaq_objects1,
    const std::vector<std::unique_ptr<FastAQ>>& fastaq_objects2,
    unsigned int k, unsigned int window_length,
    int match, int mismatch, int gap,
    unsigned int t_begin, unsigned int t_end, int k_value) {


  for (unsigned int fobj = t_begin; fobj < t_end; fobj++) {

    unsigned int sequence_length = fastaq_objects1[fobj]->sequence.size();

    if(sequence_length < (unsigned) (2 * k_value)) continue;

    std::vector<minimizer> start_minimizers = brown::minimizers(fastaq_objects1[fobj]->sequence.c_str(), k_value, k, window_length);
    std::vector<minimizer> end_minimizers = brown::minimizers(fastaq_objects1[fobj]->sequence.substr(sequence_length - k_value - 1, k_value).c_str(), k_value, k, window_length);

    std::vector<minimizer_hit_t> start_hits_same;
    std::vector<minimizer_hit_t> start_hits_rev;
    std::vector<minimizer_hit_t> end_hits_same;
    std::vector<minimizer_hit_t> end_hits_rev;

    find_minimizer_hits(ref_index, t_minimizers, start_minimizers, start_hits_same, start_hits_rev);
    find_minimizer_hits(ref_index, t_minimizers, end_minimizers, end_hits_same, end_hits_rev);

    std::unordered_map<unsigned int, int> start_hits_region = parse_hits_map(start_hits_same, k_value);
    std::unordered_map<unsigned int, int> start_hits_region_rev = parse_hits_map(start_hits_rev, k_value);
    std::unordered_map<unsigned int, int> end_hits_region = parse_hits_map(end_hits_same, k_value);
    std::unordered_map<unsigned int, int> end_hits_region_rev = parse_hits_map(end_hits_rev, k_value);

    // opcenite informacije
    // printf("Sekvenca broj %d : \n", fobj);
    // printf("Velicina sekvence:  %u\n", sequence_length);
    // printf("Start minimizers size : %lu\n", start_minimizers.size());
    // printf("End minimizers size : %lu\n", end_minimizers.size());
    // printf("start_hits_same hitova : %lu\n", start_hits_same.size());
    // printf("start_hits_rev hitova : %lu\n", start_hits_rev.size());
    // printf("end_hits_same hitova : %lu\n", end_hits_same.size());
    // printf("end_hits_rev hitova : %lu\n", end_hits_rev.size());

    // informacije o regijama
    // printf("Pocetak :\n");
    // for (auto& it: start_hits_region) {
    //   if(it.second < 5) continue;
    //   printf("Regija %d ima %d hitova\n", it.first, it.second);
    // }
    // printf("rev----------------\n");
    // for (auto& it: start_hits_region_rev) {
    //   if(it.second < 5) continue;
    //   printf("Regija %d ima %d hitova\n", it.first, it.second);
    // }
    // printf("Kraj: \n");
    // for (auto& it: end_hits_region) {
    //   if(it.second < 5) continue;
    //   printf("Regija %d ima %d hitova\n", it.first, it.second);
    // }
    // printf("rev----------------\n");
    // for (auto& it: end_hits_region_rev) {
    //   if(it.second < 5) continue;
    //   printf("Regija %d ima %d hitova\n", it.first, it.second);
    // }
    

    int region_size = sequence_length / k_value;
    int threshold = 10;

    std::vector<region_hits> start_hits_top = find_top_3(start_hits_region, threshold);
    std::vector<region_hits> start_hits_top_rev = find_top_3(start_hits_region_rev, threshold);
    std::vector<region_hits> end_hits_top = find_top_3(end_hits_region, threshold);
    std::vector<region_hits> end_hits_top_rev = find_top_3(end_hits_region_rev, threshold);

    // informacije o 3 top regije
    // printf("Pocetak :\n");
    // for(auto it : start_hits_top){
    //   printf("Regija %d ima %d hitova\n", it.first, it.second);
    // }
    // for(auto it : start_hits_top_rev){
    //   printf("Regija %d ima %d hitova, reverse\n", it.first, it.second);
    // }
    // printf("Kraj: \n");
    // for(auto it : end_hits_top){
    //   printf("Regija %d ima %d hitova\n", it.first, it.second);
    // }
    // for(auto it : end_hits_top_rev){
    //   printf("Regija %d ima %d hitova, reverse\n", it.first, it.second);
    // }

    int hits_number = 0;
    unsigned int starting_region = 0;
    unsigned int ending_region = 0;
    bool rev = false;

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

    if(hits_number == 0){
      printf("%u PRESKACEM\n", sequence_length);
      //ispis top 3 regije za pocetak i kraj
      // printf("Pocetak :\n");
      // for(auto it : start_hits_top){
      //   printf("Regija %d ima %d hitova\n", it.first, it.second);
      // }
      // for(auto it : start_hits_top_rev){
      //   printf("Regija %d ima %d hitova, reverse\n", it.first, it.second);
      // }
      // printf("Kraj: \n");
      // for(auto it : end_hits_top){
      //   printf("Regija %d ima %d hitova\n", it.first, it.second);
      // }
      // for(auto it : end_hits_top_rev){
      //   printf("Regija %d ima %d hitova, reverse\n", it.first, it.second);
      // }
      // printf("\n");
      continue;
    }

    // test za + 
    // if(!rev){
    //   int prvih10=0;
    //   std::vector<minimizer_hit_t> v;
    //   for(auto& minimizer : start_hits_same){
    //    if(minimizer.second / k_value != starting_region && minimizer.second / k_value -1 != starting_region) continue;
    //     v.push_back(minimizer);
    //   }
    //   std::sort(v.begin(), v.end(), 
    //     [] (const minimizer_hit_t& a, const minimizer_hit_t& b) {
    //     if (std::get<0>(a) == std::get<0>(b)){
    //       return std::get<1>(a) < std::get<1>(b);
    //     }
    //     return std::get<0>(a) < std::get<0>(b);});

    //   for(auto& minimizer : v){
    //     if(prvih10 > 10) break;
    //     printf("%d , %d\n", minimizer.first, minimizer.second);
    //     prvih10++;
    //   }
    //   printf("\n");
    //   prvih10= 0;
    //   std::vector<minimizer_hit_t> v1;
    //   for(auto& minimizer : end_hits_same){
    //    if(minimizer.second / k_value != ending_region && minimizer.second / k_value -1 != ending_region) continue;
    //     v1.push_back(minimizer);
    //   }
    //   std::sort(v1.begin(), v1.end(),
    //     [] (const minimizer_hit_t& a, const minimizer_hit_t& b) {
    //     if (std::get<0>(a) == std::get<0>(b)){
    //       return std::get<1>(a) < std::get<1>(b);
    //     }
    //     return std::get<0>(a) > std::get<0>(b);});

    //   for(auto& minimizer : v1){
    //     if(prvih10 > 10) break;
    //     printf("%d = %u, %d\n", minimizer.first, sequence_length - k_value - 1 + minimizer.first, minimizer.second);
    //     prvih10++;
    //   }
    // }

    unsigned int start = k_value;
    unsigned int end = 0;
    unsigned int ref_start = 0;
    unsigned int ref_end = 0;

    if(!rev){
      for(auto& minimizer : start_hits_same){
        if(minimizer.second / k_value != starting_region && minimizer.second / k_value -1 != starting_region) continue;
        if(minimizer.first < start){
          start = minimizer.first;
          ref_start = minimizer.second;
        } 
      }
      for(auto& minimizer : end_hits_same){
        if(minimizer.second / k_value != ending_region && minimizer.second / k_value -1 != ending_region) continue;
        if(minimizer.first > end){
          end = minimizer.first;
          ref_end = minimizer.second;
        } 
      }
    }else{
      for(auto& minimizer : start_hits_rev){
        if(minimizer.second / k_value != starting_region && minimizer.second / k_value -1 != starting_region) continue;
        if(minimizer.first < start){
          start = minimizer.first;
          ref_end = minimizer.second;
        } 
      }
      for(auto& minimizer : end_hits_rev){
        if(minimizer.second / k_value != ending_region && minimizer.second / k_value -1 != ending_region) continue;
        if(minimizer.first > end){
          end = minimizer.first;
          ref_start = minimizer.second;
        }
      }
    }

    end = sequence_length - k_value - 1 + end + k;
    ref_end += k;

    printf("%u\t%u\t%u\t%c\t%u\t%u\n", sequence_length, start, end, rev ? '-' : '+', ref_start, ref_end);
    // printf("\n");
  }
  return "";
}




































int main (int argc, char **argv) {
  int optchr;

  int match = 4;
  int mismatch = -1;
  int gap = -2;

  int window_length = 5;
  int k = 15;
  int t = 2;
  float f = 0.001f;

  while ((optchr = getopt_long(argc, argv, "hvm:g:M:k:w:t:", long_options, NULL)) != -1) {
    switch (optchr) {
      case 'h': {
        help();
        exit(0);
      }
      case 'v': {
        version();
        exit(0);
      }
      case 'M': {
        match = atoi(optarg);
        break;
      }
      case 'm': {
        mismatch = atoi(optarg);
        break;
      }
      case 'g': {
        gap = atoi(optarg);
        break;
      }
      case 'w': {
        window_length = atoi(optarg);
        break;
      }
      case 'k': {
        k = atoi(optarg);
        break;
      }
      case 't': {
        t = atoi(optarg);
        if (t < 1) {
          fprintf(stderr, "[mapper] error: t must be positive!\n"); 
          exit(1); 
        }
        break;
      }
      default: {
        fprintf(stderr, "[mapper] error: Unknown option. Type %s --help for usage.\n", argv[0]);
        exit(1);
      }
    }
  }

  if (argc - optind != 2) {
    fprintf(stderr, "[mapper] error: Expected 2 mapping arguments! Use --help for usage.\n");
    exit(1);
  }

  // Load files
  fprintf(stderr, "\nLoading files...");

  std::string file1 (argv[optind]);
  std::string file2 (argv[optind+1]);

  int file1_format = 0;
  int file2_format = 0;

  if (contains_extension(file1, fasta_formats)) {
    file1_format = 1;
  } else if (contains_extension(file1, fastq_formats)) {
    file1_format = 2;
  }

  file2_format = contains_extension(file2, fasta_formats);

  if ( !(file1_format && file2_format) ) {
    fprintf(stderr, "[mapper] error: Unsupported format(s)! Check --help for supported file formats.\n");
    exit(1);
  }

  fprintf(stderr, " Done!\n\nParsing files...");

  // Parse files
  std::vector<std::unique_ptr<FastAQ>> fastaq_objects1;
  std::vector<std::unique_ptr<FastAQ>> fastaq_objects2;

  FastAQ::parse(fastaq_objects1, file1, file1_format);
  FastAQ::parse(fastaq_objects2, file2, file2_format);

  fprintf(stderr, " Done!\n");

  // Print file stats
  fprintf(stderr, "\n");  
  // FastAQ::print_statistics(fastaq_objects1, file1);
  // FastAQ::print_statistics(fastaq_objects2, file2);

  // Regular version of reference minimizers search
  // std::vector<triplet_t> t_minimizers = brown::minimizers(
  //     fastaq_objects2[0]->sequence.c_str(), fastaq_objects2[0]->sequence.size(), k, window_length);


  //kod
  // int size = fastaq_objects1[0]->sequence.size();
  // printf("\nVelicina reference : %d\n", size);

  int k_value = 1000;

  











  std::vector<triplet_t> t_minimizers = brown::minimizers(
      fastaq_objects2[0]->sequence.c_str(), fastaq_objects2[0]->sequence.size(), k, window_length);

  prep_ref(t_minimizers, f);

  std::unordered_map<unsigned int, minimizer_index_t> ref_index = index_ref(t_minimizers);

  //prvih 250 - zbog testiranja
  map_paf(std::ref(t_minimizers), std::ref(ref_index), std::ref(fastaq_objects1), std::ref(fastaq_objects2),
      k, window_length, match, mismatch, gap, 0, 250, k_value);


  // kod za visedretvenost
  // std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool(t);

  // std::vector<std::future<std::string>> thread_futures;

  // for (unsigned int tasks = 0; tasks < t - 1; ++tasks) {
  //   thread_futures.emplace_back(thread_pool->submit_task(map_paf, std::ref(t_minimizers), std::ref(ref_index),
  //     std::ref(fastaq_objects1), std::ref(fastaq_objects2),
  //     k, window_length, alignment, match, mismatch, gap, tasks * fastaq_objects1.size() / t, (tasks + 1) * fastaq_objects1.size() / t, bw, c, k_value));  
  // }
  // thread_futures.emplace_back(thread_pool->submit_task(map_paf, std::ref(t_minimizers), std::ref(ref_index),
  //     std::ref(fastaq_objects1), std::ref(fastaq_objects2),
  //     k, window_length, alignment, match, mismatch, gap, (t - 1) * fastaq_objects1.size() / t, fastaq_objects1.size(), bw, c, k_value));
  
  // for (auto& it : thread_futures) {
  //   it.wait();
  //   std::cout << it.get();
  // }

  return 0;
}
