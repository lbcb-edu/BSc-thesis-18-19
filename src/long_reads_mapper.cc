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
typedef std::tuple<unsigned int, unsigned int, bool> minimizer_hit_t;
//number of hits, position
typedef std::tuple<unsigned int, unsigned int> stat;

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































bool hit_ordering(const minimizer_hit_t& a, const minimizer_hit_t& b) {
  if (std::get<2>(a) == std::get<2>(b)) {
    if(std::get<0>(a) == std::get<0>(b)){
      return std::get<1>(a) < std::get<1>(b);
    }
    return std::get<0>(a) < std::get<0>(b);
  }
  return std::get<2>(a) < std::get<2>(b);
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


// std::unordered_map<unsigned int, std::vector<minimizer_hit_t>> find_minimizer_hits(
//     const std::unordered_map<unsigned int, minimizer_index_t>& ref_index,
//     const std::vector<triplet_t>& t_minimizers,
//     const std::vector<triplet_t>& q_minimizers,
//     int k_value) {

//   std::unordered_map<unsigned int, minimizer_hit_t> hits;

//   for (const auto& minimizer : q_minimizers) {
//     auto found = ref_index.find(std::get<0>(minimizer));
//     if (found == ref_index.end()) continue;

//     for (unsigned int i = 0; i < found->second.second; i++){
//       bool srand = false;
//       if (std::get<2>(minimizer) == std::get<2>(t_minimizers[found->second.first + i])) bool = true;
//       int position = std::get<1>(t_minimizers[found->second.first + i]);

//       if(position < k_value){
//         hits[0].emplace_back(std::get<1>(minimizer),
//           std::get<1>(t_minimizers[found->second.first + i]), srand);
//       }else{
//         int key = position / k_value;
//         hits[key].emplace_back(std::get<1>(minimizer),
//           std::get<1>(t_minimizers[found->second.first + i]), srand);
//         hits[key - 1].emplace_back(std::get<1>(minimizer),
//           std::get<1>(t_minimizers[found->second.first + i]), srand);
//       }
//     }
//   }
//   return hits;
// }

std::vector<minimizer_hit_t> find_minimizer_hits(
    const std::unordered_map<unsigned int, minimizer_index_t>& ref_index,
    const std::vector<triplet_t>& t_minimizers,
    const std::vector<triplet_t>& q_minimizers) {

  std::vector<minimizer_hit_t> hits;

  for (const auto& minimizer : q_minimizers) {
    auto found = ref_index.find(std::get<0>(minimizer));
    if (found != ref_index.end()) {
      for (unsigned int i = 0; i < found->second.second; i++) {
        if (std::get<2>(minimizer) == std::get<2>(t_minimizers[found->second.first + i])) {
          hits.emplace_back(
            std::get<1>(minimizer),
            std::get<1>(t_minimizers[found->second.first + i]),
            0);  
        } else {
          hits.emplace_back(
            std::get<1>(minimizer),
            std::get<1>(t_minimizers[found->second.first + i]),
            1);
        }
        
      }
    }
  }

  return hits;
}

std::string map_paf(const std::vector<triplet_t>& t_minimizers,
    const std::unordered_map<unsigned int, minimizer_index_t>& ref_index,
    const std::vector<std::unique_ptr<FastAQ>>& fastaq_objects1,
    const std::vector<std::unique_ptr<FastAQ>>& fastaq_objects2,
    unsigned int k, unsigned int window_length, brown::AlignmentType alignment,
    int match, int mismatch, int gap,
    unsigned int t_begin, unsigned int t_end, int bw, bool c, int k_value) {

  std::string paf;


  for (unsigned int fobj = t_begin; fobj < t_end; fobj++) {

    if(fastaq_objects1[fobj]->sequence.size() < 2 * k_value) continue;

    std::vector<triplet_t> start_minimizers = brown::minimizers(
      fastaq_objects1[fobj]->sequence.c_str(), fastaq_objects1[fobj]->sequence.size(), k, window_length);

    std::vector<minimizer> start_minimizers = brown::minimizers(fastaq_objects1[fobj]->sequence.c_str(), k_value, k, window_length);
    std::vector<minimizer> end_minimizers = brown::minimizers(fastaq_objects1[fobj]->sequence.substr(sequence_length - k_value - 1, k_value).c_str(), k_value, k, window_length);

    // std::unordered_map<unsigned int, minimizer_hit_t> start_hits = find_minimizer_hits(ref_index, t_minimizers, q_minimizers, k_value);
    // std::unordered_map<unsigned int, minimizer_hit_t> end_hits = find_minimizer_hits(ref_index, t_minimizers, q_minimizers, k_value);

    std::vector<minimizer_hit_t> start_hits = find_minimizer_hits(ref_index, t_minimizers, start_minimizers);
    std::vector<minimizer_hit_t> end_hits = find_minimizer_hits(ref_index, t_minimizers, end_minimizers);

    std::sort(start_hits.begin(), start_hits.end(), hit_ordering);
    std::sort(end_hits.begin(), end_hits.end(), hit_ordering);









    std::vector<std::pair<region_t, unsigned int>> regions;

    unsigned int b = 0;
    for (unsigned int e = 0; e < hits.size(); ++e) {
      if (e == hits.size() - 1 ||
          std::get<2>(hits[e + 1]) != std::get<2>(hits[e]) ||
          std::abs(band(hits[e + 1], hits[e])) >= bw) {
        std::pair<region_t, unsigned int> temp_region = find_region(hits, b, e, k);
        if (temp_region.second > 4 * k) {
          regions.push_back(temp_region);
        }
        b = e + 1;
      }
    }

    if (regions.empty()) {
      continue;
    }

    merge_regions(regions, bw);

    std::sort(regions.begin(), regions.end(),
        [] (const std::pair<region_t, unsigned int>& a,
            const std::pair<region_t, unsigned int>& b) {
      return a.second > b.second;
    });

    for (const auto& reg : regions) {
      region_t max_region = reg.first;
      unsigned int matches = reg.second;

      std::string rel_strand;
      unsigned int total;
      unsigned int q_b, q_e, t_b, t_e;
      std::string cigar;
      unsigned int target_begin;

      t_b = std::get<1>(max_region.first);
      t_e = std::get<1>(max_region.second) + k;
      q_b = std::get<0>(max_region.first);
      q_e = std::get<0>(max_region.second) + k;
      total = (t_e - t_b);
      
      if (!std::get<2>(max_region.first)) {
        rel_strand = "+";
        if (c) {
          brown::pairwise_alignment(fastaq_objects1[fobj]->sequence.c_str() + q_b,
              q_e - q_b,
              fastaq_objects2[0]->sequence.c_str() + t_b,
              t_e - t_b,
              alignment, match, mismatch, gap, cigar, target_begin);
          matches = read_cigar(cigar, total);
        }
      } else {
        rel_strand = "-";
        if (c) {
          std::string rc = reverse_complement(fastaq_objects1[fobj]->sequence, q_b,
              q_e - q_b);
          brown::pairwise_alignment(rc.c_str(), rc.size(),
              fastaq_objects2[0]->sequence.c_str() + t_b,
              t_e - t_b,
              alignment, match, mismatch, gap, cigar, target_begin);
          matches = read_cigar(cigar, total);
        }
      }
      paf += fastaq_objects1[fobj]->name + "\t" +
             std::to_string(fastaq_objects1[fobj]->sequence.size()) + "\t" +
             std::to_string(q_b) + "\t" +
             std::to_string(q_e) + "\t" +
             rel_strand + "\t" +
             fastaq_objects2[0]->name + "\t" +
             std::to_string(fastaq_objects2[0]->sequence.size()) + "\t" +
             std::to_string(t_b) + "\t" +
             std::to_string(t_e) + "\t" +
             std::to_string(matches) + "\t" +
             std::to_string(total) + "\t" +
             std::to_string(255);
      if (c) {
        paf += "\tcg:Z:" + cigar;
      }
      paf += "\n";
    }
  }
  return paf;
}





































int main (int argc, char **argv) {
  int optchr;

  int match = 4;
  int mismatch = -1;
  int gap = -2;

  int window_length = 5;
  int k = 15;
  int t = 2;

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
  FastAQ::print_statistics(fastaq_objects1, file1);
  FastAQ::print_statistics(fastaq_objects2, file2);

  // Regular version of reference minimizers search
  // std::vector<triplet_t> t_minimizers = brown::minimizers(
  //     fastaq_objects2[0]->sequence.c_str(), fastaq_objects2[0]->sequence.size(), k, window_length);


  //kod
  int size = fastaq_objects1[0]->sequence.size();
  printf("\nVelicina reference : %d\n", size);

  int k_value = 1000;

  











  std::vector<triplet_t> t_minimizers = brown::minimizers(
      fastaq_objects2[0]->sequence.c_str(), fastaq_objects2[0]->sequence.size(), k, window_length);

  prep_ref(t_minimizers, f);

  std::unordered_map<unsigned int, minimizer_index_t> ref_index = index_ref(t_minimizers);

  std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool(t);

  std::vector<std::future<std::string>> thread_futures;

  for (unsigned int tasks = 0; tasks < t - 1; ++tasks) {
    thread_futures.emplace_back(thread_pool->submit_task(map_paf, std::ref(t_minimizers), std::ref(ref_index),
      std::ref(fastaq_objects1), std::ref(fastaq_objects2),
      k, window_length, alignment, match, mismatch, gap, tasks * fastaq_objects1.size() / t, (tasks + 1) * fastaq_objects1.size() / t, bw, c, k_value));  
  }
  thread_futures.emplace_back(thread_pool->submit_task(map_paf, std::ref(t_minimizers), std::ref(ref_index),
      std::ref(fastaq_objects1), std::ref(fastaq_objects2),
      k, window_length, alignment, match, mismatch, gap, (t - 1) * fastaq_objects1.size() / t, fastaq_objects1.size(), bw, c, k_value));
  
  for (auto& it : thread_futures) {
    it.wait();
    std::cout << it.get();
  }












  // std::multimap<unsigned int, int> rm;

  // for(int i = 0; i < size ; i = i + k_value){
  //   int length = 2 * k_value;
  //   if(i + length > size) length = size - i;
  //   std::string tmp = fastaq_objects1[0]->sequence.substr(i, length);
  //   std::vector<minimizer> minimizers = brown::minimizers(tmp.c_str(), length, k, window_length);


  //   for(auto& minimizer : minimizers){
  //     rm.insert({std::get<0>(minimizer), i / 1000});
  //   }
  // }

  // printf("Velicina rm mape : %d\n\n", rm.size());

  // for(int i = 0; i <  100 /*fastaq_objects2.size()*/; i++){
  //   int sequence_length = fastaq_objects2[i]->sequence.size();

  //   printf("Sekvenca br. %d, duljina sekvence : %d\n", i, sequence_length);
    
  //   if(sequence_length < 2 * k_value){
  //     printf("Preskacem\n");
  //     printf("----------------------------\n\n");
  //     continue;
  //   }

  //   std::vector<minimizer> start_minimizers = brown::minimizers(fastaq_objects2[i]->sequence.c_str(), k_value, k, window_length);
  //   std::vector<minimizer> end_minimizers = brown::minimizers(fastaq_objects2[i]->sequence.substr(sequence_length - k_value - 1, k_value).c_str(), k_value, k, window_length);

  //   printf("Broj start_minimizers: %d , broj end_minimizers: %d \n", start_minimizers.size(), end_minimizers.size());
  //   printf("Pocetak:\n");

  //   //kopija koda
  //   std::unordered_map<int, int> best;
  //   for(auto& minimizer : start_minimizers){

  //     auto itr1 = rm.lower_bound(std::get<0>(minimizer)); 
  //     auto itr2 = rm.upper_bound(std::get<0>(minimizer)); 
      
  //     while(itr1 != itr2){
  //       best[itr1->second]++;
  //       itr1++;
  //     } 
  //   }

  //   int max = 0;
  //   int pos = 0;
  //   for(auto& pair : best){
  //     if(pair.second > max){
  //       pos = pair.first;
  //       max = pair.second;
  //     }
  //     if(pair.second > 40) printf("Regija %d ima %d pogodaka\n", pair.first, pair.second);
  //   }
  //   printf("Najbolja je %d regija, s %d pogodaka\n", pos, max);

  //   printf("Kraj:\n");
    
  //   //kopija koda
  //   std::unordered_map<int, int> best2;
  //   for(auto& minimizer : end_minimizers){

  //     auto itr1 = rm.lower_bound(std::get<0>(minimizer)); 
  //     auto itr2 = rm.upper_bound(std::get<0>(minimizer)); 
      
  //     while(itr1 != itr2){
  //       best2[itr1->second]++;
  //       itr1++;
  //     } 
  //   }

  //   max = 0;
  //   pos = 0;
  //   for(auto& pair : best2){
  //     if(pair.second > max){
  //       pos = pair.first;
  //       max = pair.second;
  //     }
  //     if(pair.second > 40) printf("Regija %d ima %d pogodaka\n", pair.first, pair.second);
  //   }
  //   printf("Najbolja je %d regija, s %d pogodaka\n", pos, max);
  //   printf("----------------------------\n\n");
  // }


  return 0;
}
