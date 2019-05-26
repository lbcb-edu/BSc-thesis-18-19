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
#include <cstring>


#include "long_reads_mapper.hpp"
#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "brown_minimizers.hpp"
#include "ksw2.h"


#include "fastaq.hpp"
#include "index.hpp"
#include "util.hpp"

const std::set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
const std::set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};

std::map<char, char> complement_map = {{'C', 'G'}, {'A', 'T'}, {'T', 'A'}, {'U', 'A'}, {'G', 'C'}};

typedef std::tuple<unsigned int, unsigned int, bool> minimizer;
typedef std::tuple<unsigned int, unsigned int, bool> triplet_t;

// position, amount
typedef std::pair<unsigned int, unsigned int> minimizer_index_t;

// query position, reference position
typedef std::pair<unsigned int, unsigned int> minimizer_hit_t;

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
         "  -K  or  --kvalue         <int>\n"
         "                             default: 1000\n"
         "                             number of bases to take from start and end\n"
         "  -M  or  --match          <int>\n"
         "                             default: 2\n"
         "                             match number\n"
         "  -m  or  --mismatch       <int>\n"
         "                             default: -4\n"
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
         "  -c  or  --cigar          enable cigar string\n"
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

std::string reverse_complement(const std::string& original, unsigned int pos, unsigned int length) {
  std::string rc(original.begin() + pos, original.begin() + pos + length);
  unsigned int j = pos + length - 1;
  for (unsigned int i = 0; i < length; ++i) {
    rc[i] = complement_map[original[j--]];
  }
  return rc;
}


std::string map_paf(const std::vector<triplet_t>& t_minimizers,
    const std::unordered_map<unsigned int, minimizer_index_t>& ref_index,
    const std::vector<std::unique_ptr<FastAQ>>& data_set,
    const std::vector<std::unique_ptr<FastAQ>>& reference,
    unsigned int k, unsigned int window_length,
    int match, int mismatch, int gap,
    unsigned int t_begin, unsigned int t_end, int k_value,
    bool bool_cigar, unsigned int ref_num) {

  std::string paf;
  std::string sam;
  double percentage = 0.5;
  int threshold = 5;

  for (unsigned int fobj = t_begin; fobj < t_end; fobj++) {
    unsigned int sequence_length = data_set[fobj]->sequence.size();

    if(sequence_length < (unsigned) (2 * k_value)) continue; //ovo netreba jer su file-ovi vec filtrirani

    std::vector<minimizer> start_minimizers = brown::minimizers(data_set[fobj]->sequence.c_str(), k_value, k, window_length);
    std::vector<minimizer> end_minimizers = brown::minimizers(data_set[fobj]->sequence.substr(sequence_length - k_value - 1, k_value).c_str(), k_value, k, window_length);

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

    unsigned int region_size = sequence_length / k_value;    

    std::vector<region_hits> start_hits_top = find_top_3(start_hits_region, threshold);
    std::vector<region_hits> start_hits_top_rev = find_top_3(start_hits_region_rev, threshold);
    std::vector<region_hits> end_hits_top = find_top_3(end_hits_region, threshold);
    std::vector<region_hits> end_hits_top_rev = find_top_3(end_hits_region_rev, threshold);

    int hits_number = 0;
    unsigned int starting_region = 0;
    unsigned int ending_region = 0;
    bool rev = false;

    find_regions(start_hits_top, start_hits_top_rev, end_hits_top, end_hits_top_rev, hits_number, starting_region, ending_region, rev, region_size);

    bool bool_start = false;
    bool bool_end = false;

    bool ispis = false; //izbrisi -> koristi se zbog testiranja

    if(hits_number == 0){

      find_what_exists(start_hits_top, start_hits_top_rev, end_hits_top, end_hits_top_rev, bool_start, bool_end);

      //ako nema pocetka ni kraja preskoci
      if(!bool_start && !bool_end){
        // printf("%s PRESKACEM, nema niceg\n", data_set[fobj]->name.c_str());
        continue;
      }

      //ima pocetak
      if(bool_start){

        bool bool_found_something = false;
        int stop = region_size * percentage;

        for(int i = 1; i <= stop + 1; i++){
          end_hits_same.clear();
          end_hits_rev.clear();
          end_minimizers = brown::minimizers(data_set[fobj]->sequence.substr(sequence_length - (i * k_value) - 1, k_value).c_str(), k_value, k, window_length);
          find_minimizer_hits(ref_index, t_minimizers, end_minimizers, end_hits_same, end_hits_rev);
          end_hits_region = parse_hits_map(end_hits_same, k_value);
          end_hits_region_rev = parse_hits_map(end_hits_rev, k_value);
          end_hits_top = find_top_all(end_hits_region);
          end_hits_top_rev = find_top_all(end_hits_region_rev);
          find_regions(start_hits_top, start_hits_top_rev, end_hits_top, end_hits_top_rev, hits_number, starting_region, ending_region, rev, region_size - (i - 1));

          if (hits_number != 0){
            bool_found_something = true;
            region_size = region_size - (i - 1);
            break;
          }
        }

        //ako nije nista nadeno preskoci
        if(!bool_found_something){
          // printf("%s PRESKACEM, ima pocetak, ali nije nadeno\n", data_set[fobj]->name.c_str());
          continue;
        }
        ispis = true;  //za testiranje
      }

      //ima kraj
      if(bool_end){
        
        bool bool_found_something = false;
        int stop = region_size * percentage;

        for(int i =  0; i <= stop; i++){
          start_hits_same.clear();
          start_hits_rev.clear();
          start_minimizers = brown::minimizers(data_set[fobj]->sequence.substr((i * k_value), k_value).c_str(), k_value, k, window_length);
          find_minimizer_hits(ref_index, t_minimizers, start_minimizers, start_hits_same, start_hits_rev);
          start_hits_region = parse_hits_map(start_hits_same, k_value);
          start_hits_region_rev = parse_hits_map(start_hits_rev, k_value);
          start_hits_top = find_top_all(start_hits_region);
          start_hits_top_rev = find_top_all(start_hits_region_rev);
          find_regions(start_hits_top, start_hits_top_rev, end_hits_top, end_hits_top_rev, hits_number, starting_region, ending_region, rev, region_size - i);

          if (hits_number != 0){
            bool_found_something = true;
            region_size = region_size - i;
            break;
          }
        }

        //ako nije nista nadeno preskoci
        if(!bool_found_something){
          // printf("%s PRESKACEM, ima kraj, ali nije nadeno\n", data_set[fobj]->name.c_str());
          continue;
        }
      }
    }

    unsigned int start = k_value;
    unsigned int end = 0;
    unsigned int ref_start = 0;
    unsigned int ref_end = 0;

    find_positions(start_hits_same, start_hits_rev, end_hits_same,end_hits_rev, rev, starting_region, ending_region, start, end, ref_start, ref_end, k_value);

    unsigned int help_end = end;
    ref_end += k;
    end = sequence_length - k_value - 1 + end + k;

    //ako se trazio samo kraj ili samo pocetak
    if(region_size != sequence_length / k_value){
      if(bool_start){
        int diff = (sequence_length / k_value - region_size) + 1;
        end = sequence_length - k_value * diff - 1 + help_end + k;
      }
      if(bool_end){
        int diff = sequence_length / k_value - region_size;
        start += diff * k_value;
      }
    }

    std::string cigar ("NO CIGAR");

    std::string seq;
    if(rev){
      seq = reverse_complement(data_set[fobj]->sequence, start, end - start);
    }else{
      seq = data_set[fobj]->sequence.substr(start, end - start);
    }

    //ksw2
    if(bool_cigar){
      int8_t a = match;
      int8_t b = mismatch;
      int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
      int tl = ref_end - ref_start;
      int ql = end - start;
      uint8_t *ts, *qs, c[256];
      ksw_extz_t ez;
      memset(&ez, 0, sizeof(ksw_extz_t));
      memset(c, 5, 256);
      // c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
      // c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; c['U'] = c['u'] = 3;

      c['C'] = c['c'] = 0; c['A'] = c['a'] = 1;
      c['T'] = c['t'] = c['U'] = c['u'] = 2; c['G'] = c['g'] = 3;

      ts = (uint8_t*)malloc(tl);
      qs = (uint8_t*)malloc(ql);

      for (int i = 0; i < tl; ++i) ts[i] = c[(uint8_t)seq[i]]; // encode to 0/1/2/3
      for (int i = 0; i < ql; ++i) qs[i] = c[(uint8_t)reference[ref_num]->sequence[i + ref_start]]; // encode to 0/1/2/3      

      int gapo = 4;
      int gapo2 = 24;
      int gape = 1;
      int gape2 = 2;
      int w = 500;
      int z1 = 400;
      int z2 = 200;
      int s = 80;
      int noncan = 0;

      // printf("%d\t%d\n", tl, ql);
      // for(int i = 0; i < 20; i++){
      //   std::cout << data_set[fobj]->sequence[i + start];
      // }
      // printf("\n");
      // for(int i = 0; i < 20; i++){
      //   std::cout << reference[ref_num]->sequence[i + ref_start];
      // }
      // printf("\n");

      ksw_extz2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, w, z1, 0, 0, &ez);
      // ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, gapo2, gape2, w, z1, 0, 0, &ez);
      // ksw_exts2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, gapo2, noncan, z1,  0, &ez); //ispis??

      //ove tri minimap2 ne implementira
      // ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, w, z1, 0, &ez);
      // ksw_extd(0, ql, qs, tl, ts, 5, mat, gapo, gape, gapo2, gape2, w, z1, 0, &ez);
      // ksw_extf2_sse(0, ql, qs, tl, ts, a, b, gape, w, z1, &ez);


      cigar.clear();
      for (int i = 0; i < ez.n_cigar; ++i) cigar += std::to_string(ez.cigar[i]>>4) + "MID"[ez.cigar[i]&0xf];

      free(ez.cigar); free(ts); free(qs);
    }

    // SAM format
    sam += data_set[fobj]->name + "\t" +
          "0" + "\t" + //FLAG
          reference[ref_num]->name + "\t" +
          std::to_string(ref_start + 1) + "\t" +
          "255" + "\t" + //MAPQ
          cigar + "\t" +
          "*" + "\t" + //RNEXT
          "0" + "\t" + //PNEXT
          "0" + "\t" +  //TLEN
          seq + "\t" +
          "*\n"; //QUAL

    // PAF format
    paf += data_set[fobj]->name + "\t" +
            std::to_string(sequence_length) + "\t" + 
            std::to_string(start) + "\t" + 
            std::to_string(end) + "\t" + 
            (rev ? "-" : "+") + "\t" +
            reference[ref_num]->name + "\t" +
            std::to_string(reference[ref_num]->sequence.size()) + "\t" +
            std::to_string(ref_start) + "\t" +
            std::to_string(ref_end) + "\t" +
            "*" + "\t" + //Number of residue matches
            std::to_string(end - start) + "\t" + //Alignment block length
            "255\n"; //Mapping quality 


  }
  // return sam;
  return paf;
  // return "";
}

int main (int argc, char **argv) {
  int optchr;

  int match = 2;
  int mismatch = -4;
  int gap = -2;

  int window_length = 10;
  int k = 15;
  int t = 2;
  float f = 0.0002f;
  bool bool_cigar = false;
  int k_value = 1000;

  while ((optchr = getopt_long(argc, argv, "hvm:g:M:k:w:t:cK:", long_options, NULL)) != -1) {
    switch (optchr) {
      case 'h': {
        help();
        exit(0);
      }
      case 'v': {
        version();
        exit(0);
      }
      case 'K': {
        k_value = atoi(optarg);
        if(k_value <= 0){
          fprintf(stderr, "[mapper] error: kvalue must be positive!\n"); 
          exit(1); 
        }
        break;
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
      case 'c': {
        bool_cigar = true;
        break;
      }
      case 't': {
        t = atoi(optarg);
        if(t < 1){
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
  std::vector<std::unique_ptr<FastAQ>> data_set;
  std::vector<std::unique_ptr<FastAQ>> reference;

  FastAQ::parse(data_set, file1, file1_format);
  FastAQ::parse(reference, file2, file2_format);

  fprintf(stderr, " Done!\n");

  // Print file stats
  // FastAQ::print_statistics(data_set, file1);
  // FastAQ::print_statistics(reference, file2);

  for(unsigned int ref_num = 0; ref_num < reference.size(); ref_num++){

    // Parallelized version of reference minimizers search
    std::shared_ptr<thread_pool::ThreadPool> thread_pool_ref = thread_pool::createThreadPool(t);
    std::vector<std::future<std::vector<triplet_t>>> thread_futures_ref;

    for (int tasks = 0; tasks < t - 1; ++tasks) {
      thread_futures_ref.emplace_back(thread_pool_ref->submit_task(brown::minimizers,
          reference[ref_num]->sequence.c_str() + tasks * reference[ref_num]->sequence.size() / t,
          reference[ref_num]->sequence.size() / t - 1 + window_length + k - 1,
          k, window_length));
    }
    thread_futures_ref.emplace_back(thread_pool_ref->submit_task(brown::minimizers,
          reference[ref_num]->sequence.c_str() + (t - 1) * reference[ref_num]->sequence.size() / t,
          reference[ref_num]->sequence.size() - (t - 1) * reference[ref_num]->sequence.size() / t,
          k, window_length));

    std::vector<triplet_t> t_minimizers;
    for (int i = 0; i < t; ++i) {
      thread_futures_ref[i].wait();
      unsigned int offset = i * reference[ref_num]->sequence.size() / t;
      for (auto& el : thread_futures_ref[i].get()) {
        std::get<1>(el) += offset;
        t_minimizers.push_back(el);
      }
    }

    // Regular version of reference minimizers search
    // std::vector<triplet_t> t_minimizers = brown::minimizers(
    //     reference[ref_num]->sequence.c_str(), reference[ref_num]->sequence.size(), k, window_length);

    prep_ref(t_minimizers, f);
    std::unordered_map<unsigned int, minimizer_index_t> ref_index = index_ref(t_minimizers);

    // prvih 250 - zbog testiranja
    // std::cout << map_paf(std::ref(t_minimizers), std::ref(ref_index), std::ref(data_set), std::ref(reference),
    //     k, window_length, match, mismatch, gap, 0, 250, k_value, bool_cigar, ref_num);


    // kod za visedretvenost
    std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool(t);

    std::vector<std::future<std::string>> thread_futures;

    for (int tasks = 0; tasks < t - 1; ++tasks) {
      thread_futures.emplace_back(thread_pool->submit_task(map_paf, std::ref(t_minimizers), std::ref(ref_index),
        std::ref(data_set), std::ref(reference),
        k, window_length, match, mismatch, gap, tasks * data_set.size() / t, (tasks + 1) * data_set.size() / t, k_value, bool_cigar, ref_num));  
    }
    thread_futures.emplace_back(thread_pool->submit_task(map_paf, std::ref(t_minimizers), std::ref(ref_index),
        std::ref(data_set), std::ref(reference),
        k, window_length, match, mismatch, gap, (t - 1) * data_set.size() / t, data_set.size(), k_value, bool_cigar, ref_num));
    
    for (auto& it : thread_futures) {
      it.wait();
      std::cout << it.get();
    }
  }
  return 0;
}