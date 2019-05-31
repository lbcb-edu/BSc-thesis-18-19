#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <cctype>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
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
#include <chrono>

#include "long_reads_mapper.hpp"
#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "ksw2.h"

#include "brown_minimizers.hpp"
#include "fastaq.hpp"
#include "index.hpp"
#include "util.hpp"

typedef std::tuple<unsigned int, unsigned int, bool> minimizer; // value, position, origin
typedef std::pair<unsigned int, unsigned int> minimizer_index_t; // position, amount
typedef std::pair<unsigned int, unsigned int> minimizer_hit_t; // query position, reference position
typedef std::pair<unsigned int, int> region_hits; //region number, number of hits

// bioparser ptr
typedef std::unique_ptr<bioparser::Parser<FastAQ>> parser_ptr_t;

const std::unordered_set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
const std::unordered_set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};

std::map<char, char> complement_map = {{'C', 'G'}, {'A', 'T'}, {'T', 'A'}, {'U', 'A'}, {'G', 'C'}};

// Sample size in bytes used for insert size inferrence
constexpr uint32_t sample_bytes = 512 * 1024 * 1024;

int k_value = 1000;
float f = 0.001f; // za odbacivanje minimizera
float percentage = 0.5; // koliki postotak sekvence gledati ako nema dovoljan broj hitova

int window_length = 5;
int k = 15;
int bandwidth = 500;
int threshold = 5;

int match = 2;
int mismatch = -4;
int gap_open = 4;
int gap_extension = 2;
int z_drop = 400;

bool bool_cigar = false;
bool sam_format = false;
int t = 3;

static struct option long_options[] = {
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},
  {"kvalue", required_argument, NULL, 'K'},
  {"kmers", required_argument, NULL, 'k'},
  {"window_length", required_argument, NULL, 'w'},
  {"bandwidth", required_argument, NULL, 'r'},
  {"min_hits", required_argument, NULL, 'n'},
  {"match", required_argument, NULL, 'A'},
  {"mismatch", required_argument, NULL, 'B'},
  {"gap_open", required_argument, NULL, 'O'},
  {"gap_extension", required_argument, NULL, 'E'},
  {"z_drop", required_argument, NULL, 'z'},
  {"sam", no_argument, NULL, 'a'},
  {"cigar", no_argument, NULL, 'c'},
  {"threads", required_argument, NULL, 't'},
  {NULL, no_argument, NULL, 0}
};

void help(void) {
  printf("long_reads_mapper - for mapping long fragments on reference genome.\n\n"
         "Usage: long_reads_mapper [OPTIONS] [file1 file2]\n"
         "file1 - FASTA file containing reference genome\n\n"
         "file2 - FASTA/FASTQ file containing a set of fragments\n"
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
         "  -w  or  --window_length  <int>\n"
         "                             default: 5\n"
         "                             length of window\n"
         "  -k  or  --kmers          <int>\n"
         "                             default: 15\n"
         "                             constraints: largest supported is 16\n"
         "                             number of letters in substrings\n"
         "  -r  or  --bandwidth      <int>\n"
         "                             default: 500\n"
         "                             bandwidth used in chaining and DP-based alignment\n"
         "  -n  or  --min_hits       <int>\n"
         "                             default: 5\n"
         "                             minimum hits when finding regions\n"
         "  -A  or  --match          <int>\n"
         "                             default: 2\n"
         "                             match number\n"
         "  -B  or  --mismatch       <int>\n"
         "                             default: -4\n"
         "                             mismatch number\n"
         "  -O  or  --gap_open       <int>\n"
         "                             default: 4\n"
         "                             gap open penalty\n"
         "  -E  or  --gap_extension  <int>\n"
         "                             default: 2\n"
         "                             gap extension penalty\n"
         "  -z  or  --z_drop         <int>\n"
         "                             default: 400\n"
         "                             Z-drop score\n"
         "  -c  or  --cigar          output CIGAR in PAF\n"
         "  -a  or  --sam            output in the SAM format (PAF by default)\n"
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

// Check file extension
// Args: filename   - name of file to be checked for extension
//       extensions - set of accepted extensions
// Return: extension accepted or not accepted
bool check_extension(const std::string& filename, const std::unordered_set<std::string>& extensions) {
  for (const auto& it : extensions) {
    if (filename.size() > it.size()) {
      if (filename.compare(filename.size()-it.size(), std::string::npos, it) == 0) {
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

std::string map_paf(
    const std::vector<minimizer>& t_minimizers,
    const std::unordered_map<unsigned int, minimizer_index_t>& ref_index,
    const std::vector<std::unique_ptr<FastAQ>>& data_set,
    const std::vector<std::unique_ptr<FastAQ>>& reference,
    unsigned int t_begin, unsigned int t_end, unsigned int ref_num) {

  std::string result;

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

    std::string seq;
    if(rev){
      seq = reverse_complement(data_set[fobj]->sequence, start, end - start);
    }else{
      seq = data_set[fobj]->sequence.substr(start, end - start);
    }

    //ksw2
    std::string cigar;
    if(bool_cigar || sam_format){
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

      ksw_extz2_sse(0, ql, qs, tl, ts, 5, mat, gap_open, gap_extension, bandwidth, z_drop, 0, 0, &ez);

      // int gapo = 4;
      // int gapo2 = 24;
      // int gape = 1;
      // int gape2 = 2;
      // int w = 500;
      // int z1 = 400;
      // int z2 = 200;
      // int s = 80;
      // int noncan = 0;
      // ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, gapo2, gape2, w, z1, 0, 0, &ez);
      // ksw_exts2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, gapo2, noncan, z1,  0, &ez); //ispis??

      for (int i = 0; i < ez.n_cigar; ++i) cigar += std::to_string(ez.cigar[i]>>4) + "MID"[ez.cigar[i]&0xf];
      free(ez.cigar); free(ts); free(qs);
    }

    //PAF format
    if(!sam_format){ 
      result += data_set[fobj]->name + "\t" +
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
            "255" + "\t" + //Mapping quality
            (bool_cigar ? ("cg:Z:" + cigar + "\n") : "\n");
    //SAM format
    }else{ 
      result += data_set[fobj]->name + "\t" +
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
    }
  }
  return result;
}

int main (int argc, char **argv) {
  int optchr;

  while ((optchr = getopt_long(argc, argv, "hvK:w:k:r:n:A:B:O:E:z:t:ca", long_options, NULL)) != -1) {
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
          fprintf(stderr, "[ERROR]: kvalue must be positive!\n"); 
          exit(1); 
        }
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
      case 'r': {
        bandwidth = atoi(optarg);
        break;
      }
      case 'n': {
        threshold = atoi(optarg);
        break;
      }
      case 'A': {
        match = atoi(optarg);
        break;
      }
      case 'B': {
        mismatch = atoi(optarg);
        break;
      }
      case 'O': {
        gap_open = atoi(optarg);
        break;
      }
      case 'E': {
        gap_extension = atoi(optarg);
        break;
      }
      case 'z': {
        z_drop = atoi(optarg);
        break;
      }
      case 'c': {
        bool_cigar = true;
        break;
      }
      case 'a': {
        sam_format = true;
        break;
      }
      case 't': {
        t = atoi(optarg);
        if(t < 1){
          fprintf(stderr, "[ERROR]: t must be positive!\n"); 
          exit(1); 
        }
        break;
      }
      default: {
        fprintf(stderr, "[ERROR]: Unknown option. Type %s --help for usage.\n", argv[0]);
        exit(1);
      }
    }
  }

  if (argc - optind != 2) {
    fprintf(stderr, "[ERROR]: Expected 2 arguments! Use --help for usage.\n");
    exit(1);
  }

  auto start_time = std::chrono::steady_clock::now();

  std::string reference_file(argv[optind]);
  std::string data_file(argv[optind + 1]);

  if (!check_extension(reference_file, fasta_formats)) {
      fprintf(stderr, "[ERROR]: Unsupported reference file format.\n");
      exit(1);
  }
    if (!check_extension(data_file, fasta_formats) && !check_extension(data_file, fastq_formats)) {
      fprintf(stderr, "[ERROR]: Unsupported data file format.\n");
      exit(1);
  }

  std::vector<std::unique_ptr<FastAQ>> reference;
  parser_ptr_t ref_parser = bioparser::createParser<bioparser::FastaParser, FastAQ>(reference_file);
  ref_parser->parse_objects(reference, -1);

  std::vector<std::unique_ptr<FastAQ>> data_set;
  parser_ptr_t data_parser;
  if(check_extension(data_file, fasta_formats)) {
    data_parser = bioparser::createParser<bioparser::FastaParser, FastAQ>(data_file);
  }else {
    data_parser = bioparser::createParser<bioparser::FastqParser, FastAQ>(data_file);
  }
  data_parser->parse_objects(data_set, -1);

  auto loading_time = std::chrono::steady_clock::now();

  FastAQ::print_statistics(reference, reference_file);
  fprintf(stderr,"\n");
  FastAQ::print_statistics(data_set, data_file);
  fprintf(stderr,"\n");

  fprintf(stderr, "Mapping process started with parameters:\n"
                  "  k value        = %d\n"
                  "  k              = %d\n"
                  "  window length  = %d\n"
                  "  bandwidth      = %d\n"
                  "  min hits       = %d\n"
                  "  threads        = %d\n",
                  k_value, k, window_length, bandwidth, threshold, t);
  if (bool_cigar || sam_format) {
    fprintf(stderr, " Alignment\n"
                    "  match          = %d\n"
                    "  mismatch       = %d\n"
                    "  gap open       = %d\n"
                    "  gap extension  = %d\n"
                    "  z drop         = %d\n",
                    match, mismatch, gap_open, gap_extension, z_drop);
  }
  fprintf(stderr,"\n");

  for(unsigned int ref_num = 0; ref_num < reference.size(); ref_num++){

    std::shared_ptr<thread_pool::ThreadPool> thread_pool_ref = thread_pool::createThreadPool(t);
    std::vector<std::future<std::vector<minimizer>>> thread_futures_ref;

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

    std::vector<minimizer> t_minimizers;
    for (int i = 0; i < t; ++i) {
      thread_futures_ref[i].wait();
      unsigned int offset = i * reference[ref_num]->sequence.size() / t;
      for (auto& el : thread_futures_ref[i].get()) {
        std::get<1>(el) += offset;
        t_minimizers.push_back(el);
      }
    }

    prep_ref(t_minimizers, f);
    std::unordered_map<unsigned int, minimizer_index_t> ref_index = index_ref(t_minimizers);

    std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool(t);
    std::vector<std::future<std::string>> thread_futures;

    for (int tasks = 0; tasks < t - 1; ++tasks) {
      thread_futures.emplace_back(thread_pool->submit_task(map_paf, std::ref(t_minimizers), std::ref(ref_index),
        std::ref(data_set), std::ref(reference), tasks * data_set.size() / t, (tasks + 1) * data_set.size() / t, ref_num));  
    }
    thread_futures.emplace_back(thread_pool->submit_task(map_paf, std::ref(t_minimizers), std::ref(ref_index),
        std::ref(data_set), std::ref(reference), (t - 1) * data_set.size() / t, data_set.size(), ref_num));
    
    for (auto& it : thread_futures) {
      it.wait();
      std::cout << it.get();
    }
  }

  auto end_time = std::chrono::steady_clock::now();
  auto interval1 = std::chrono::duration_cast<std::chrono::duration<double>>(loading_time - start_time);
  auto interval2 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - loading_time);
  fprintf(stderr, "Time spend on loading files: %.2f sec\n", interval1.count());
  fprintf(stderr, "Time spend on mapping: %.2f sec\n", interval2.count());
  
  return 0;
}