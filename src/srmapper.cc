#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <utility>
#include <chrono>

#include "srmapper.hpp"
#include "fastaq.hpp"
#include "mapping_params.hpp"
#include "index.hpp"
#include "map.hpp"
#include "brown_minimizers.hpp"
#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"

// Minimizer: value, position, origin
typedef std::tuple<uint64_t, uint32_t, bool> minimizer_t;
// Index: position, range
typedef std::pair<uint32_t, uint32_t> index_pos_t;
// Paired read
typedef std::pair<std::vector<std::unique_ptr<fastaq::FastAQ>>, std::vector<std::unique_ptr<fastaq::FastAQ>>> paired_reads_t;

// Accepted file formats
const std::unordered_set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
const std::unordered_set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};

static struct option long_options[] = {
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},
  {"paired", no_argument, NULL, 'p'},
  {"all", no_argument, NULL, 'a'},
  {"match", required_argument, NULL, 'm'},
  {"mismatch", required_argument, NULL, 'M'},
  {"gap-open", required_argument, NULL, 'o'},
  {"gap-extend", required_argument, NULL, 'e'},
  {"band", required_argument, NULL, 'b'},
  {"kmers", required_argument, NULL, 'k'},
  {"window_length", required_argument, NULL, 'w'},
  {"frequency", required_argument, NULL, 'f'},
  {"insert_size", required_argument, NULL, 'i'},
  {"threads", required_argument, NULL, 't'},
  {"threshold", required_argument, NULL, 'T'},
  {NULL, no_argument, NULL, 0}
};


void help(void) {
  printf("srmapper - tool for mapping short reads to reference genome.\n\n"

         "Usage: srmapper [OPTIONS] reference [reads]\n"
         "  reference - FASTA file containing reference genome\n"
         "  reads     - one or two FASTA/FASTQ file containing a set of fragments\n\n"

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
         "  -p  or  --paired         use pairing information (needs insert size information)\n"
         "  -a  or  --all            output all found mappings\n"
         "  -m  or  --match          <int>\n"
         "                             default: 1\n"
         "                             match value\n"
         "  -M  or  --mismatch       <int>\n"
         "                             default: -2\n"
         "                             mismatch value\n"
         "  -o  or  --gap-open       <int>\n"
         "                             default: 2\n"
         "                             gap open value\n"
         "  -e  or  --gap-extend     <int>\n"
         "                             default: 1\n"
         "                             gap extend value\n"
         "  -b  or  --band           <int>\n"
         "                             default: -1\n"
         "                             option: -2  => set band to half of QLEN\n"
         "                             ksw2 alignment band, band < 0 => disabled\n"
         "  -k  or  --kmers          <uint>\n"
         "                             default: 18\n"
         "                             constraints: largest supported is 32\n"
         "                             number of letters in substrings\n"
         "  -w  or  --window_length  <uint>\n"
         "                             default: 3\n"
         "                             length of window\n"
         "  -f  or  --frequency      <float>\n"
         "                             default: 0.001\n"
         "                             constraints: must be from [0, 1]\n"
         "                             number of most frequent minimizers that\n"
         "                             are not taken into account\n"
         "  -i  or  --insert_size    <uint>\n"
         "                             default: 215\n"
         "                             fragment insert size\n"
         "  -t  or  --threads        <uint>\n"
         "                             default: 3\n"
         "                             number of threads\n"
         "  -T  or  --threshold      <uint>\n"
         "                             default: 2\n"
         "                             number of hits needed in order to consider\n"
         "                             a region a candidate for mapping\n"
  );
}

void version(void) {
  printf("srmapper %d.%d\n",
    srmapper_VERSION_MAJOR,
    srmapper_VERSION_MINOR
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

int main(int argc, char **argv) {
  int optchr;
  mapping_params_t parameters;
  parameters.all = false;
  parameters.mch = 1;
  parameters.mis = -2;
  parameters.gapo = 2;
  parameters.gape = 1;
  parameters.band = -1;
  parameters.k = 18;
  parameters.w = 3;
  parameters.f = 0.001f;
  parameters.insert_size = 215;
  parameters.threshold = 2;
  bool paired = false;
  uint32_t t = 3;
  std::string cl_flags;

  while ((optchr = getopt_long(argc, argv, "hvpam:M:o:e:b:k:w:f:i:t:T:", long_options, NULL)) != -1) {
    cl_flags += "-", cl_flags += optchr, cl_flags += " ";
    if (optarg != nullptr) cl_flags += optarg, cl_flags += " ";
    switch (optchr) {
      case 'h': {
        help();
        exit(0);
      }
      case 'v': {
        version();
        exit(0);
      }
      case 'p': {
        paired = true;
        break;
      }
      case 'a': {
        parameters.all = true;
        break;
      }
      case 'm': {
        parameters.mch = atoi(optarg);
        break;
      }
      case 'M': {
        parameters.mis = atoi(optarg);
        break;
      }
      case 'o': {
        parameters.gapo = atoi(optarg);
        break;
      }
      case 'e': {
        parameters.gape = atoi(optarg);
        break;
      }
      case 'b': {
        parameters.band = atoi(optarg);
        break;
      }
      case 'k': {
        parameters.k = atoi(optarg);
        break;
      }
      case 'w': {
        parameters.w = atoi(optarg);
        break;
      }
      case 'f': {
        parameters.f = atof(optarg);
        if (parameters.f < 0.0f || parameters.f > 1.0f) {
          fprintf(stderr, "[srmapper] error: f must be from [0, 1].\n"); 
          exit(1); 
        }
        break;
      }
      case 'i': {
        parameters.insert_size = atoi(optarg);
        break;
      }
      case 't': {
        t = atoi(optarg);
        break;
      }
      case 'T': {
        parameters.threshold = atoi(optarg);
        break;
      }
      default: {
        fprintf(stderr, "[srmapper] error: Unknown option. Type %s --help for usage.\n", argv[0]);
        exit(1);
      }
    }
  }

  if (argc - optind > 3 || argc - optind < 2) {
    fprintf(stderr, "[srmapper] error: Expected read(s) and reference. Use --help for usage.\n");
    exit(1);
  }
  if (argc - optind == 2 && paired) {
    fprintf(stderr, "[srmapper] error: Expected paired reads in order to use pairing information (option -p).\n");
    exit(1);
  }

  auto i_start = std::chrono::steady_clock::now();

  fprintf(stderr, "[srmapper-load] loading reference... ");

  std::string reference_file(argv[optind]);
  if (!check_extension(reference_file, fasta_formats)) {
      fprintf(stderr, "[srmapper] error: Unsupported reference file format. Check --help for supported file formats.\n");
      exit(1);
  }
  std::vector<std::unique_ptr<fastaq::FastAQ>> reference;
  fastaq::FastAQ::parse(reference, reference_file, check_extension(reference_file, fasta_formats));

  fprintf(stderr, "\r[srmapper-load] loaded reference           \n"
                  "[srmapper-index] indexing reference... ");

  std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool(t);

  std::vector<std::future<std::vector<minimizer_t>>> thread_futures_ref;

  for (uint32_t tasks = 0; tasks < t - 1; ++tasks) {
    thread_futures_ref.emplace_back(thread_pool->submit_task(brown::minimizers,
        reference[0]->sequence.c_str() + tasks * reference[0]->sequence.size() / t,
        reference[0]->sequence.size() / t + parameters.w + parameters.k - 1,
        parameters.k, parameters.w));
  }
  thread_futures_ref.emplace_back(thread_pool->submit_task(brown::minimizers,
        reference[0]->sequence.c_str() + (t - 1) * reference[0]->sequence.size() / t,
        reference[0]->sequence.size() - (t - 1) * reference[0]->sequence.size() / t,
        parameters.k, parameters.w));

  std::vector<minimizer_t> t_minimizers;
  for (uint32_t i = 0; i < t; ++i) {
    thread_futures_ref[i].wait();
    uint32_t offset = i * reference[0]->sequence.size() / t;
    for (auto& el : thread_futures_ref[i].get()) {
      std::get<1>(el) += offset;
      t_minimizers.push_back(el);
    }
  }
  prep_ref(t_minimizers, parameters.f);
  std::unordered_map<uint64_t, index_pos_t> ref_index = index_ref(t_minimizers);

  fprintf(stderr, "\r[srmapper-index] indexed reference        \n");
  fastaq::FastAQ::print_statistics(reference, reference_file);

  auto i_end = std::chrono::steady_clock::now();

  auto i_interval = std::chrono::duration_cast<std::chrono::duration<double>>(i_end - i_start);

  std::cerr << "[srmapper-index] indexing time: " << i_interval.count() << " sec" << std::endl;
  
  auto m_start = std::chrono::steady_clock::now();

  std::cout << "@HD\tVN:1.6\n"
               "@SQ\tSN:" << reference[0]->name << "\tLN:" << reference[0]->sequence.size() << "\n"
               "@PG\tID:srmapper\tPN:srmapper\tCL:" << argv[0] << " " << cl_flags << argv[optind] << " ";

  if (argc - optind == 3) {
    fprintf(stderr, "[srmapper-load] loading paired-end reads... ");

    std::string reads_file1(argv[optind + 1]);
    std::string reads_file2(argv[optind + 2]);
    if (!(check_extension(reads_file1, fasta_formats) || check_extension(reads_file1, fastq_formats))
        || !(check_extension(reads_file2, fasta_formats) || check_extension(reads_file2, fastq_formats))) {
      fprintf(stderr, "[srmapper] error: Unsupported paired-end reads formats. Check --help for supported file formats.\n");
      exit(1);
    }
    paired_reads_t paired_reads;
    fastaq::FastAQ::parse(paired_reads.first, reads_file1, check_extension(reads_file1, fasta_formats));
    fastaq::FastAQ::parse(paired_reads.second, reads_file2, check_extension(reads_file2, fasta_formats));

    fprintf(stderr, "\r[srmapper-load] loaded paired-end reads        \n");
    fastaq::stats pr1_stats = fastaq::FastAQ::print_statistics(paired_reads.first, reads_file1);
    fastaq::stats pr2_stats = fastaq::FastAQ::print_statistics(paired_reads.second, reads_file2);

    if (paired) {
      fprintf(stderr, "[srmapper-map] using insert size information\n");
      if (pr1_stats.num != pr2_stats.num) {
        fprintf(stderr, "[srmapper] error: Paired-end read files must have equal number of reads (pairs).\n");
        exit(1);
      }
      if (pr1_stats.max - pr1_stats.min > 0 || pr2_stats.max - pr2_stats.min > 0) {
        fprintf(stderr, "[srmapper] warning: Reads are not of fixed size.\n");
      }
    }
    std::cout << argv[optind + 1] << " " << argv[optind + 2] << "\n";
    std::vector<std::future<std::string>> thread_futures;
      for (unsigned int tasks = 0; tasks < t - 1; ++tasks) {
        thread_futures.emplace_back(thread_pool->submit_task(paired ? map_paired : map_as_single, 
                std::ref(ref_index), std::ref(t_minimizers), std::ref(reference[0]), std::ref(paired_reads),
                std::ref(parameters), tasks * paired_reads.first.size() / t, (tasks + 1) * paired_reads.first.size() / t));  
      }
      thread_futures.emplace_back(thread_pool->submit_task(paired ? map_paired : map_as_single, 
                std::ref(ref_index), std::ref(t_minimizers), std::ref(reference[0]), std::ref(paired_reads),
                std::ref(parameters), (t - 1) * paired_reads.first.size() / t, paired_reads.first.size()));
      for (auto& it : thread_futures) {
        it.wait();
        std::cout << it.get();
      }
  } else {
    fprintf(stderr, "[srmapper-load] loading reads... ");
    std::string reads_file(argv[optind + 1]);
    if (!(check_extension(reads_file, fasta_formats) || check_extension(reads_file, fastq_formats))) {
      fprintf(stderr, "[srmapper] error: Unsupported format. Check --help for supported file formats.\n");
      exit(1);
    }
    std::vector<std::unique_ptr<fastaq::FastAQ>> reads;
    fastaq::FastAQ::parse(reads, reads_file, check_extension(reads_file, fasta_formats));

    fprintf(stderr, "\r[srmapper-load] loaded reads        \n");
    fastaq::FastAQ::print_statistics(reads, reads_file);

    std::cout << argv[optind + 1] << "\n";

    std::vector<std::future<std::string>> thread_futures;
    for (unsigned int tasks = 0; tasks < t - 1; ++tasks) {
      thread_futures.emplace_back(thread_pool->submit_task(map_single, std::ref(ref_index), std::ref(t_minimizers),
              std::ref(reference[0]), std::ref(reads),
              std::ref(parameters), tasks * reads.size() / t, (tasks + 1) * reads.size() / t));  
    }
    thread_futures.emplace_back(thread_pool->submit_task(map_single, std::ref(ref_index), std::ref(t_minimizers),
              std::ref(reference[0]), std::ref(reads),
              std::ref(parameters), (t - 1) * reads.size() / t, reads.size()));
    
    for (auto& it : thread_futures) {
      it.wait();
      std::cout << it.get();
    }
  }

  auto m_end = std::chrono::steady_clock::now();

  auto m_interval = std::chrono::duration_cast<std::chrono::duration<double>>(m_end - m_start);

  std::cerr << "[srmapper-map] mapping time: " << m_interval.count() << " sec" << std::endl;

  return 0;
}