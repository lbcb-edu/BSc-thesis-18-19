#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <cctype>
#include <cstdint>
#include <iostream>
#include <string>
#include <set>
#include <tuple>

#include "srmapper.hpp"
#include "fastaq.hpp"
#include "brown_minimizers.hpp"
#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"

typedef std::tuple<uint64_t, uint32_t, bool> minimizer_t;

const std::set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
const std::set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};

static struct option long_options[] = {
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},
  {NULL, no_argument, NULL, 0}
};

void help(void) {
  printf("srmapper - tool for mapping short reads to reference genome.\n\n"

         "Usage: srmapper [OPTIONS] [reads reference]   start mapper\n"
         "  reads     - FASTA/FASTQ file containing a set of fragments\n"
         "  reference - FASTA file containing reference genome\n\n"

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
  );
}

void version(void) {
  printf("srmapper %d.%d\n",
    srmapper_VERSION_MAJOR,
    srmapper_VERSION_MINOR
  );
}

bool check_extension(const std::string& filename, const std::set<std::string>& extensions) {
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
  uint32_t k = 20;
  uint32_t w = 5;

  while ((optchr = getopt_long(argc, argv, "hvk:w:", long_options, NULL)) != -1) {
    switch (optchr) {
      case 'h': {
        help();
        exit(0);
      }
      case 'v': {
        version();
        exit(0);
      }
      case 'k': {
        k = atoi(optarg);
        break;
      }
      case 'w': {
        w = atoi(optarg);
        break;
      }
      default: {
        fprintf(stderr, "[mapper] error: Unknown option. Type %s --help for usage.\n", argv[0]);
        exit(1);
      }
    }
  }

  if (argc - optind != 3) {
    fprintf(stderr, "[srmapper] error: Expected 3 mapping arguments! Use --help for usage.\n");
    exit(1);
  }

  fprintf(stderr, "Loading reads and reference...\n");

  std::string reads_file1(argv[optind]);
  std::string reads_file2(argv[optind + 1]);
  std::string reference_file(argv[optind + 2]);

  if (!(check_extension(reads_file1, fasta_formats) || check_extension(reads_file1, fastq_formats))
      || !(check_extension(reads_file2, fasta_formats) || check_extension(reads_file2, fastq_formats))
      || !check_extension(reference_file, fasta_formats)) {
    fprintf(stderr, "[srmapper] error: Unsupported format(s)! Check --help for supported file formats.\n");
    exit(1);
  }


  std::vector<std::unique_ptr<fastaq::FastAQ>> reads1;
  std::vector<std::unique_ptr<fastaq::FastAQ>> reads2;
  std::vector<std::unique_ptr<fastaq::FastAQ>> reference;
  // fastaq::FastAQ::parse(reads1, reads_file1, check_extension(reads_file1, fasta_formats));
  // fastaq::FastAQ::parse(reads2, reads_file2, check_extension(reads_file2, fasta_formats));
  fastaq::FastAQ::parse(reference, reference_file, check_extension(reference_file, fasta_formats));

  fastaq::FastAQ::print_statistics(reference, reference_file);

  std::vector<minimizer_t> t_minimizers = brown::minimizers(reference[0]->sequence.c_str(), 
                                                            reference[0]->sequence.size(), 
                                                            k, w);

  std::cerr << t_minimizers.size() << std::endl;

  return 0;
}