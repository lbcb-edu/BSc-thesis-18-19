#include <getopt.h>
#include <stdio.h>
#include <string>
#include <vector>

#include "bioparser/bioparser.hpp"

#include "OSALG_lib.h"

#include "OptSeqAlignmentLongGapsConfig.h"

struct option options[] = {
	{"version", no_argument, 0, 'v'},
	{"help", no_argument, 0, 'h'}
};

class FASTAQEntity {
	
	public:
		std::string name;
		std::string sequence;
		std::string quality;
		
		FASTAQEntity(
			const char *name, uint32_t name_length,
			const char *sequence, uint32_t sequence_length) {

			(this->name).assign(name, name_length);

			(this->sequence).assign(sequence, sequence_length);
		}

		FASTAQEntity(
			const char* name, uint32_t name_length,
			const char* sequence, uint32_t sequence_length,
			const char* quality, uint32_t quality_length) : FASTAQEntity(name, name_length, sequence, sequence_length) {
			
			(this->quality).assign(quality, quality_length);
		}
};

std::vector<std::string> FASTAExtensionVector{".fasta", ".fa", ".fasta.gz", ".fa.gz"};
std::vector<std::string> FASTQExtensionVector{".fastq", ".fq", ".fastq.gz", ".fq.gz"};

void version() {
	printf("v%d.%d.%d\n",
		OptSeqAlignmentLongGaps_VERSION_MAJOR,
		OptSeqAlignmentLongGaps_VERSION_MINOR,
		OptSeqAlignmentLongGaps_VERSION_PATCH);
}

void help() {
	printf(
		"\n"
		"==================\n"
		"The program accepts two files as floating arguments and outputs CIGAR strings gotten as result of sequence alignment allowing for long gaps.\n"
		"Program approximates concave function with 2 linear functions(u1=4, v1=1, u2=2, v2=13)\n"
		"Supported formats are .fasta, .fa, .fastq, .fq, .fasta.gz, .fa.gz, .fastq.gz, .fq.gz\n"
		"\n"
		"Usage: OptSeqAlignmentLongGaps <query_file_name> <reference_file_name>\n"
		"\n"
		"options: \n"
        	"	-h, --help  -  prints help menu (currently displaying)\n"
        	"	-v, --version  -  prints program version\n"
		"==================\n"
		"\n"
	);
}

bool endsWith(std::string const &fullString, std::string const &ending) {
	if(fullString.length() < ending.length()) return false;
	
	return (fullString.compare(fullString.length() - ending.length(), ending.length(), ending) == 0);
}

bool isExtensionMemberOfVector(std::string const &str, std::vector<std::string> const vec) {
	for(auto const& s : vec) {
		if(endsWith(str, s)) return true;
	}

	return false;
}

std::vector<std::unique_ptr<FASTAQEntity>> readFASTQFile(std::string const &filePath) {
	std::vector<std::unique_ptr<FASTAQEntity>> fastq_objects;
	auto fastq_parser = bioparser::createParser<bioparser::FastqParser, FASTAQEntity>(filePath);

	uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
	while (true) {
		auto status = fastq_parser->parse(fastq_objects, size_in_bytes);
		if (status == false) {
			break;
		}
	}

	return fastq_objects;
}

std::vector<std::unique_ptr<FASTAQEntity>> readFASTAFile(std::string const &filePath) {
	std::vector<std::unique_ptr<FASTAQEntity>> fasta_objects;
	auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAQEntity>(filePath);
	
	fasta_parser->parse(fasta_objects, -1);
	
	return fasta_objects;
}

int main(int argc, char **argv) {
	
	char optchr;
	int option_index = 0;
	while((optchr = getopt_long(argc, argv, "hv", options, &option_index)) != -1) {
		switch(optchr) {
		
			case 0:
				break;
			case 'h':
				help();
				break;
			case 'v':
				version();
				break;
			default:
				fprintf(stderr, "Entered option is not valid.\n");
				fprintf(stderr, "Use \"-h\" or \"--help\" for more information.\n");
				return 1;
		}
	}

	if(argc - optind < 2) {
		fprintf(stderr, "Program requires more than one argument!\n");
		fprintf(stderr, "Use \"-h\" or \"--help\" for more information.\n");
		return 1;
	}

	std::string firstFilePath = argv[optind];
	std::string secondFilePath = argv[optind + 1];

	bool isFirstFASTA = isExtensionMemberOfVector(firstFilePath, FASTAExtensionVector);

	if(!((isFirstFASTA || isExtensionMemberOfVector(firstFilePath, FASTQExtensionVector))
		&& isExtensionMemberOfVector(secondFilePath, FASTAExtensionVector))) {
		fprintf(stderr, "One or more given arguments are not valid file format!\n");
		fprintf(stderr, "Use \"-h\" or \"--help\" for more information.\n");
		return 1;
	}


	std::vector<std::unique_ptr<FASTAQEntity>> query_entries;
	if(isFirstFASTA) {
		query_entries = readFASTAFile(firstFilePath);
	} else {
		query_entries = readFASTQFile(firstFilePath);
	}

	std::vector<std::unique_ptr<FASTAQEntity>> ref_entries = readFASTAFile(secondFilePath);

	for(auto const &query : query_entries) {
		for(auto const &ref : ref_entries) {
			std::vector<std::string> cigars;
			OSALG::long_gaps_alignment(query->sequence, ref->sequence, cigars);

			printf("----------------------\n");
			printf("reference: %s\n", (ref->name).c_str());
			printf("query: %s\n", (query->name).c_str());
			printf("cigars: \n");
			for(auto const &cig : cigars) {
				printf("--$$$$$--\n");
				printf("%s\n", cig.c_str());
				printf("--$$$$$--\n");
			}
			printf("----------------------\n");
		}
	}

	//std::vector<std::string> cigars;
	//OSALG::long_gaps_alignment("AGTGATCCCTGAAGTTGC", "TGAGTTTCAGCTCAAAGGGTCTCGATCTCC", cigars);
	//OSALG::long_gaps_alignment("AGTGATCC", "TGAGTTTCAG", cigars);
	//OSALG::long_gaps_alignment("GATCC", "TTCAG", cigars);
	//OSALG::long_gaps_alignment("AGTGAC", "TGAGCTT", cigars);

	//for(int i = 0; i < cigars.size(); ++i) {
	//	printf("%s\n", cigars[i].c_str());
	//}

	return 0;
}
