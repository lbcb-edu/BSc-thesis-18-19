#include <getopt.h>
#include <stdio.h>

#include "OptSeqAlignmentLongGapsConfig.h"

struct option options[] = {
	{"version", no_argument, 0, 'v'},
	{"help", no_argument, 0, 'h'}
};

void version() {
	printf("v%d.%d.%d\n",
		OptSeqAlignmentLongGaps_VERSION_MAJOR,
		OptSeqAlignmentLongGaps_VERSION_MINOR,
		OptSeqAlignmentLongGaps_VERSION_PATCH);
}

void help() {
	//TO DO
	printf("help\n");
}

int main(int argc, char **argv) {
	
	char optchr;
	int option_index = 0;
	while((optchr = getopt_long(argc, argv, "hv", options, &option_index)) != -1) {
		printf("%c\n", optchr);
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

	return 0;
}
