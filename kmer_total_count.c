// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kmer_utils.h"

int main(int argc, char **argv) {

  unsigned long long i = 0;
  long long j = 0;

  if(argc <  3) {
    fprintf(stderr, "Please supply a filename and a kmer, then as many kmers as you want (in ACGT format)\n");
    exit(EXIT_FAILURE);
  }

	unsigned long kmer = atoi(argv[2]);
	unsigned long long width = pow_four(kmer);

	if(kmer == 0) { 
		fprintf(stderr, "Error: invalid kmer.\n");
		exit(EXIT_FAILURE);
	}

	unsigned long long *counts = get_kmer_counts_from_file(argv[1], atoi(argv[2]));

	// print out our counts arrray
	// manually unrolled 4 loops to reduce fprintf calls
	
	if(argc == 3) {
		for(i = 0; i < width; i=i+4)
			fprintf(stdout, "%llu\n%llu\n%llu\n%llu\n", counts[i], counts[i+1], counts[i+2], counts[i+3]);
	} 
	else {
		for(i = 3; i < argc; i++) { 
			unsigned long k = width;	
			size_t len = 0;
			fprintf(stdout, "%s\t", argv[i]);
			
			len = strlen(argv[i]);
			if(len != kmer) { 
				fprintf(stdout, "ERR\n");
				continue;
			}
			for(j = 0; j < len; j++)
				argv[i][j] = alpha[argv[i][j]];

			k = num_to_index(argv[i], kmer, width);

			if(k == width)  {
				fprintf(stdout, "ERR\n");
			}
			else {
				fprintf(stdout, "%llu\n", counts[k]);
			}
		}
	}
	return EXIT_SUCCESS;
}
