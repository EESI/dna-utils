// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "kmer_utils.h"

int main(int argc, char **argv) {


	char *filename = NULL;

	unsigned int kmer = NULL;

	bool nonzero = 0;
	bool label = 0;
	bool kmer_set = 0;

	unsigned long long width = 0;

  unsigned long long i = 0;

  long long j = 0;

	static struct option long_options[] = {
		{"input", required_argument, 0, 'i'},
		{"kmer",  required_argument, 0, 'k'},
		{"nonzero", no_argument, 0, 'n'},
		{"label", no_argument, 0, 'l'},
		{"help", no_argument, 0, 'h'},
		{0, 0, 0, 0}
	};

	while (1) {

		int option_index = 0;
		int c = 0;

		c = getopt_long (argc, argv, "i:k:nlvh", long_options, &option_index);

		if (c == -1)
			break;

		switch (c) {
			case 'i':
				filename = optarg;
				break;
			case 'k':
				kmer = atoi(optarg);
				kmer_set = true;
				break;
			case 'n':
				nonzero = true; 
				break;
			case 'l':
				label = true;
				break;
			case 'h':
				printf("help-text\n");
				exit(EXIT_SUCCESS);
			default:
				break;
		}
	}
	if(filename == NULL) {
		fprintf(stderr, "Error: filename (-i) must be supplied\n");
		exit(EXIT_FAILURE);
	}
	if(kmer == NULL && !kmer_set) {
		fprintf(stderr, "Error: kmer (-k) must be supplied\n");
		exit(EXIT_FAILURE);
	}
	if(kmer == 0) { 
		fprintf(stderr, "Error: invalid kmer - '%d'.\n", kmer);
		exit(EXIT_FAILURE);
	}

	width = pow_four(kmer);

	unsigned long long *counts = get_kmer_counts_from_file(filename, kmer);

	// If nonzero is set, only print non zeros
	if(nonzero) {
		// if labels is set, print out our labels
		if(label) {
			for(i = 0; i < width; i++)
				if(counts[i] != 0) 
					fprintf(stdout, "%s\t%llu\n", index_to_kmer(i, kmer),i, counts[i]);
				
		}
		else {
			for(i = 0; i < width; i++)
				if(counts[i] != 0) 
					fprintf(stdout, "%llu\t%llu\n", i, counts[i]);
				
		}
	}
	// If we aren't printing nonzeros print everything
	else {
		if(label) {
			for(i = 0; i < width; i++)
				fprintf(stdout, "%s\t%llu\n", index_to_kmer(i, kmer), counts[i]);
		} 
		else {
		for(i = 0; i < width; i=i+4)
			fprintf(stdout, "%llu\n%llu\n%llu\n%llu\n", counts[i], counts[i+1], counts[i+2], counts[i+3]);
		}
	}
	return EXIT_SUCCESS;
}
