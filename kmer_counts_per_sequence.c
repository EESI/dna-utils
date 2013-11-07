// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kmer_utils.h"

unsigned long position = 0;
int main(int argc, char **argv) {

	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	if(argc != 3) {
		printf("Please supply a filename and a kmer\n");
		exit(EXIT_FAILURE);
	}

	FILE *fh = fopen(argv[1], "r" );
	if(fh == NULL) {
		fprintf(stderr, "Error opening %s - %s\n", argv[1], strerror(errno));
		exit(EXIT_FAILURE);
	}

	unsigned long kmer = atoi(argv[2]);
	if(kmer == 0) { 
		fprintf(stderr, "Error: invalid kmer.\n");
		exit(EXIT_FAILURE);
	}

	const unsigned long width = (unsigned long)1 << (kmer * 2);

	unsigned long long *counts = malloc((width+ 1) * sizeof(unsigned long long));
	if(counts == NULL) 
		exit(EXIT_FAILURE);

	while ((read = getline(&line, &len, fh)) != -1) {
		if(line[0] != '>' && (read > kmer)) {

			unsigned int i = 0;
			unsigned long total = 0;

			// reset our count matrix to zero
			memset(counts, 0, width * sizeof(unsigned long long));

			for(i = 0; i < read - kmer; i++) {
				line[i] = alpha[(int)line[i]];
			}

			for(i = 0; i < read - kmer; i++) {
				counts[num_to_index(&line[i],kmer, width)]++;
			}

			for(i = 0; i < width; i++)
				total += counts[i];

			for(i = 0; i < width - 1; i++)
				printf("%llu\t", counts[i]);
			printf("%llu\n", counts[width - 1]);

		}
	}

	free(counts);
	free(line);


	return EXIT_SUCCESS;
}
