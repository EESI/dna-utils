// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "kmer_utils.h"

void help() {
	printf("usage: kmer_counts_per_sequence input_file kmer [kmer-file] ...\n\n"
				 "count mers in each sequence of size k from a fasta file\n"
				 "\n"
				 "  --input    -i  input fasta file to count\n"
				 "  --kmer     -k  size of mers to count\n"
				 "  --mer-file -m  a file containing a list of mers you are interested\n"
				 "                 in opening. this will enable output your results in\n"
				 "                 a sparse format \n"
				 "  --sparse   -s  output values in a sparse format. output is in the\n"
				 "                 order sequence_number, mer_index, value\n"
				 "\n"
				 "Report all bugs to mutantturkey@gmail.com\n"
				 "\n"
				 "Copyright 2014 Calvin Morrison, Drexel University.\n"
				 "\n"
				 "If you are using any dna-utils tool for a publication\n"
				 "please cite your usage:\n\n"
				 "dna-utils. Drexel University, Philadelphia USA, 2014;\n"
				 "software available at www.github.com/EESI/dna-utils/\n");
}

int main(int argc, char **argv) {

	// getdelim variables
	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	// specific mer variables
	bool specific_mers = false;
	size_t num_desired_indicies = 0;
	size_t *desired_indicies = NULL;


	if(argc < 3) {
		printf("Please supply a filename and a kmer\n");
		exit(EXIT_FAILURE);
	}

	unsigned long kmer = atoi(argv[2]);
	if(kmer == 0) { 
		fprintf(stderr, "Error: invalid kmer.\n");
		exit(EXIT_FAILURE);
	}

	const unsigned long long width = (unsigned long long)1 << (kmer * 2);

	if(argc == 4) {
		specific_mers = true;
		desired_indicies = malloc((width) * sizeof(size_t));
		if(desired_indicies == NULL) 
			exit(EXIT_FAILURE);
		num_desired_indicies = load_specific_mers_from_file(argv[3], kmer, width, desired_indicies);

	}

	FILE *fh = fopen(argv[1], "r" );
	if(fh == NULL) {
		fprintf(stderr, "Error opening %s - %s\n", argv[1], strerror(errno));
		exit(EXIT_FAILURE);
	}
	unsigned long long *counts = malloc((width+ 1) * sizeof(unsigned long long));
	if(counts == NULL) 
		exit(EXIT_FAILURE);

	char *str = malloc(BUFSIZ);
	if(str == NULL) { 
		fprintf(stderr, strerror(errno));
		exit(EXIT_FAILURE);
	}

	unsigned long long str_size = BUFSIZ;
	unsigned long long sequence = 0;
	while ((read = getdelim(&line, &len, '>', fh)) != -1) {
		size_t i = 0;

		memset(counts, 0, width * sizeof(unsigned long long));

		// find our first \n, this should be the end of the header
		char *start = strchr(line, '\n');	
		if(start == NULL) 
			continue;

		// point to one past that.
		start = start + 1;

		size_t start_len = strlen(start);


		// if our current str buffer isn't big enough, realloc
		if(start_len + 1 > str_size + 1) { 
			str = realloc(str, start_len + 1);
			if(str == NULL) { 
				exit(EXIT_FAILURE);
				fprintf(stderr, strerror(errno));
			}
		}


		// strip out all other newlines to handle multiline sequences
		str = strnstrip(start, str, '\n',start_len);
		size_t seq_length = strlen(str);

		// reset our count matrix to zero

		for(i = 0; i < seq_length; i++) {
			str[i] = alpha[(int)str[i]];
		}

		for(i = 0; i < (seq_length - kmer + 1); i++) {
			size_t mer = num_to_index(&str[i],kmer, width, &i);
			counts[mer]++;
		}
		

		if(specific_mers) {
			for(i = 0; i < num_desired_indicies; i++)
				printf("%llu\t%zu\t%llu\n", sequence, desired_indicies[i], counts[desired_indicies[i]]);
			sequence++;
		} 
		else {
			for(i = 0; i < width - 1; i++)
				printf("%llu\t", counts[i]);
			printf("%llu\n", counts[width - 1]);
		}
	}

free(counts);
free(line);


return EXIT_SUCCESS;
}
