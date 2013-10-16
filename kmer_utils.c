// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kmer_total_count.h"

const unsigned char alpha[256] = 
{5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 5, 1, 5, 5, 5, 2,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

inline unsigned long long pow_four(unsigned long long x) {
	return (unsigned long long)1 << (x * 2);
}

// convert a string of k-mer size base-4 values  into a
// base-10 index
inline unsigned long num_to_index(const char *str, const int kmer, const long error_pos) {

  int i = 0;
  unsigned long out = 0;
  unsigned long multiply = 1;

  for(i = kmer - 1; i >= 0; i--){

		if(str[i] >> 2) { 
			#ifndef SHARED
			position += i;
			#endif
			return error_pos;
		}

    out += str[i] * multiply;
    multiply = multiply << 2;
  }

  return out;
}

void convert_kmer_to_num(char *str, const unsigned long length) {

  unsigned long i = 0;

  for(i = 0; i < length; i++) {
    // this is _not_ portable, only works with ASCII values.
    switch(str[i] | 0x20 ) {
      case 'a':
        str[i] = 0;
        break;
      case 'c':
        str[i] = 1;
        break;
      case 'g':
        str[i] = 2;
        break;
      case 't':
        str[i] = 3;
        break;
      default:
        // Error Checking: use 4 so we can shift right twice
        // to check quickly is there is an non ACGT character 
        str[i] = 5;
    }

  }

}

char *strnstrip(const char *s, char *dest, int c, unsigned long long len) {

	unsigned long long i = 0;
	unsigned long long j = 0;

	for(i = 0; i < len; i++) {
		if(s[i] != c) {
			dest[j] = s[i];
			j++;
		}
	}

	dest[j] = '\0';

	return dest;
}

unsigned long long * get_kmer_counts_from_file(const char *fn, const unsigned int kmer) {
	
  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  long long i = 0;
	long long posistion = 0;

  FILE * const fh = fopen(fn, "r");
  if(fh == NULL) {
    fprintf(stderr, "Error opening %s - %s\n", fn, strerror(errno));
    exit(EXIT_FAILURE);
  }

	// width is 4^kmer  
	// there's a sneaky bitshift to avoid pow dependency
	const unsigned long width = pow_four(kmer); 

	// malloc our return array
  unsigned long long * counts = malloc((width+ 1) * sizeof(unsigned long long));
  if(counts == NULL)  {
		fprintf(stderr, strerror(errno));
    exit(EXIT_FAILURE);
	}

	char *str = malloc(4096);
	if(str == NULL) { 
		fprintf(stderr, strerror(errno));
		exit(EXIT_FAILURE);
	}

	unsigned long long str_size = 4096;

	while ((read = getdelim(&line, &len, '>', fh)) != -1) {

		// find our first \n, this should be the end of the header
		char *start = strchr(line, '\n');	
		if(start == NULL) 
			continue;

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

		// relace A, C, G and T with 0, 1, 2, 3 respectively
		// everything else is 5 
		for(i = 0; i < seq_length; i++) {
			str[i] = alpha[(int)str[i]];
		}

		// loop through our string to process each k-mer
		for(posistion = 0; posistion < (seq_length - kmer + 1); posistion++) {
			unsigned long mer = 0;
			unsigned long multiply = 1;

			// for each char in the k-mer check if it is an error char
			for(i = posistion + kmer - 1; i >= posistion; i--){
				if(str[i] >> 2) { 
					mer = width;
					posistion = i;
					goto next;
				}

				// multiply this char in the mer by the multiply
				// and bitshift the multiply for the next round
				mer += str[i] * multiply;
				multiply = multiply << 2;
			}
			// use this point to get mer of our loop
			next:
			// bump up the mer value in the counts array
			counts[mer]++;
		}
	} 

	return counts;
}
