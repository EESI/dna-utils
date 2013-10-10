// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kmer_utils.h"
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

int main(int argc, char **argv) {

  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  long i = 0;
	long posistion = 0;

  if(argc != 3) {
    printf("Please supply a filename and a kmer\n");
    exit(EXIT_FAILURE);
  }

  FILE * const fh = fopen(argv[1], "r");
  if(fh == NULL) {
    fprintf(stderr, "Error opening %s - %s\n", argv[1], strerror(errno));
    exit(EXIT_FAILURE);
  }

  // second argument is the kmer
  const unsigned int kmer = atoi(argv[2]);

  // width is 4^kmer 
	const unsigned long width = (unsigned long)1 << (kmer * 2);

  // malloc our counts matrix
  unsigned long long * const counts = malloc((width+ 1) * sizeof(unsigned long long));

  if(counts == NULL) 
    exit(EXIT_FAILURE);
	
  while ((read = getline(&line, &len, fh)) != -1) {
    if(line[0] != '>' && read > kmer) {

  		for(i = 0; i < read; i++) {
				line[i] = alpha[line[i]];
			}

      for(posistion = 0; posistion < (read - kmer); posistion++) {
				unsigned long out = 0;
				unsigned long multiply = 1;
				

				for(i = posistion + kmer - 1; i >= posistion; i--){
					if(line[i] >> 2) { 
						out = width;
						posistion = i;
						goto next;
					}

					out += line[i] * multiply;
					multiply = multiply << 2;
				}
				next:
        counts[out]++;
      }
    } 
  }

  for(i = 0; i < (unsigned)width; i=i+4)
    printf("%llu\n%llu\n%llu\n%llu\n", counts[i], counts[i+1], counts[i+2], counts[i+3]);


  return EXIT_SUCCESS;
}
