// Copyright 2013 Calvin Morrison
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "kmer_utils.h"

long position = 0;

int main(int argc, char **argv) {

  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  long i = 0;

  unsigned long long *counts;

  if(argc != 3) {
    printf("Please supply a filename and a kmer\n");
    exit(EXIT_FAILURE);
  }

  FILE *fh = fopen(argv[1], "r" );
  if(fh == NULL) {
    fprintf(stderr, "Couldn't open: %s\n", argv[1]);
    exit(EXIT_FAILURE);
  }

  // second argument is the kmer
  const unsigned int kmer = atoi(argv[2]);

  // width is 4^kmer 
  const unsigned long width = (int)pow(4, kmer);

  // malloc our counts matrix
  counts = malloc((width+ 1) * sizeof(unsigned long long));
  if(counts == NULL) 
    exit(EXIT_FAILURE);

  while ((read = getline(&line, &len, fh)) != -1) {
    if(line[0] != '>' && read > kmer) {
      convert_kmer_to_num(line, read);

      for(position = 0; position < (read - kmer); position++) {
        counts[num_to_index(&line[position],kmer, width)]++;
      }
    } 
  }

  for(i = 0; i < (unsigned)width; i++)
    printf("%llu\n", counts[i]);


  return EXIT_SUCCESS;
}
