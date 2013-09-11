// Copyright 2013 Calvin Morrison
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "kmer_utils.h"

int main(int argc, char **argv) {

  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  unsigned int i = 0;

  if(argc != 3) {
    printf("Please supply a filename, and only a filename\n");
    exit(EXIT_FAILURE);
  }

  FILE *fh = fopen(argv[1], "r" );
  if(fh == NULL) {
    fprintf(stderr, "Couldn't open: %s\n", argv[1]);
    exit(EXIT_FAILURE);
  }
  
  int kmer = atoi(argv[2]);
  int width = (int)pow(4, kmer);

  unsigned long long *counts = malloc((width+ 1) * sizeof(unsigned long long));
  if(counts == NULL) 
    exit(EXIT_FAILURE);

  while ((read = getline(&line, &len, fh)) != -1) {
    if(line[0] != '>')   {

      for(i = 0; i < strlen(line) - kmer; i++) {
        counts[convert_kmer_to_index(&line[i],kmer, width)]++;
      }
    }
  }

  for(i = 0; i < width; i++)
    printf("%llu\n", counts[i]);


  return EXIT_SUCCESS;
}
