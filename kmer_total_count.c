// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kmer_utils.h"

int main(int argc, char **argv) {

  unsigned long long i = 0;

  if(argc !=  3) {
    printf("Please supply a filename and a kmer\n");
    exit(EXIT_FAILURE);
  }

	unsigned long long *counts = get_kmer_counts_from_file(argv[1], atoi(argv[2]));

	// print out our counts arrray
	// manually unrolled 4 loops to reduce fprintf calls
  for(i = 0; i < pow_four(atoi(argv[2])); i=i+4)
    fprintf(stdout, "%llu\n%llu\n%llu\n%llu\n", counts[i], counts[i+1], counts[i+2], counts[i+3]);

  return EXIT_SUCCESS;
}
