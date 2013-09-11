#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// This function takes a char array containing sequences and converts it
// into a kmer index (essentially a base 4 radix converstion to base 10, with
// some mumbo-jumbo of taking each of the characters we want (ACGT) and turning
// them into a radix-4 number we can use.
//
// kind of convoluted but it works.
//
// Arguemnts: 
//  char *str - a NULL terminated character array 
//  long kmer - how long of a index value you want to return
//  long error_pos - what index to return for a non ACGT character
//
long convert_kmer_to_index(char *str, long kmer, long error_pos) {

  char **ptr = NULL;
  char vals[kmer];
  long i = 0;
  long ret = 0;

  for(i = 0; i < kmer; i++) {
    int val = 0;
    switch(str[i]) {
      case 'a':
      case 'A':
        val = 48;
        break;
      case 'c':
      case 'C':
        val = 49;
        break;
      case 'g':
      case 'G':
        val = 50;
        break;
      case 't':
      case 'T':
        val = 51;
        break;
      default:
        return error_pos; 
    }


    vals[i] = val; 
  }

  ret = strtol(vals, ptr, 4);
  return ret;
}
