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
inline long convert_kmer_to_index(char *str, unsigned long kmer, long error_pos) {

  char **ptr = NULL;
  char vals[kmer];
  unsigned long i = 0;

  for(i = 0; i < kmer; i++) {
    switch(str[i]) {
      case 'a':
      case 'A':
        vals[i] = '0';
        break;
      case 'c':
      case 'C':
        vals[i] = '1';
        break;
      case 'g':
      case 'G':
        vals[i] = '2';
        break;
      case 't':
      case 'T':
        vals[i] = '3';
        break;
      default:
        return error_pos; 
    }
  }

  return strtol(vals, ptr, 4);
}
