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
inline long convert_kmer_to_index(const char *str, long kmer, long error_pos) {

  long i = 0;
  long out = 0;
  long multiply = 1;

  for(i = kmer - 1; i >= 0; i--){
    switch(str[i]) {
      case 'a':
      case 'A':
        break;
      case 'c':
      case 'C':
        out += 1 * multiply;
        break;
      case 'g':
      case 'G':
        out += 2 * multiply;
        break;
      case 't':
      case 'T':
        out += 3 * multiply;
        break;
      default:
        return error_pos; 
    }

    multiply = multiply << 2;
  }

  return out;
}
