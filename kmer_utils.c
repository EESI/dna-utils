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
inline long convert_kmer_to_index(const char *str, int kmer, long error_pos) {

  int i = 0;
  unsigned long out = 0;
  unsigned long multiply = 1;

  for(i = kmer - 1; i >= 0; i--){
    switch(str[i] | 0x20 ) {
      case 'a':
        break;
      case 'c':
        out += 1 * multiply;
        break;
      case 'g':
        out += 2 * multiply;
        break;
      case 't':
        out += 3 * multiply;
        break;
      default:
        return error_pos; 
    }

    multiply = multiply << 2;
  }

  return out;
}
