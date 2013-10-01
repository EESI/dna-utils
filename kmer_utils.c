#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// convert a string of k-mer size base-4 values  into a
// base-10 index
long num_to_index(const char *str, const int kmer, const long error_pos) {

  int i = 0;
  unsigned long out = 0;
  unsigned long multiply = 1;

  for(i = kmer - 1; i >= 0; i--){

    if(str[i] >> 2)
      return error_pos;

    out += str[i] * multiply;
    multiply = multiply << 2;
  }

  return out;
}

// replaces values in a char array of ACGT's and others with
// values that correspond to their base 4 value to be used in
// num_to_index.
void convert_kmer_to_num(char *str, const long length) {

  long i = 0;

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
        str[i] = 4;
    }

  }

}
