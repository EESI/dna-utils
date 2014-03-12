// dna-util's function library.

unsigned long num_to_index(const char *str, const int kmer, const long error_pos, long long *current_position);
char *index_to_kmer(unsigned long long index, long kmer);

// Utility functions


// strip char 'c' out of char array *s of length len
size_t strnstrip(char *s, int c, size_t len);

// reverse char arry *s of length len
void reverse_string(char *s, size_t len);

// quicky calculate 4^x 
unsigned long long pow_four(unsigned long long x);

// check if pointer is null. a helper for dealing with NULL
// return values as errors. Calls strerror and quits if 
// ptr is null, optionally takes *error char array as 
// a error to output
void check_null_ptr(void *ptr, const char *error);

void count_sequence(const char *seq, const size_t seq_length, const unsigned int kmer, unsigned long long *counts);

// Variables
//
const unsigned char alpha[256]; 
const unsigned char reverse_alpha[4];
const unsigned char compliment[5];


// file loading functions

// open file from filename in char array *fn, and try and parse in one mer per
// line, of size kmer, and store the indicies of those mers in the *arr
// pointer;
unsigned long long load_specific_mers_from_file(const char *fn, unsigned int kmer, size_t width, size_t *arr);

unsigned long long * get_kmer_counts_from_filename(const char *fn, const unsigned int kmer, const bool count_compliment);
unsigned long long * get_kmer_counts_from_file(FILE *fh, const int kmer, const bool count_compliment);

unsigned long long * get_continuous_kmer_counts_from_filename(const char *fn, const unsigned int kmer, const bool count_compliment);
unsigned long long * get_continuous_kmer_counts_from_file(FILE *fh, const unsigned int kmer, const bool count_compliment);

