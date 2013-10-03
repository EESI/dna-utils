VERSION=\"v0.0.1\"
CC = gcc
CFLAGS = -O3 -s -mtune=native -Wall -lm -DVERSION=$(VERSION) -Wextra

all: libkmer.so kmer_total_count kmer_frequency_per_sequence

libkmer.so:
	$(CC) kmer_utils.c -o libkmer.so $(CFLAGS) -shared -FPIC -DSHARED=0
kmer_total_count:
	$(CC) kmer_total_count.c kmer_utils.c -o kmer_total_count $(CFLAGS)
kmer_frequency_per_sequence:
	$(CC) kmer_frequency_per_sequence.c kmer_utils.c -o kmer_frequency_per_sequence $(CFLAGS)

clean:
	rm -vf kmer_total_count kmer_frequency_per_sequence libkmer.so
