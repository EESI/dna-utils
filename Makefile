VERSION=\"v0.0.1\"
CC = gcc
CFLAGS = -O3 -s -mtune=native -Wall -DVERSION=$(VERSION) -Wextra
CLIBS = -I. -L. -lkmer

all: libkmer.so libkmer.a kmer_total_count kmer_frequency_per_sequence

libkmer.so:
	$(CC) kmer_utils.c -o libkmer.so $(CFLAGS) -shared -fPIC -DSHARED=0
libkmer.a:
	ar rcs libkmer.a libkmer.so
	chmod +x libkmer.a
kmer_total_count:
	$(CC) kmer_total_count.c -o kmer_total_count $(CLIBS) $(CFLAGS) -DSHARED=0
kmer_frequency_per_sequence:
	$(CC) kmer_frequency_per_sequence.c -o kmer_frequency_per_sequence $(CLIBS) $(CFLAGS)

clean:
	rm -vf kmer_total_count kmer_frequency_per_sequence libkmer.so libkmer.a
