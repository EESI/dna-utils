VERSION=\"v0.0.1\"
CC = gcc
CFLAGS = -O3 -s -mtune=native -Wall -DVERSION=$(VERSION) -Wextra
CLIBS = libkmer.a

all: libkmer.o libkmer.so libkmer.a kmer_total_count kmer_counts_per_sequence

libkmer.o: kmer_utils.c
	$(CC) -c kmer_utils.c -o libkmer.o $(CFLAGS) -fPIC -DSHARED=0
libkmer.so: libkmer.o
	$(CC) kmer_utils.c -o libkmer.so $(CFLAGS) -shared -fPIC -DSHARED=0
libkmer.a: libkmer.o
	ar rcs libkmer.a libkmer.o
	chmod +x libkmer.a
kmer_total_count: libkmer.a kmer_total_count.c kmer_utils.h
	$(CC) kmer_total_count.c -o kmer_total_count $(CLIBS) $(CFLAGS) -DSHARED=0
kmer_counts_per_sequence: libkmer.a kmer_counts_per_sequence.c kmer_utils.h
	$(CC) kmer_counts_per_sequence.c -o kmer_counts_per_sequence $(CLIBS) $(CFLAGS)

clean:
	rm -vf kmer_total_count kmer_counts_per_sequence libkmer.so libkmer.a libkmer.o
