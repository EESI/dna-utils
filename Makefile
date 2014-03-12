VERSION=\"0.0.4\"
CC = gcc
CFLAGS = -O3 -s -mtune=native -Wall -DVERSION=$(VERSION) -Wextra
DESTDIR = /usr/local/


all: libkmer.o libkmer.so kmer_total_count kmer_counts_per_sequence kmer_continuous_count

libkmer.o: kmer_utils.c
	$(CC) -c kmer_utils.c -o libkmer.o $(CFLAGS) -fPIC
libkmer.so: libkmer.o
	$(CC) kmer_utils.c -o libkmer.so $(CFLAGS) -shared -fPIC
kmer_total_count: libkmer.o kmer_total_count.c kmer_utils.h
	$(CC) libkmer.o kmer_total_count.c -o kmer_total_count $(CLIBS) $(CFLAGS)
kmer_counts_per_sequence: libkmer.o kmer_counts_per_sequence.c kmer_utils.h
	$(CC) libkmer.o kmer_counts_per_sequence.c -o kmer_counts_per_sequence $(CLIBS) $(CFLAGS)
kmer_continuous_count: libkmer.o kmer_continuous_count.c kmer_utils.h
	$(CC) libkmer.o kmer_continuous_count.c -o kmer_continuous_count $(CLIBS) $(CFLAGS)

clean:
	rm -vf kmer_total_count kmer_counts_per_sequence kmer_continuous_count libkmer.so libkmer.o

debug: CFLAGS = -ggdb -Wall -Wextra -DVERSION=$(VERSION)\"-debug\"
debug: all

install: all
	@cp -vf kmer_counts_per_sequence kmer_total_count kmer_continuous_count $(DESTDIR)/bin/
	@cp -vf  libkmer.so  $(DESTDIR)/lib/

