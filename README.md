# dna-utils 


This repository contains general utilities for processing sequences in fasta files.


### Tools included ###


kmer_total_count - this program will count each kmer in a fasta file, and print to standard out.

#### Usage

    kmer_total_count filename kmer_size
    
    # Example
    calvin@barnabas:~/dna-utils$ ./kmer_total_count SuperManSequences.fasta 8 
 
the order of the array corresponds to AAAA, AAAC, AAAG, AAAT and so on

usage: kmer_total_count -i input_file -k kmer [-n] [-l] ...

#### Licensing and Citation
Report all bugs to mutantturkey@gmail.com

Copyright 2014 Calvin Morrison, Drexel University.

If you are using any dna-utils tool for a publication
please cite your usage:

dna-utils. Drexel University, Philadelphia USA, 2014;
software available at www.github.com/EESI/dna-utils/

