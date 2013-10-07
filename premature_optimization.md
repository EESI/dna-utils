# Premature Optimization #

Often you hear the ringing refrain that 'premature opimization is the root of all evil' in one form or another. In this post I will take you down the path of 'Optimization is the result of bordom.'

Here's the best I could do with 5 minutes of google search:
![dilbert](http://joshreads.com/images/07/09/i070903dilbert.png "Man Hours")

## How it all started ##
During a vacation I recently took, I spent quite a bit of time in the airports and flying. During that time I was unable to access our servers, so I started to hack out a replacement for a piece of Standard ML code that was bugging me.

The reason I wanted to replace it was because I wanted something that was more maintainable and smaller.  Having a glaring dependency on an entire compiler suite for a relatively small job didn't make sense to me, and neither did parsing large matricies from the original binary into the rest of our software. I didn't write it, and nobody left at my work wrote it, and none of us (sorry!) are comfortable with using OCaML, or any ML derivatives. 

Essentially, it was a "black box" with a very simple purpose, and so it never got replaced. Now that I had a bit of freetime, I could pick away at it.

## Problem Scope ##
The original program 'count-kmers' counted k-mers up in a file full of DNA sequences (each base is denoted as a A, C, G, T). A k-mer is a combination of bases of size k. So, a 6-mer for example would be 'ACGAAC'. 

If you wanted to know how many ACGAAC's in the file, you'd have to do an exhaustive search. Here's an easy way to do that:

    calvin@barnabas:~/> grep -c 'ACGAAC' ~/sample/sample-13959.fasta 

Now this is ignoring headers and other problems you might encounter with the FASTA format.


That's not difficult, but what if we wanted it for each kmer? AAA AAC AAG AAT and so on. That's what count-kmers did. ( i added the k-mer comments)

    calvin@barnabas:~/quikr/src/nbc/> ./count kmers -r 6 -1 -U ./test.fasta | head -n 5
    2462 # AAAAAA
    169123 # AAAAAC
    1264616 # AAAAAG
    38421 # AAAAAT
    82301 # AAAACA
    
    
From here on out I will be using dna sequences from the EMP dataset "study960.fna" as my reference file. The time it took to count a 2.0Gb fasta file, on a standard Intel i5. (fasta files have a header for every sequence) with the original program i'll refer to as the epoch.

    time count-kmers -r 6 -1 -u ~/emp/data/sequences/study_940_split_library_seqs.fna > /dev/null

    real	0m24.801s
    user	0m24.098s
    sys     0m0.624s

 



## First Attempt ##
So the problem is small, but CPU intensive so I gave it a try. My first attempt performed pretty terribly, being 40 seconds slower.

The basic format of my program was to iterate through the lines in a file using getline, convert each kmer I saw into a base-10 index so that convert_mmer_to_index'ing AAA would return 0  and AAC would return 1, and adding 1 to that index in our counting array.

```C
    
// snippet from my main function.

while ((read = getline(&line, &len, fh)) != -1) {
    if(line[0] != '>')   {
        for(i = 0; i < strlen(line) - kmer; i++) {
            counts[convert_kmer_to_index(&line[i],kmer, width)]++;
        }
    }
}

    
long convert_kmer_to_index(char *str, long kmer, long error_pos) {

  char **ptr = NULL;
  char vals[kmer];
  long i = 0;
  long ret = 0;

  for(i = 0; i < kmer; i++) {
    int val = 0;
    switch(str[i]) {
      case 'a':
      case 'A':
        val = 48;
        break;
      case 'c':
      case 'C':
        val = 49;
        break;
      case 'g':
      case 'G':
        val = 50;
        break;
      case 't':
      case 'T':
        val = 51;
        break;
      default:
        return error_pos; 
    }


    vals[i] = val; 
  }

  ret = strtol(vals, ptr, 4);
  return ret;
}

```

```

time ./kmer_count_total ~/emp/data/sequences/study_940_split_library_seqs.fna 6 > /dev/null

real	1m4.811s
user	1m3.512s
sys	0m0.628s
```
    
Why was it so slow? Valgrind revealed some telling issues with my conversion function. Our switch statement was slow, and our strtol was taking up a bulk of our code.

Plenty of room for improvement!

## Attempt #2 ##

I am going to skip forward to commit #dc659709ec94c9d414df8b04f62322f74e1e4a0d.  Between that commit and the earlier one was a bit of clarification of code, a few declaration changes, but nothing that would redefine the nature of the code




