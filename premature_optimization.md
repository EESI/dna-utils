# Premature Optimization #

Often you hear the ringing refrain that 'premature opimization is the root of all evil' in one form or another. In this post I will take you down the path of 'Optimization is the result of bordom.'

Here's the best I could do with 5 minutes of google search:
![dilbert](http://joshreads.com/images/07/09/i070903dilbert.png "Man Hours")

## How it all started ##
During a vacation I recently took, I spent quite a bit of time in the airports and flying. During that time I was unable to access our servers, so I started to hack out a replacement for a piece of Standard ML code that was bugging me.

The reason I wanted to replace it was because I wanted something that was more maintainable and smaller.  Having a glaring dependency on an entire compiler suite for a relatively small job didn't make sense to me, and neither did parsing large matricies from the original binary into the rest of our software. I didn't write it, and nobody left at my work wrote it, and none of us (sorry!) are comfortable with using OCaML, or any ML derivatives. 

Essentially, it was a "black box" with a very simple purpose, and so it never got replaced. Now that I had a bit of freetime, I could pick away at it.

## Problem Scope ##
The original program 'count-kmers' did that - It counted k-mers up in a file full of DNA sequences (each base is denoted as a A, C, G, T or another). A k-mer is a combination of bases of size k. So, a 6-mer for example would be 'ACGAAC'. 

If you wanted to know how many ACGAAC's in the file, you'd have to do an exhaustive search. Here's an easy way to do that:

    calvin@barnabas:~/> grep -c 'ACGAAC' ~/sample/sample-13959.fasta 


That's not difficult, but what if we wanted it for each kmer? AAA AAC AAG AAT and so on. That's what count-kmers did.

    calvin@barnabas:~/quikr/src/nbc/> ./count kmers -r 6 -1 -U ./test.fasta | head -n 5
    2462
    169123
    1264616
    38421
    82301
    
    
From here on out I will be using dna sequences from the EMP dataset "study960" as my reference file.

    
So the problem is small, so I took a try. My first attempt didn't perform terribly.
