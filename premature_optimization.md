= Premature Optimization =

Often you hear the ringing refrain that 'premature opimization is the root of all evil" in one form or another. In this ost I will take you down the path of Optimization being the root of all bordom.

During a vacation I recently had, I spent quite a bit of time in the airports and flying. At that time I was unable to access our servers, so I started to hack out a replacement for a piece of Standard ML code that I wanted to get rid of .

The reason I wanted to replace it was because I wanted something that was more maintainable and smaller.  Having a glaring dependency on an entire compiler suite for a relatively small job didn't make sense anymore, and neither did parsing large matricies from the original into the rest of our software. I didn't write it, and nobody left at my work wrote it, and none of us (sorry!) are comfortable with using OCaML, or any ML derivatives. Essentially, it was a "black box" with a very simple purpose, and so it never got replaced. Now that I had a bit of freetime, I could pick away at it.

The original program 'count-kmers' did that - It counted kmers up in a file full of DNA sequences (each base is denoted as a A, C, G, T or another). A kmer is a combination of bases. A 6mer for example would be 'ACGAAC'. Now if you wanted to know how many ACGAAC's in the file, you'd have to do an exhaustive search. That's not difficult, but we want it for each kmer, AAA AAC AAG AAT and so on.

    calvin@barnabas:~/quikr/src/nbc/> ./count kmers -r 6 -1 -U ./test.fasta | head -n 5
    2462
    169123
    1264616
    38421
    82301
    
And so on and so forth
