'''

The kmer module provides an interface between the C-library
for quikr kmer counting in Python 

'''
__author__ = "Calvin Morrison"
__copyright__ = "Copyright 2013, EESI Laboratory"
__credits__ = ["Calvin Morrison"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Calvin Morrison"
__email__ = "mutantturkey@gmail.com"
__status__ = "development"

import numpy as np
import ctypes as c
try:
	libkmer = c.CDLL("libkmer.so"); 
except:
	raise Exception("Error: could not load libkmer.so")

def get_counts_from_file(fn, kmer):
	'''
	load kmer counts from a filename (not a handle), of size kmer.

returns an array size 4^k of counts

	>>> import kmer
	>>> counts = kmer.get_counts_from_file("test.fa", 6)

	'''

	if type(fn) is not str:
		raise TypeError("fn must be a str");
	if type(kmer) is not int:
		raise TypeError("fn must be int");

	ret = []
	width = (kmer ** 4)
	libkmer.get_kmer_counts_from_filename.restype = c.POINTER(c.c_ulonglong * width )
	counts = libkmer.get_kmer_counts_from_filename(fn, kmer);

	if counts:
		for i in counts.contents:
			ret.append(i)
	else:
		ret = 'error could not count mers on' + str(fn)

	return ret 

