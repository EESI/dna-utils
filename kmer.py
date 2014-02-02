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

def load_kmer_counts_from_file(fh, kmer):
  '''
  '''
  width = (kmer ** 4) + 1
  libkmer.get_kmer_counts_from_file.restype = c.POINTER(c.c_ulonglong * width )

  counts = libkmer.get_kmer_counts_from_filename(fh, kmer);

  return counts

