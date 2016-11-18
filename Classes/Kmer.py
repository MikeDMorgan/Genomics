'''
Functions for handling, generating and processing kmers
'''

import random
import itertools


def kmerGen(n):
    '''
    Generate all overlapping Kmers of length n
    from a sequence of length L
    '''
    
    alphabet = ['A', 'T', 'C', 'G']

    for each in itertools.product(alphabet, repeat=n):
        yield ''.join(each)

def randomKmer(L):
    '''
    Generate a random sequence of length L
    '''
    
    seq = []
    alphabet = ['A', 'T', 'C', 'G']
    while len(seq) < 50:
        idx = random.randint(0, 3)
        seq.append(alphabet[idx])
    
    return ''.join(seq)


def findKmers(seq, n):
    '''
    find all overlapping kmers of length n in sequence
    '''
    
    kmer_set = set()
    L = len(seq)
    if n > L:
        return seq
    for k in range(0, (L - n) + 1):
        kmer = seq[k:k+n]
        kmer_set.add(kmer)
        
    return kmer_set


