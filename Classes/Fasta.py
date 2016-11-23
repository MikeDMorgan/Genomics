import re
import gzip


class FastaRecord(object):
    '''
    Container for a single FASTA sequence record

    Attributes
    ----------
    header: str
        The header for a FASTA record, preceded by '>'
        
    sequence: str
        The FASTA record sequence
        
    gc: float
        The GC content of the FASTA sequence
    
    length: int
        The total number of nucleotides in the FASTA record
    '''
    
    def __init__(self, header, sequence):
        self.header = header.rstrip("\n").lstrip(">")
        self.sequence = sequence
        self.length = self.len()
        self.gc = self._gcContent()
        
    def __len__(self):
        return len(self.sequence)
    
    def len(self):
        return self.__len__()
    
    def _gcContent(self):
        # what is quicker? remove all non-GC then count these
        # or iterate over total sequence?
        gc = len([g for g in self.sequence if g.lower() in ('g', 'c')])
        return gc/float(self.length)


class Fasta(object):
    '''
    Container for fasta file object,
    methods for parsing fasta file
    '''
    
    def __init__(self, file_handle):
        self.seqs = []

        if type(file_handle) == gzip.GzipFile:
            self.file_handle = file_handle

        elif type(file_handle) == str:
            suff = file_handle[-2:]
            if suff in ['gz', '.z']:
                self.file_handle = gzip.open(file_handle, 'rt')
            else:
                self.file_handle = open(file_handle, 'r')
                    
        self.iterator = self.parse(self.file_handle)
    
    def __iter__(self):
        return self
    
    def __next__(self):
        return next(self.iterator)
    
    def next(self):
        return next(self.iterator)
            
    def parse(self, infile):
        '''
        Parse the input fasta file, assign attributes from
        record header and sequence content
        
        Code from CGAT.FastaIterator.py
        '''
        
        
        h = infile.readline()[:-1]
            
        if not h:
            infile.close()
            raise StopIteration
        # skip everything until first fasta entry starts
        while h[0] != ">":
            h = infile.readline()[:-1]
            if not h:
                infile.close()
                raise StopIteration
            continue

        h = h[1:]
        seq = []
        for line in infile.readlines():
            if line.startswith('#'):
                continue

            if line.startswith('>'):
                yield FastaRecord(line, ''.join(seq))
                    
                h = line[1:-1]
                seq = []
                continue
            
            seq.append(line[:-1])
        
        yield FastaRecord(h, ''.join(seq))
