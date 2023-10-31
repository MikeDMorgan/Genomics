import re
import gzip
import Genomics.Classes.Fasta as Fasta

class FastqRecord(Fasta.FastaRecord):
    '''
    A container for a single FASTQ record - an extension of FASTA

    Attributes
    ----
    qual: string
       The quality string giving the read quality
      
    '''

    def __init__(self, header, sequence, qual):
        self.header = header
        self.parseHeader()
        #self.header = header.rstrip("\n").lstrip("@")
        self.sequence = sequence
        self.length = self.len()
        self.gc = self._gcContent()
        self.qual = qual

    def parseHeader(self):
        '''
        Parse the header information for a single FASTQ record
        '''

        components = self.header.split(" ")
        self.machine, self.runid, self.flowcell, self.lane, self.tilenum, self.xcoord, self.ycoord = components[0].split(":")
        self.pair, self.filtered, self.bits, indices = components[1].split(":")
        self.R1idx, self.R2idx = indices.split("+")

    def record(self):
        '''
        Method to output a string format of the record to write to file
        '''

        return("{}\n{}\n+\n{}\n".format(self.header, self.sequence, self.qual))
        

        
class Fastq(Fasta.Fasta):
    '''
    Container class for FASTQ files, including
    methods for parsing.
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
                self.file_handle = open(file_handle, 'rt')

        self.iterator = self.parse(self.file_handle)
        
    def __iter__(self):
        return self

    def __next__(self):
        return next(self.iterator)

    def next(self):
        return next(self.iterator)

    def parse(self, infile):
        '''
        Parse the input FASTQ file, assign attributes from record 
        header, sequence and qual score

        Based on Fasta record
        '''

        ln = 1 # FASTQ files are in 4-line blocks
        
        seq = []
        is_qual = False
        for line in infile:
            if ln % 4 == 0:
                seq.append(line.rstrip("\n"))
                yield FastqRecord(seq[0], seq[1], seq[3])
                seq = []
            else:
                seq.append(line.rstrip("\n"))

            ln += 1
