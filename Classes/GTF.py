from collections import OrderedDict
import re
import gzip


class Entry(object):
    '''
    A class for representing a GFF/GTF entry
    
    Attributes
    ----------
    contig: str
        chromosome which the feature is mapped to
    
    start: int
        bp starting position 1-based [) half-open of feature
        
    end: int
        bp ending position 1-based [), i.e. does not include this actual position
        
    length: int
        1-based length of feature
        
    feature_type: str
        feature_type described by the ensembl bio_type
        
    source: str
        source from which the feature was derived/annotated
        
    strand: str
        strand on which the feature is oriented (+ or -).  If == . then
        the feature orientation is irrelevant
        
    gene_id: str
        gene identifier for features belonging to the same structural gene
        
    transcript_id: str
        transcrip identifier for features belonging to the same transcribed unit
        
    attributes: dict
        a dictionary of additional attributes as key, value pairs
    '''
    
    # gtf attributes are positional

    def __init__(self, entry):
        self.contig = entry[0]
        self.source = entry[1]
        self.feature_type = entry[2]
        self.start = int(entry[3])
        self.end = int(entry[4])
        self.strand = entry[6]
        
        attributes = entry[8].split(";")
        self.gene_id = attributes.pop(0).split(" ")[1].strip('"')
        try:
            self.transcript_id = [tx for tx in attributes if re.search("transcrip_id", tx)][0].strip('"')
        except IndexError:
            self.transcript_id = self.gene_id
            
        self.length = self.end - self.start
        self.attributes = OrderedDict()
        for q in attributes:
            # splitting on white space may result in additional component in list
            # where there is a leading white, will result in index error
            # single trailing whitespaces should be ignored
            if len(q.split(" ")) == 2:
                x, y = q.split(" ")
            elif len(q.split(" ")) > 2:
                x, y = q.split(" ")[1:3]
            else:
                pass
            self.attributes[x] = y


class GTF(object):
    '''
    A class for iterating over a GFF/GTF file
    '''
    
    def __init__(self, file_handle):
        self.data = []
        self.file_handle = file_handle
            
    def __iter__(self):
        if type(self.file_handle) == gzip.GzipFile:
            for x in self.file_handle.readlines():
                line = x.rstrip("\n")
                if line[0] != '#':
                    entry = Entry(line.split("\t"))
                    yield entry

        elif type(self.file_handle) == str:
            if self.file_handle[-2:] == "gz":
                with gzip.open(self.file_handle, "rt") as gfile:
                    for x in gfile.readlines():
                        line = x.rstrip("\n")
                        if line[0] != '#':
                            entry = Entry(line.split("\t"))
                            yield entry
            else:
                with open(selffile_handle, "r") as gfile:
                    for x in gfile.readlines():
                        line = x.rstrip("\n")
                        if line[0] != '#':
                            entry = Entry(line.split("\t"))
                            yield entry
