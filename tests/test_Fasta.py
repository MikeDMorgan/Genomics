import unittest
import hashlib
import os
import sys
import Fasta
import gzip


class TestFastaClass(unittest.TestCase):
    
    def test_string_filename(self):
        known_seq = ['GGGACAGGGGGC', 'GGGACTGGGGGGGC',
                     'ATGGCATAT', 'ATGGCATATCA',
                     'ATCGGAGGGATACGAG']
        ifile = os.path.join(os.getcwd(),
                             "_test_fasta.fa")
        _test_fasta = Fasta.Fasta(ifile)
        for _fas in _test_fasta:
            self.assertTrue(_fas.sequence in known_seq)

    def test_gzip_filename(self):
        known_seq = ['GGGACAGGGGGC', 'GGGACTGGGGGGGC',
                     'ATGGCATAT', 'ATGGCATATCA',
                     'ATCGGAGGGATACGAG']
        ifile = os.path.join(os.getcwd(),
                             "_test_fasta.fa.gz")
        _test_fasta = Fasta.Fasta(ifile)
        for _fas in _test_fasta:
            self.assertTrue(_fas.sequence in known_seq)

    def test_open_filename(self):
        known_seq = ['GGGACAGGGGGC', 'GGGACTGGGGGGGC',
                     'ATGGCATAT', 'ATGGCATATCA',
                     'ATCGGAGGGATACGAG']
        _file = os.path.join(os.getcwd(),
                             "_test_fasta.fa.gz")
        ifile = gzip.open(_file, 'rt')
        _test_fasta = Fasta.Fasta(ifile)
        for _fas in _test_fasta:
            self.assertTrue(_fas.sequence in known_seq)

    def test_seq_len(self):
        seq_lens = [12, 14, 9, 11, 16]
        _file = os.path.join(os.getcwd(),
                             "_test_fasta.fa.gz")
        ifile = gzip.open(_file, 'rt')
        _test_fasta = Fasta.Fasta(ifile)
        for _fas in _test_fasta:
            self.assertTrue(_fas.length in seq_lens)

    def test_gc_content(self):
        known_seq = ['GGGACAGGGGGC', 'GGGACTGGGGGGGC',
                     'ATGGCATAT', 'ATGGCATATCA',
                     'ATCGGAGGGATACGAG']
        gc_ = [0.833333333333, 0.857142857143, 0.333333333333,
               0.363636363636, 0.5625]
        gc_dict = dict(zip(known_seq, gc_))
        ifile = os.path.join(os.getcwd(),
                             "_test_fasta.fa.gz")
        _test_fasta = Fasta.Fasta(ifile)
        for _fas in _test_fasta:
            self.assertAlmostEqual(_fas.gc,
                                   gc_dict[_fas.sequence],
                                   places=4)

        


suite = unittest.TestLoader().loadTestsFromTestCase(TestFastaClass)
unittest.TextTestRunner(verbosity=2).run(suite)
