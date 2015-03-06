from nose.tools import raises
import os

import propex

class TestFastaIndex:

    def setUp(self):
        testdir = os.path.dirname(__file__)
        self.valid_index = os.path.join(testdir, 'data', 'valid_index.fasta.fai')
        self.invalid_index = os.path.join(testdir, 'data', 'uneven.fasta.fai')
        self.empty_index = os.path.join(testdir, 'data', 'empty.fasta.fai')

    def test_faidx_length(self):
        faidx = propex.FastaIndex(self.valid_index)
        assert len(faidx) == 4

    def test_faidx_order(self):
        faidx = propex.FastaIndex(self.valid_index)
        assert faidx[0]['name'] == 'seq1'
        assert faidx[1]['name'] == 'seq2'
        assert faidx[2]['name'] == 'aaa'
        assert faidx[3]['name'] == 'bbb'

    def test_str(self):
        faidx = propex.FastaIndex(self.valid_index)
        assert str(faidx) == '\n'.join(['seq1\t78\t6\t28\t29',
                                        'seq2\t28\t93\t28\t29',
                                        'aaa\t44\t127\t28\t29',
                                        'bbb\t73\t178\t28\t29'])

    def test_keys(self):
        faidx = propex.FastaIndex(self.valid_index)
        assert faidx.keys() == ['seq1', 'seq2', 'aaa', 'bbb']

    def test_iter(self):
        faidx = propex.FastaIndex(self.valid_index)
        for k, v in faidx:
            pass

    def test_repr(self):
        faidx = propex.FastaIndex(self.valid_index)
        assert repr(faidx) == '<FastaIndex for {0}>' \
            .format(os.path.splitext(self.valid_index)[0])

    @raises(ValueError)
    def test_nonexisting_index(self):
        faidx = propex.FastaIndex(self.invalid_index)

    @raises(ValueError)
    def test_empty_index(self):
        faidx = propex.FastaIndex(self.empty_index)

    @raises(ValueError)
    def test_incorrect_filetype(self):
        faidx = propex.FastaIndex(os.path.splitext(self.invalid_index)[0])

class TestFasta:

    def setUp(self):
        testdir = os.path.dirname(__file__)
        self.valid_index = os.path.join(testdir, 'data', 'valid_index.fasta')
        self.threechrs_fname = os.path.join(testdir, 'data', 'uneven.fasta')
        self.dups_fname = os.path.join(testdir, 'data', 'dups.fasta')

    def test_fasta_length(self):
        fasta = propex.Fasta(self.valid_index)
        assert len(fasta) == 4, 'unexpected number of sequences'

    def test_fasta_headers(self):
        fasta = propex.Fasta(self.valid_index)
        headers = ['seq1', 'seq2', 'aaa', 'bbb']
        for i, record in enumerate(fasta):
            assert record.header == headers[i], \
                'header {0} is not {1}'.format(record.header, headers[i])
            assert len(record.seq.replace(' ', '')) == len(record.seq), \
                'spaces in sequence'

    def test_sequence_length(self):
        fasta = propex.Fasta(self.valid_index)
        lens = [78, 28, 44, 73]
        for i, record in enumerate(fasta):
            assert len(record) == lens[i], \
                'sequence length ({0}) is not {1}'.format(len(record), lens[i])

    def test_indexing(self):
        fasta = propex.Fasta(self.valid_index)
        assert len(fasta[1]) == 28
        assert fasta[1].seq == 'cacaggaggatagaccagatgacagata'
        assert repr(fasta[1]) == '<FastaRecord \'seq2\': cacag... (28 nt)>'

    @raises(IndexError)
    def test_invalid_index(self):
        fasta = propex.Fasta(self.valid_index)
        fasta[4]

    @raises(ValueError)
    def test_parse_duplicate_fasta(self):
        fasta = propex.Fasta(self.dups_fname)

class TestFastaWithoutIndex:

    def setUp(self):
        testdir = os.path.dirname(__file__)
        self.valid_noindex = os.path.join(testdir, 'data', 'valid_noindex.fasta')
        self.dups_noindex = os.path.join(testdir, 'data', 'dups_noindex.fasta')

    def tearDown(self):
        if os.path.isfile(self.valid_noindex + '.fai'):
            os.unlink(self.valid_noindex + '.fai')

    def test_fasta_length(self):
        fasta = propex.Fasta(self.valid_noindex)
        assert len(fasta) == 4, 'found {0} seqs, expected 4'.format(len(fasta))

    def test_fasta_headers(self):
        fasta = propex.Fasta(self.valid_noindex)
        headers = ['seq1', 'seq2', 'aaa', 'bbb']
        for i, record in enumerate(fasta):
            assert record.header == headers[i], \
                'header {0} is not {1}'.format(record.header, headers[i])
            assert len(record.seq.replace(' ', '')) == len(record.seq), \
                'spaces in sequence'

    def test_sequence_length(self):
        fasta = propex.Fasta(self.valid_noindex)
        lens = [78, 28, 44, 73]
        for i, record in enumerate(fasta):
            assert len(record) == lens[i], \
                'sequence length ({0}) is not {1}'.format(len(record), lens[i])

    @raises(ValueError)
    def test_duplicate_headers(self):
        fasta = propex.Fasta(self.dups_noindex)

class TestInvalidFasta:

    def setUp(self):
        testdir = os.path.dirname(__file__)
        self.empty_sequence = os.path.join(testdir, 'data', 'empty_sequence.fasta')
        self.uneven = os.path.join(testdir, 'data', 'uneven.fasta')

    def test_empty_sequence(self):
        fasta = propex.Fasta(self.empty_sequence)

    def test_fasta_headers(self):
        fasta = propex.Fasta(self.empty_sequence)
        headers = ['seq1', 'seq2', 'empty', 'aaa', 'bbb']
        for i, record in enumerate(fasta):
            assert record.header == headers[i], \
                'header {0} is not {1}'.format(record.header, headers[i])
            assert len(record.seq.replace(' ', '')) == len(record.seq), \
                'spaces in sequence'

    def test_sequence_length(self):
        fasta = propex.Fasta(self.empty_sequence)
        lens = [78, 28, 0, 44, 73]
        for i, record in enumerate(fasta):
            assert len(record) == lens[i], \
                'sequence length ({0}) is not {1}'.format(len(record), lens[i])

    @raises(ValueError)
    def test_uneven_rows(self):
        fasta = propex.Fasta(self.uneven)
