from nose.tools import raises
from nose.plugins.skip import SkipTest
import os

import propex

class TestGenBank:

    def setUp(self):
        self.testdir = os.path.dirname(__file__)
        self.genbankdir = os.path.join(os.path.expanduser('~'), 'Dropbox',
                                       'operon_extractor', 'data')

        if not os.path.isdir(self.genbankdir):
            raise SkipTest

        self.lmg718 = os.path.join(self.genbankdir, 'LMG718-cremoris.gb')
        self.lmga18 = os.path.join(self.genbankdir, 'LMGA18-cremoris.gb')

    def test_index_length(self):
        gb = propex.GenBank(self.lmg718)
        assert len(gb) == 251, 'unexpected number of loci: {0}'.format(len(gb))

    def test_duplicate_locus_length(self):
        gb = propex.GenBank(self.lmga18)
        assert len(gb) == 231, 'unexpected number of loci: {0}'.format(len(gb))

    def test_sequence_length(self):
        gb = propex.GenBank(self.lmg718)
        assert len(gb.get_locus(0).seq) == 1522

    def test_iteration(self):
        gb = propex.GenBank(self.lmg718)
        for locus in gb:
            pass

    def test_load_directory(self):
        gbs = [propex.GenBank(os.path.join(self.genbankdir, x)) \
            for x in os.listdir(self.genbankdir)]
