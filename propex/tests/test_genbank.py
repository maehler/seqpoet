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

class TestGenBankFeature:

    def test_qualifier_names(self):
        f = {'name': 'lalala'}
        gbf = propex.GenBankFeature('CDS', '123..679', f)
        assert gbf.get_qualifier('name') == f['name'], \
            'wrong name: {0}'.format(gbf.get_qualifier('name'))

    def test_parse_feature(self):
        feature = '''     CDS             complement(52625..53704)
                     /gene="recF"
                     /locus_tag="LMG718_02589"
                     /inference="ab initio prediction:Prodigal:2.60"
                     /inference="similar to AA sequence:UniProtKB:Q9RVE0"
                     /codon_start=1
                     /transl_table=11
                     /product="DNA replication and repair protein RecF"
                     /translation="MKLKQIELKNFRNYEDLKLDFHPNLNIFLGQNAQGKTNILEAIH
                     FLALTRSHRTSHDKELICWSGQEMKVSGLVEKAHVNVPLEVQLSSKGRIAKANHLKEN
                     RLADYIGQLKILMFAPENLELVKGSPATRRRFMDIELGQIHAVYLYDSMRYNRALKER
                     NAYLKFDQAKIDKNFLTVLDEQLAEHGNKIMFERKTFIEKLEIHAKKIHEQLTHGLET
                     LKITYNQNVKTDFSKELLSRQDHDIFRHQTTVGPHRDDLQFFINEINVADFGSQGQQR
                     TVTLSIKLAEIDLIFEETGEYPILLLDDVMSELDNHRQLDLIETSLGKTQTFITTTTL
                     DHLKNLPENLSIFHVTDGTIEKEKE"'''
        gbf = propex.GenBankFeature.from_string(feature)

        assert gbf.feature_type == 'CDS'
        gbf_gene = gbf.get_qualifier('gene')
        assert gbf_gene == 'recF', 'gene name ({0}) is wrong'.format(gbf_gene)
        gbf_inference = gbf.get_qualifier('inference')
        assert len(gbf_inference) == 2
        assert gbf_inference == ['ab initio prediction:Prodigal:2.60',
                                 'similar to AA sequence:UniProtKB:Q9RVE0']

    def test_empty_qualifiers(self):
        feature = '''     CDS             complement(52625..53704)
                     /gene="recF"
                     /locus_tag=
                     /note
                     /random=""'''
        gbf = propex.GenBankFeature.from_string(feature)

        assert gbf.get_qualifier('locus_tag') == ''
        assert gbf.get_qualifier('note') is None
        assert gbf.get_qualifier('random') == ''
