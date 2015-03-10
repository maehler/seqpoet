from nose.tools import raises
from nose.plugins.skip import SkipTest
import os

import propex
from propex.genbank import Location

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

    @raises(KeyError)
    def test_missing_qualifier(self):
        feature = '''     CDS             complement(52625..53704)
                     /gene="recF"'''
        gbf = propex.GenBankFeature.from_string(feature)
        gbf.get_qualifier('locus_tag')

class TestLocationRegex:

    def setUp(self):
        # Examples taken from http://www.insdc.org/files/feature_table.html#3.4.3
        self.single = '467'
        self.range = '340..565'
        self.lower_unknown = '<345..500'
        self.lower_unknown2 = '<1..888'
        self.upper_unknown = '1..>888'
        self.one_of = '102.110'
        self.complement = 'complement(340..565)'
        self.complement2 = 'complement(467)'

    def test_complement(self):
        match = Location.loc_complement.match(self.complement)
        assert match
        assert match.group(1) == '340..565'
        match = Location.loc_complement.match(self.complement2)
        assert match
        assert match.group(1) == '467'
        match = Location.loc_complement.match(self.one_of)
        assert match is None

    def test_range_regex(self):
        match = Location.loc_range.match(self.range)
        assert Location.loc_one_of.match(self.range) is None
        assert match
        assert match.group(1) == '340'
        assert match.group(2) == '565'

    def test_single_regex(self):
        match = Location.loc_single.match(self.single)
        assert match
        assert match.group(1) == '467'

    def test_lower_unknown(self):
        match1 = Location.loc_lower_unknown.match(self.lower_unknown)
        match2 = Location.loc_lower_unknown.match(self.lower_unknown2)
        assert match1
        assert match1.group(1) == '345'
        assert match1.group(2) == '500'
        assert match2
        assert match2.group(1) == '1'
        assert match2.group(2) == '888'

    def test_upper_unknown(self):
        match = Location.loc_upper_unknown.match(self.upper_unknown)
        assert match
        assert match.group(1) == '1'
        assert match.group(2) == '888'

    def test_one_of(self):
        match = Location.loc_one_of.match(self.one_of)
        assert match
        assert match.group(1) == '102'
        assert match.group(2) == '110'

class TestLocation:

    def setUp(self):
        # Examples taken from http://www.insdc.org/files/feature_table.html#3.4.3
        self.single = '467'
        self.range = '340..565'
        self.lower_unknown = '<345..500'
        self.lower_unknown2 = '<1..888'
        self.upper_unknown = '1..>888'
        self.one_of = '102.110'
        self.complement = 'complement(340..565)'
        self.complement2 = 'complement(467)'

    def test_single(self):
        loc = Location(self.single)
        assert loc.start == 467
        assert loc.stop == 467
        assert not loc.is_complement

    def test_range(self):
        loc = Location(self.range)
        assert loc.start == 340
        assert loc.stop == 565
        assert not loc.is_complement

    def test_lower_unknown(self):
        loc = Location(self.lower_unknown)
        assert loc.start == 345
        assert loc.stop == 500
        assert not loc.is_complement
        loc = Location(self.lower_unknown2)
        assert loc.start == 1
        assert loc.stop == 888
        assert not loc.is_complement

    def test_upper_unknown(self):
        loc = Location(self.upper_unknown)
        assert loc.start == 1
        assert loc.stop == 888
        assert not loc.is_complement

    def test_one_of(self):
        loc = Location(self.one_of)
        assert loc.start == 102
        assert loc.stop == 110
        assert not loc.is_complement

    def test_complement(self):
        loc = Location(self.complement)
        assert loc.start == 340
        assert loc.stop == 565
        assert loc.is_complement
        loc = Location(self.complement2)
        assert loc.start == 467
        assert loc.stop == 467
        assert loc.is_complement

    @raises(ValueError)
    def test_invalid_location(self):
        loc = Location('123..noloc')
