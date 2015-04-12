from nose.tools import raises
from nose.plugins.skip import SkipTest
import os

import propex
from propex.genbank import Location

class TestGenBank:

    def setUp(self):
        self.testdir = os.path.dirname(__file__)
        self.genbankdir = os.path.join(os.path.expanduser('~'), 'Dropbox',
                                       'operon_extractor', 'data_genbank')

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
        assert len(gb[0].seq) == 1522

    def test_iteration(self):
        gb = propex.GenBank(self.lmg718)
        for locus in gb:
            pass

    def test_load_directory(self):
        gbs = [propex.GenBank(os.path.join(self.genbankdir, x)) \
            for x in os.listdir(self.genbankdir)]

    def test_features_at_location(self):
        gb = propex.GenBank(self.lmg718)
        locus = gb.get_locus_from_name('718_Contig_100_c')[0]
        f = locus.features_at_location(Location('800'))
        assert len(f) == 1, 'found {0} features, expected 1'.format(len(f))
        assert f[0].get_qualifier('locus_tag') == 'LMG718_00002'

        f = locus.features_at_location(Location('450..720'))
        assert len(f) == 1, 'found {0} features, expected 1'.format(len(f))
        assert f[0].get_qualifier('locus_tag') == 'LMG718_00001'

        locus = gb.get_locus_from_name('718_Contig_102_c')[0]
        f = locus.features_at_location(Location('8800..8900'))
        assert len(f) == 2, 'found {0} features, expected 2'.format(len(f))
        assert f[0].get_qualifier('locus_tag') == 'LMG718_00019'
        assert f[1].get_qualifier('locus_tag') == 'LMG718_00020'

    def test_get_locus_from_name(self):
        gb = propex.GenBank(self.lmg718)
        loci = gb.get_locus_from_name('718_Contig_106_c')
        assert len(loci) > 0
        assert len(loci[0].seq) == 8967

    @raises(propex.genbank.ParsingError)
    def test_parse_fasta(self):
        gb = propex.GenBank(os.path.join(self.genbankdir, '..', 'data_fasta',
            'LMG718-cremoris.fasta'))

    def test_next_downstream(self):
        gb = propex.GenBank(self.lmg718)
        locus = gb.get_locus_from_name('718_Contig_10_co')[0]
        gbf = locus.features_at_location(Location('1355'))[0]
        next = locus.next_downstream(gbf)
        assert str(next.location) == '2532..2819'

    def test_next_downstream_duplicate_loci(self):
        gb = propex.GenBank(self.lmga18)
        locus = gb.get_locus_from_name('LMGA18_Contig_10')[1]
        gbf = locus.features_at_location(Location('301'))[0]
        next = locus.next_downstream(gbf)
        assert str(next.location) == '3180..3404'

    def test_next_downstream_last(self):
        gb = propex.GenBank(self.lmg718)
        locus = gb.get_locus_from_name('718_Contig_102_c')[0]
        gbf = locus.features_at_location(Location('9765'))[0]
        next = locus.next_downstream(gbf)
        assert next is None

    def test_next_downstream_complement(self):
        gb = propex.GenBank(self.lmg718)
        locus = gb.get_locus_from_name('718_Contig_101_c')[0]
        gbf = locus.features_at_location(Location('7664'))[0]
        next = locus.next_downstream(gbf)
        assert str(next.location) == 'complement(7271..7543)'

        gbf = locus.features_at_location(Location('5718'))[0]
        next = locus.next_downstream(gbf)
        assert str(next.location) == 'complement(2752..5457)'

    def test_next_upstream(self):
        gb = propex.GenBank(self.lmg718)
        locus = gb.get_locus_from_name('718_Contig_106_c')[0]
        gbf = locus.features_at_location(Location('754'))[0]
        next = locus.next_upstream(gbf)
        assert str(next.location) == '58..747'

    def test_next_upstream_duplicate_loci(self):
        gb = propex.GenBank(self.lmga18)
        locus = gb.get_locus_from_name('LMGA18_Contig_10')[1]
        gbf = locus.features_at_location(Location('3180'))[0]
        next = locus.next_upstream(gbf)
        assert str(next.location) == '301..1245'

    def test_next_upstream_last(self):
        gb = propex.GenBank(self.lmg718)
        locus = gb.get_locus_from_name('718_Contig_106_c')[0]
        gbf = locus.features_at_location(Location('58'))[0]
        next = locus.next_upstream(gbf)
        assert next is None

    def test_next_upstream_complement(self):
        gb = propex.GenBank(self.lmg718)
        locus = gb.get_locus_from_name('718_Contig_106_c')[0]
        gbf = locus.features_at_location(Location('7161'))[0]
        next = locus.next_upstream(gbf)
        assert str(next.location) == 'complement(8696..8953)'

        gbf = locus.features_at_location(Location('1945'))[0]
        next = locus.next_upstream(gbf)
        assert str(next.location) == 'complement(5718..6176)'

class TestGenBankFeature:

    def test_qualifier_names(self):
        f = {'name': 'lalala'}
        gbf = propex.GenBankFeature('testlocus', 'CDS', '123..679', f)
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
        gbf = propex.GenBankFeature.from_string('testlocus', feature)

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
        gbf = propex.GenBankFeature.from_string('testlocus', feature)

        assert gbf.get_qualifier('locus_tag') == ''
        assert gbf.get_qualifier('note') is None
        assert gbf.get_qualifier('random') == ''

    @raises(KeyError)
    def test_missing_qualifier(self):
        feature = '''     CDS             complement(52625..53704)
                     /gene="recF"'''
        gbf = propex.GenBankFeature.from_string('testlocus', feature)
        gbf.get_qualifier('locus_tag')

    def test_empty_qualifiers(self):
        gbf = propex.GenBankFeature('testlocus', 'CDS', '123..679')
        assert isinstance(gbf.qualifiers, list)
        assert len(gbf.qualifiers) == 0

    def test_equality(self):
        gbf1 = propex.GenBankFeature('testlocus', 'CDS',
            Location('123..679'), {'name': 'randomname'})
        gbf2 = propex.GenBankFeature('testlocus', 'CDS',
            Location('123..679'), {'name': 'randomname'})
        gbf3 = propex.GenBankFeature('testlocus', 'CDS',
            Location('123..679'), {'name': 'otherrandomname'})
        gbf4 = propex.GenBankFeature('testlocus', 'CDS',
            Location('120..679'), {'name': 'randomname'})

        assert gbf1 == gbf2
        assert not gbf1 != gbf2
        assert gbf3 != gbf1 and gbf3 != gbf2
        assert gbf1 != gbf4 and not gbf1 == gbf4

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
        match = Location._re_complement.match(self.complement)
        assert match
        assert match.group(1) == '340..565'
        match = Location._re_complement.match(self.complement2)
        assert match
        assert match.group(1) == '467'
        match = Location._re_complement.match(self.one_of)
        assert match is None

    def test_range_regex(self):
        match = Location._re_range.match(self.range)
        assert Location._re_one_of.match(self.range) is None
        assert match
        assert match.group(1) == '340'
        assert match.group(2) == '565'

    def test_single_regex(self):
        match = Location._re_single.match(self.single)
        assert match
        assert match.group(1) == '467'

    def test_lower_unknown(self):
        match1 = Location._re_lower_unknown.match(self.lower_unknown)
        match2 = Location._re_lower_unknown.match(self.lower_unknown2)
        assert match1
        assert match1.group(1) == '345'
        assert match1.group(2) == '500'
        assert match2
        assert match2.group(1) == '1'
        assert match2.group(2) == '888'

    def test_upper_unknown(self):
        match = Location._re_upper_unknown.match(self.upper_unknown)
        assert match
        assert match.group(1) == '1'
        assert match.group(2) == '888'

    def test_one_of(self):
        match = Location._re_one_of.match(self.one_of)
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
        assert loc.start == 466
        assert loc.end == 466
        assert not loc.is_complement

    def test_range(self):
        loc = Location(self.range)
        assert loc.start == 339
        assert loc.end == 564
        assert not loc.is_complement

    def test_lower_unknown(self):
        loc = Location(self.lower_unknown)
        assert loc.start == 344
        assert loc.end == 499
        assert not loc.is_complement
        loc = Location(self.lower_unknown2)
        assert loc.start == 0
        assert loc.end == 887
        assert not loc.is_complement

    def test_upper_unknown(self):
        loc = Location(self.upper_unknown)
        assert loc.start == 0
        assert loc.end == 887
        assert not loc.is_complement

    def test_one_of(self):
        loc = Location(self.one_of)
        assert loc.start == 101
        assert loc.end == 109
        assert not loc.is_complement

    def test_complement(self):
        loc = Location(self.complement)
        assert loc.start == 339
        assert loc.end == 564
        assert loc.is_complement
        loc = Location(self.complement2)
        assert loc.start == 466
        assert loc.end == 466
        assert loc.is_complement

    def test_overlap(self):
        loc1 = Location('100..200')
        loc2 = Location('150..250')
        loc3 = Location('150')
        loc4 = Location('200..300')
        loc5 = Location('201..301')
        assert loc1.overlaps(loc2)
        assert loc1.overlaps(loc3)
        assert loc1.overlaps(loc4)
        assert not loc3.overlaps(loc4)
        assert not loc4.overlaps(loc3)
        assert not loc1.overlaps(loc5)

    @raises(propex.genbank.LocationError)
    def test_invalid_location(self):
        loc = Location('123..noloc')

    def test_equality(self):
        loc1 = Location('100..200')
        loc2 = Location('100..200')
        loc3 = Location('100..201')
        loc4 = Location('complement(100..200)')

        assert loc1 == loc2
        assert not loc1 != loc2
        assert loc1 != loc4
        assert not loc2 == loc4
        assert loc1 != loc3

    def test_min_distance(self):
        loc1 = Location('100..200')
        loc2 = Location('100..200')
        loc3 = Location('250..300')
        loc4 = Location('1..50')

        assert loc1.min_distance(loc2) == 0
        assert loc1.min_distance(loc3) == 50
        assert loc3.min_distance(loc4) == 200
        assert loc1.min_distance(loc4) == 50

    def test_from_int(self):
        assert str(Location.from_int(100)) == '100'
        assert str(Location.from_int(100, 200)) == '100..200'
        assert str(Location.from_int(100, 200, '-')) == 'complement(100..200)'
        assert str(Location.from_int(100, strand='-')) == 'complement(100)'
