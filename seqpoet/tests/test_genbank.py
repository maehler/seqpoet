from nose.tools import raises
from nose.plugins.skip import SkipTest
import os
import tempfile

import seqpoet
from seqpoet.genbank import Location, JoinLocation

def temp_gbfile(gbstring):
    with tempfile.NamedTemporaryFile(delete=False) as temp:
        temp.write(gbstring)

    return temp.name

class TestGenBank:

    def setUp(self):
        self.testdir = os.path.dirname(__file__)
        self.sc = os.path.join(self.testdir, 'data', 'U49845.gb')
        self.gb = seqpoet.GenBank(self.sc)

    def test_sequence_length(self):
        assert len(self.gb[0].seq) == 5028

    def test_mRNA(self):
        assert len(self.gb[0].features['mRNA']) == 3

    def test_neighbors(self):
        locus = self.gb[0]
        gbf = locus.features['mRNA'][0]
        assert gbf is not None
        assert str(gbf.location) == '<1..>206'
        next = locus.next_downstream(gbf)
        assert next is not None
        assert str(next.location) == '<687..>3158'

        # Weird issue of alternating results when selecting next
        # downstream
        gbstring = '\n'.join(['LOCUS testlocus 5758 bp  DNA linear  12-APR-2015',
            'FEATURES            Location/qualifiers',
            '    source          1..5758',
            '    CDS             7..693',
            '    CDS             697..3303',
            '    CDS             complement(3381..4166)',
            '    CDS             complement(4167..5516)',
            'ORIGIN',
            '//'])

        gbfile = temp_gbfile(gbstring)

        gb = seqpoet.GenBank(gbfile)

        locus = gb[0]

        gbf = locus.features_at_location(Location('4170'))[0]
        assert gbf.location.is_complement
        assert str(gbf.location) == 'complement(4167..5516)'

        ds = locus.next_downstream(gbf)
        assert ds.location.is_complement
        assert str(ds.location) == 'complement(3381..4166)'

        ds = locus.next_downstream(ds)
        assert ds is None, 'should be None, found feature at {0}' \
            .format(ds.location)

        os.unlink(gbfile)

    def test_header(self):
        header = self.gb[0].header

        assert all(x in header for x in ['LOCUS', 'DEFINITION',
            'ACCESSION', 'VERSION', 'KEYWORDS', 'SOURCE', 'REFERENCE'])

        assert header['LOCUS']['molecule'] == 'DNA'

        assert header['ACCESSION'] == 'U49845'

        assert len(header['REFERENCE']) == 2
        assert header['REFERENCE'][0][0] == '1  (bases 1 to 5028)'
        assert all(x in header['REFERENCE'][0][1] for x in ['AUTHORS',
            'TITLE', 'JOURNAL', 'PUBMED'])

    def test_parse_header(self):
        '''Parse malformed locus line'''
        gbstring = '\n'.join(['LOCUS       NODE_18 673 bp   DNA linear',
            '03-FEB-2015',
            'FEATURES             Location/Qualifiers',
            'ORIGIN',
            '//'])

        gbfile = temp_gbfile(gbstring)

        gb = seqpoet.GenBank(gbfile)

        header = gb[0].header

        assert header['LOCUS']['name'] == 'NODE_18'
        assert header['LOCUS']['length'] == '673 bp'

        os.unlink(gbfile)

    @raises(seqpoet.genbank.ParsingError)
    def test_invalid_header(self):
        '''Invalid LOCUS line raises seqpoet.genbank.ParsingError'''
        gbstring = '\n'.join(['LOCUS       NODE_18 673',
            'FEATURES             Location/Qualifiers',
            'ORIGIN',
            '//'])

        gbfile = temp_gbfile(gbstring)

        gb = seqpoet.GenBank(gbfile)

        try:
            locus = gb[0]
        except seqpoet.genbank.ParsingError:
            raise
        finally:
            os.unlink(gbfile)

class TestGenBankLocal:

    def setUp(self):
        self.testdir = os.path.dirname(__file__)
        self.genbankdir = os.path.join(os.path.expanduser('~'), 'Dropbox',
                                       'operon_extractor', 'data_genbank')

        if not os.path.isdir(self.genbankdir):
            raise SkipTest

        self.lmg718 = os.path.join(self.genbankdir, 'LMG718-cremoris.gb')
        self.lmga18 = os.path.join(self.genbankdir, 'LMGA18-cremoris.gb')

    def test_index_length(self):
        gb = seqpoet.GenBank(self.lmg718)
        assert len(gb) == 251, 'unexpected number of loci: {0}'.format(len(gb))

    def test_duplicate_locus_length(self):
        gb = seqpoet.GenBank(self.lmga18)
        assert len(gb) == 231, 'unexpected number of loci: {0}'.format(len(gb))

    def test_sequence_length(self):
        gb = seqpoet.GenBank(self.lmg718)
        assert len(gb[0].seq) == 1522

    def test_iteration(self):
        gb = seqpoet.GenBank(self.lmg718)
        for locus in gb:
            pass

    def test_load_directory(self):
        gbs = [seqpoet.GenBank(os.path.join(self.genbankdir, x)) \
            for x in os.listdir(self.genbankdir) \
            if os.path.isfile(os.path.join(self.genbankdir, x))]

    def test_features_at_location(self):
        gb = seqpoet.GenBank(self.lmg718)
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
        gb = seqpoet.GenBank(self.lmg718)
        loci = gb.get_locus_from_name('718_Contig_106_c')
        assert len(loci) > 0
        assert len(loci[0].seq) == 8967

    @raises(seqpoet.genbank.ParsingError)
    def test_parse_fasta(self):
        gb = seqpoet.GenBank(os.path.join(self.genbankdir, '..', 'data_fasta',
            'LMG718-cremoris.fasta'))

    def test_next_downstream(self):
        gb = seqpoet.GenBank(self.lmg718)
        locus = gb.get_locus_from_name('718_Contig_10_co')[0]
        gbf = locus.features_at_location(Location('1355'))[0]
        next = locus.next_downstream(gbf)
        assert str(next.location) == '2532..2819'

    def test_next_downstream_duplicate_loci(self):
        gb = seqpoet.GenBank(self.lmga18)
        locus = gb.get_locus_from_name('LMGA18_Contig_10')[1]
        gbf = locus.features_at_location(Location('301'))[0]
        next = locus.next_downstream(gbf)
        assert str(next.location) == '3180..3404'

    def test_next_downstream_last(self):
        gb = seqpoet.GenBank(self.lmg718)
        locus = gb.get_locus_from_name('718_Contig_102_c')[0]
        gbf = locus.features_at_location(Location('9765'))[0]
        next = locus.next_downstream(gbf)
        assert next is None

    def test_next_downstream_complement(self):
        gb = seqpoet.GenBank(self.lmg718)
        locus = gb.get_locus_from_name('718_Contig_101_c')[0]
        gbf = locus.features_at_location(Location('7664'))[0]
        next = locus.next_downstream(gbf)
        assert str(next.location) == 'complement(7271..7543)'

        gbf = locus.features_at_location(Location('5718'))[0]
        next = locus.next_downstream(gbf)
        assert str(next.location) == 'complement(2752..5457)'

    def test_next_upstream(self):
        gb = seqpoet.GenBank(self.lmg718)
        locus = gb.get_locus_from_name('718_Contig_106_c')[0]
        gbf = locus.features_at_location(Location('754'))[0]
        next = locus.next_upstream(gbf)
        assert str(next.location) == '58..747'

    def test_next_upstream_duplicate_loci(self):
        gb = seqpoet.GenBank(self.lmga18)
        locus = gb.get_locus_from_name('LMGA18_Contig_10')[1]
        gbf = locus.features_at_location(Location('3180'))[0]
        next = locus.next_upstream(gbf)
        assert str(next.location) == '301..1245'

    def test_next_upstream_last(self):
        gb = seqpoet.GenBank(self.lmg718)
        locus = gb.get_locus_from_name('718_Contig_106_c')[0]
        gbf = locus.features_at_location(Location('58'))[0]
        next = locus.next_upstream(gbf)
        assert next is None

    def test_next_upstream_complement(self):
        gb = seqpoet.GenBank(self.lmg718)
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
        gbf = seqpoet.GenBankFeature('testlocus', 'CDS', '123..679', f)
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
        gbf = seqpoet.GenBankFeature.from_string('testlocus', feature)

        assert gbf.feature_type == 'CDS'
        gbf_gene = gbf.get_qualifier('gene')
        assert gbf_gene == 'recF', 'gene name ({0}) is wrong'.format(gbf_gene)
        gbf_inference = gbf.get_qualifier('inference')
        assert len(gbf_inference) == 2
        assert gbf_inference == ['ab initio prediction:Prodigal:2.60',
                                 'similar to AA sequence:UniProtKB:Q9RVE0']

    def test_multiple_qualifiers(self):
        feature = '''     ncRNA           476448..476561
                     /ncRNA_class="SRP_RNA"
                     /gene="ffs"
                     /locus_tag="b0455"
                     /gene_synonym="ECK0449"
                     /gene_synonym="JWR0009"
                     /product="4.5S sRNA component of Signal Recognition
                     Particle (SRP)"
                     /note="4.5S RNA; component of ribonucleoprotein particle;
                     works with the Ffh protein;
                     adjusted endpoints to reflect the mature 4.5S RNA (114
                     nt)"
                     /function="2.2.6 information transfer; RNA related; rRNA,
                     stable RNA"
                     /function="2.3.2 information transfer; protein related;
                     translation"
                     /function="7.1 location of gene products; cytoplasm"
                     /function="component of Signal Recognition Particle (SRP)
                     with the Ffh protein; involved in co-translational
                     targeting of proteins to membranes"
                     /function="RNA; Ribosomal and stable RNAs"
                     /db_xref="ASAP:ABE-0001579"
                     /db_xref="EcoGene:EG30027"'''
        gbf = seqpoet.GenBankFeature.from_string('testlocus', feature)
        func = gbf.get_qualifier('function')
        assert len(func) == 5

    def test_feature_join_location(self):
        feature = '''     CDS             join(52625..53704,54000..55000)
                     /gene="recF"
                     /locus_tag="LMG718_02589"
                     /inference="ab initio prediction:Prodigal:2.60"
                     /inference="similar to AA sequence:UniProtKB:Q9RVE0"
                     /codon_start=1
                     /transl_table=11'''
        gbf = seqpoet.GenBankFeature.from_string('testlocus', feature)

        assert gbf.feature_type == 'CDS'
        assert gbf.get_qualifier('gene') == 'recF'
        assert len(gbf.get_qualifier('inference')) == 2

    def test_multiline_location(self):
        feature = '''     CDS             complement(join(1294426..1294992,1294992..1295141,
                     1295140..1295322))
                     /gene="insZ"
                     /locus_tag="b4573"'''
        gbf = seqpoet.GenBankFeature.from_string('testlocus', feature)

        assert gbf.feature_type == 'CDS'
        assert isinstance(gbf.location, JoinLocation)
        assert len(gbf.location.locations) == 3
        assert gbf.get_qualifier('gene') == 'insZ'
        assert gbf.get_qualifier('locus_tag') == 'b4573'

        feature = '''     CDS             complement(join(1294426..1294992,
                     1294992..1295141,
                     1295140..1295322))
                     /gene="insZ"
                     /locus_tag="b4573"'''

        assert gbf.feature_type == 'CDS'
        assert isinstance(gbf.location, JoinLocation)
        assert len(gbf.location.locations) == 3
        assert gbf.get_qualifier('gene') == 'insZ'
        assert gbf.get_qualifier('locus_tag') == 'b4573'

    def test_empty_qualifiers(self):
        feature = '''     CDS             complement(52625..53704)
                     /gene="recF"
                     /locus_tag=
                     /note
                     /random=""'''
        gbf = seqpoet.GenBankFeature.from_string('testlocus', feature)

        assert gbf.get_qualifier('locus_tag') == ''
        assert gbf.get_qualifier('note') is None
        assert gbf.get_qualifier('random') == ''

    def test_minimal_feature(self):
        feature = '     CDS             complement(52625..53704)'
        gbf = seqpoet.GenBankFeature.from_string('testlocus', feature)

        assert gbf.feature_type == 'CDS'
        assert str(gbf.location) == 'complement(52625..53704)'

    @raises(KeyError)
    def test_missing_qualifier(self):
        feature = '''     CDS             complement(52625..53704)
                     /gene="recF"'''
        gbf = seqpoet.GenBankFeature.from_string('testlocus', feature)
        gbf.get_qualifier('locus_tag')

    def test_empty_qualifiers(self):
        gbf = seqpoet.GenBankFeature('testlocus', 'CDS', '123..679')
        assert isinstance(gbf.qualifiers, list)
        assert len(gbf.qualifiers) == 0

    def test_equality(self):
        gbf1 = seqpoet.GenBankFeature('testlocus', 'CDS',
            Location('123..679'), {'name': 'randomname'})
        gbf2 = seqpoet.GenBankFeature('testlocus', 'CDS',
            Location('123..679'), {'name': 'randomname'})
        gbf3 = seqpoet.GenBankFeature('testlocus', 'CDS',
            Location('123..679'), {'name': 'otherrandomname'})
        gbf4 = seqpoet.GenBankFeature('testlocus', 'CDS',
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
        self.lower_upper_unkown = '<1..>888'
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

    def test_lower_upper_unknown(self):
        match = Location._re_lower_upper_unknown.match(self.lower_upper_unkown)
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

    @raises(seqpoet.genbank.LocationError)
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

class TestJoinLocation:

    def setUp(self):
        self.jloc1 = JoinLocation('join(1..200,300..400)')
        self.jloc2 = JoinLocation('join(1..100, 200..300)')
        self.jloc3 = JoinLocation('join(150..175,180..190,310..320)')
        self.jloc4 = JoinLocation('join(complement(100),complement(200))')

    def test_instance(self):
        assert isinstance(self.jloc1, JoinLocation)
        assert isinstance(self.jloc2, JoinLocation)

    @raises(seqpoet.genbank.LocationError)
    def test_invalid_location(self):
        jloc = JoinLocation('join(1..200,300..400')

    @raises(seqpoet.genbank.LocationError)
    def test_invalid_strands(self):
        jloc = JoinLocation('join(complement(100),200)')

    @raises(seqpoet.genbank.LocationError)
    def test_messed_up_location(self):
        jloc = JoinLocation('complement(join(687..700,800..900,1000..1100))mRNA            <687..>3158')

    def test_start(self):
        # Remember, 0-indexed
        assert self.jloc1.start == 0
        assert self.jloc2.start == 0
        assert self.jloc3.start == 149

    def test_end(self):
        assert self.jloc1.end == 399
        assert self.jloc2.end == 299
        assert self.jloc3.end == 319

    def test_loctype(self):
        assert self.jloc1.loctype == 'join'
        assert self.jloc2.loctype == 'join'
        assert self.jloc3.loctype == 'join'

    def test_complement(self):
        assert not self.jloc1.is_complement
        assert self.jloc4.is_complement

    def test_complement_wrap(self):
        jloc = JoinLocation('complement(join(380844..381260,382591..382872))')
        assert isinstance(jloc, JoinLocation)
        assert jloc.is_complement
        assert jloc.start == 380843
        assert jloc.end == 382871

    def test_overlap(self):
        assert self.jloc1.overlaps(self.jloc2)
        assert self.jloc1.overlaps(Location('200..300'))
        assert self.jloc1.overlaps(Location('150..250'))
        assert not self.jloc3.overlaps(self.jloc2)
        assert self.jloc3.overlaps(self.jloc1)
        assert Location('150..250').overlaps(self.jloc1)

    def test_str(self):
        assert str(self.jloc1) == 'join(1..200,300..400)'
        assert str(self.jloc2) == 'join(1..100, 200..300)'

    def test_repr(self):
        assert repr(self.jloc1) == '<JoinLocation: \'join(1..200,300..400)\'>'
        assert repr(self.jloc2) == '<JoinLocation: \'join(1..100, 200..300)\'>'

    def test_min_distance(self):
        assert self.jloc1.min_distance(self.jloc2) == 0
        assert self.jloc1.min_distance(self.jloc3) == 0
        assert self.jloc2.min_distance(self.jloc3) == 10
        assert self.jloc1.min_distance(Location('250')) == 50
        assert Location('250').min_distance(self.jloc1) == 50
