import imp
import os

from nose.plugins.skip import SkipTest

import seqpoet

currentdir = os.path.dirname(__file__)
rootdir = os.path.dirname(os.path.dirname(currentdir))
bindir = os.path.join(rootdir, 'bin')

seqpoet_script = imp.load_source('seqpoet_script',
	os.path.join(bindir, 'seqpoet'))

class TestFindOperon:

	def setup(self):
		gb_dir = os.path.expanduser('~/Dropbox/operon_extractor/data_genbank')
		gb_fname = os.path.join(gb_dir, 'LMG718-cremoris.gb')
		if not os.path.exists(gb_fname):
			raise SkipTest

		gb = seqpoet.genbank.GenBank(gb_fname)
		self.seqs = {
			gb_fname: gb
		}
		self.matches = [{
			'filename': gb_fname,
			'hitend': 3360,
			'hitstart': 3311,
			'length': 50,
			'seq': seqpoet.sequence.Sequence(
				'aattttactgatagctttttaaaaaataaaaaaaattactgacagaaatt'),
			'seqindex': 61,
			'seqname': '718_Contig_156_c',
			'strand': '+'
		}]
		self.minus_matches = [{
			'filename': gb_fname,
			'hitend': 3360,
			'hitstart': 3311,
			'length': 50,
			'seq': seqpoet.sequence.Sequence(
				'aatttctgtcagtaattttttttattttttaaaaagctatcagtaaaatt'),
			'seqindex': 61,
			'seqname': '718_Contig_156_c',
			'strand': '-'
		}]

	def test_operon_find(self):
		res = seqpoet_script.find_operon(self.matches, self.seqs,
			max_distance=500, extend_downstream=0, extend_upstream=0)
		assert len(res) == 0

	def test_operon_find_extend_upstream(self):
		res = seqpoet_script.find_operon(self.matches, self.seqs,
			max_distance=500, extend_downstream=0, extend_upstream=10)
		assert len(res) == 1, 'expected 1 result, got {0}'.format(len(res))
		assert len(res[0]['operon']) == 2

		operon_len = len(res[0]['seq'])
		assert operon_len == 3296, 'length is {0}'.format(operon_len)

	def test_operon_find_extend_downstream(self):
		res = seqpoet_script.find_operon(self.matches, self.seqs,
			max_distance=500, extend_downstream=100, extend_upstream=0)
		assert len(res) == 0, 'expected no results, got {0}'.format(len(res))

	def test_revcomp_operon_find(self):
		res = seqpoet_script.find_operon(self.matches, self.seqs,
			max_distance=500, extend_downstream=0, extend_upstream=0)
		assert len(res) == 0

	def test_revcomp_operon_find_extend_downstream(self):
		res = seqpoet_script.find_operon(self.minus_matches, self.seqs,
			max_distance=500, extend_downstream=100, extend_upstream=0)
		assert len(res) == 0

	def test_revcomp_operon_find_extend_upstream(self):
		res = seqpoet_script.find_operon(self.minus_matches, self.seqs,
			max_distance=500, extend_downstream=0, extend_upstream=100)
		assert len(res) == 1
		assert len(res[0]['operon']) == 2

		assert all(x.location.is_complement for x in res[0]['operon'])

		operon_len = len(res[0]['seq'])
		assert operon_len == 2135, 'length is {0}'.format(operon_len)
