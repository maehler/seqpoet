import collections
import itertools
import os

class FastaIndex(object):

    def __init__(self, fname):
        self.filename = fname
        if not os.path.isfile(fname):
            self.index = self.create_index()
        else:
            self.index = self.parse_index()

    def parse_index(self):
        index = collections.OrderedDict()
        with open(self.filename) as f:
            for line in f:
                line = line.strip().split('\t')
                try:
                    length, offset, linenbase, linelen = map(int, line[1:])
                except:
                    raise ValueError('index file of incorrect format')
                name = line[0]
                if name in index:
                    raise ValueError('fasta contains duplicate headers')
                index[name] = {'name': name, 'length': length,
                    'offset': offset, 'nbase': linenbase, 'linelen': linelen}
        if len(index) == 0:
            raise ValueError('fasta index is empty')
        return index

    def create_index(self):
        fasta_fname = os.path.splitext(self.filename)[0]
        index = collections.OrderedDict()
        header_offset = 0
        with open(fasta_fname) as f:
            while True:
                line = f.readline()
                if not line:
                    break
                if len(line.strip()) == 0:
                    line = f.readline()
                if line[0] == '>':
                    header = line[1:].strip()
                    linelens = []
                    baselens = []
                    offset = f.tell()
                    line = f.readline()
                    while len(line) > 0 and line[0] != '>':
                        linelens.append(len(line))
                        baselens.append(len(line.strip()))
                        header_offset = f.tell()
                        line = f.readline()
                    if header in index:
                        raise ValueError('fasta contains duplicate headers')
                    if len(linelens) > 0:
                        if not all(x == linelens[0] for x in linelens[:-1]) \
                                or not all(x == baselens[0] for x in baselens[:-1]):
                            raise ValueError('fasta has lines of different '
                                             'length for the same sequence: '
                                             '{0}'.format(header))
                        index[header] = {
                            'name': header, 'length': sum(baselens),
                            'offset': offset,
                            'nbase': baselens[0],
                            'linelen': linelens[0]
                        }
                    else:
                        # In case of empty sequence
                        header_offset = offset
                        index[header] = {
                            'name': header, 'length': sum(baselens),
                            'offset': offset,
                            'nbase': 0,
                            'linelen': 0
                        }
                f.seek(header_offset)
        return index

    def keys(self):
        return self.index.keys()

    def __len__(self):
        return len(self.index)

    def __getitem__(self, key):
        return self.index[self.index.keys()[key]]

    def __iter__(self):
        return self.index.iteritems()

    def __repr__(self):
        return '<FastaIndex for {0}>'.format(os.path.splitext(self.filename)[0])

    def __str__(self):
        return '\n'.join('{name}\t{length}\t{offset}\t{nbase}\t{linelen}'
            .format(**x) for x in self.index.itervalues())

class FastaRecord(object):

    def __init__(self, seq, header):
        self.seq = seq
        self.header = header

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        return '<FastaRecord \'{0}\': {1}... ({2} nt)>'.format(self.header,
            self.seq[:5], len(self))

class Fasta(object):

    def __init__(self, fname):
        self.filename = fname
        self.index = FastaIndex(self.filename + '.fai')

    def get_record(self, key):
        indexdict = self.index[key]
        if indexdict['length'] == 0:
            return FastaRecord('', indexdict['name'])
        seq = ''
        with open(self.filename) as f:
            f.seek(indexdict['offset'])
            ncompletelines = indexdict['length'] // indexdict['nbase']
            ncompletebytes = ncompletelines * indexdict['linelen']
            for lineno in xrange(ncompletelines):
                seq += f.read(indexdict['linelen']).strip()
            restbytes = indexdict['length'] % len(seq)
            seq += f.read(restbytes).strip()
        return FastaRecord(seq, indexdict['name'])

    def generate_records(self):
        for i in xrange(len(self.index)):
            yield self.get_record(i)

    def __getitem__(self, key):
        return self.get_record(key)

    def __iter__(self):
        return self.generate_records()

    def __len__(self):
        return len(self.index)
