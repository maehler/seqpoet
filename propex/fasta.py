"""The FASTA module contains classes for accessing FASTA files and
FASTA index files.
"""

import collections
import itertools
import os

class FastaIndex(object):
    """Represents an index for a FASTA file."""

    def __init__(self, fname):
        """FastaIndex constructor.

        Args:
            fname: filename of the FASTA index
        """
        self.filename = fname
        if not os.path.isfile(fname):
            self.index = self.create_index()
        else:
            self.index = self.parse_index()

    def parse_index(self):
        """Parse a FASTA index file

        Parse a FASTA index file to match the format of samtools faidx.

        Returns:
            an OrderedDict with sequence names (headers) as keys and
            dicts as values. The value dicts have the following members:

            name:    the sequence name (FASTA header line)
            length:  sequence length
            offset:  the byte offset of the first base of the sequence
            nbase:   number of bases per line of sequence
            linelen: number of bytes per line of sequence
        Raises:
            ValueError: if the file cannot be parsed, if the file contains
                        duplicated headers or if the file is empty.
        """
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
        """Generate a FASTA index from a FASTA file.

        This function assumes that the fasta index filename in the constructor
        has a corresponding FASTA file without the ".fai" (or any other)
        extension.

        Returns:
            see parse_index
        Raises:
            ValueError: if the FASTA file contains duplicate headers or
                        if the FASTA file has sequence entries where lines
                        have different lengths (not counting the last line).
        """
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
    """Represent a FASTA record."""

    def __init__(self, seq, header):
        """FastaRecord constructor.

        Args:
            seq: the sequence of the record
            header: the name (or FASTA header) of the record
        """
        self.seq = seq
        self.header = header

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        return '<FastaRecord \'{0}\': {1}... ({2} nt)>'.format(self.header,
            self.seq[:5], len(self))

class Fasta(object):
    """Represent a FASTA file."""

    def __init__(self, fname):
        """Fasta constructor.

        Args:
            fname: filename of the FASTA file
        """
        self.filename = fname
        self.index = FastaIndex(self.filename + '.fai')

    def get_record(self, key):
        """Get a single FASTA record.

        Args:
            key: an integer.
        Returns:
            the FastaRecord stored at key.
        """
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
        """FastaRecord generator.

        Returns:
            a FastaRecord generator.
        """
        for i in xrange(len(self.index)):
            yield self.get_record(i)

    def __getitem__(self, key):
        return self.get_record(key)

    def __iter__(self):
        return self.generate_records()

    def __len__(self):
        return len(self.index)
