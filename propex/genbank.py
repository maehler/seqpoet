import collections
import re

from propex.sequence import Sequence

class Location(object):

    loc_complement = re.compile(r'^complement\((.+)\)$')

    loc_single = re.compile(r'^(\d+)$')
    loc_range = re.compile(r'^(\d+)\.\.(\d+)$')
    loc_lower_unknown = re.compile(r'^<(\d+)\.\.(\d+)$')
    loc_upper_unknown = re.compile(r'^(\d+)\.\.>(\d+)$')
    loc_one_of = re.compile(r'^(\d+)\.(\d+)$')

    def __init__(self, locstring):
        self.locstring = locstring
        self.loctype, self.start, self.stop, self.is_complement = self._parse()

    def _regex_dict(self):
        return {
            'single': Location.loc_single,
            'range': Location.loc_range,
            'upper_unknown': Location.loc_lower_unknown,
            'lower_unknown': Location.loc_upper_unknown,
            'one_of': Location.loc_one_of
        }

    def _parse(self):
        locstring = self.locstring
        re_name = None
        regex = None
        is_complement = False
        if Location.loc_complement.match(locstring):
            is_complement = True
            locstring = Location.loc_complement.match(locstring).group(1)
        for name, r in self._regex_dict().iteritems():
            if r.match(locstring) is not None:
                re_name = name
                regex = r
        if re_name is None:
            raise ValueError('unknown location string: {0}'.format(self.locstring))

        if re_name == 'single':
            start = stop = int(regex.match(locstring).group(1))
        else:
            start, stop = map(int, regex.match(locstring).groups())

        return re_name, start, stop, is_complement

    def overlaps(self, loc2):
        return self.start <= loc2.stop and loc2.start <= self.stop

    def __str__(self):
        return self.locstring

    def __repr__(self):
        return '<Location: {0}>'.format(repr(self.locstring))

class GenBankFeature(object):

    def __init__(self, feature_type, location, qualifiers=None):
        self.feature_type = feature_type
        self.location = location
        self.qualifiers = qualifiers

    @classmethod
    def from_string(cls, feature_string):
        lines = [x.strip() for x in feature_string.splitlines()]
        ftype, location = lines[0].strip().split()

        qualifiers = []

        for line in lines[1:]:
            if line.startswith('/'):
                # New qualifier
                i = line.find('=')
                key = line[1:i]
                value = line[i + 1:].strip('"')
                if i == -1:
                    key = line[1:]
                    value = None
                elif not value:
                    value = ''

                if len(qualifiers) > 0 and key == qualifiers[-1][0]:
                    # Multiple qualifiers with the same key
                    qualifiers[-1] = (key, [qualifiers[-1][1], value])
                else:
                    qualifiers.append((key, value))
            else:
                # Continuation of qualifier
                key = qualifiers[-1][0]
                value = qualifiers[-1][1] + line.strip('"')
                qualifiers[-1] = (key, value)

        return cls(ftype, Location(location), dict(qualifiers))

    def get_qualifier(self, qualifier_name):
        if qualifier_name not in self.qualifiers:
            raise KeyError('{0} is not a feature'.format(qualifier_name))
        return self.qualifiers[qualifier_name]

class GenBankLocus(object):

    def __init__(self, name, seq, features=None):
        self.name = name
        self.seq = seq

class GenBank(object):

    def __init__(self, fname):
        self.filename = fname
        self.index = self.index()

    def index(self):
        indexdicts = []
        with open(self.filename) as f:
            offset = 0
            for lineno, line in enumerate(f):
                if line.strip().split()[0] == 'LOCUS':
                    current_locus = line.strip().split()[1]
                    indexdicts.append({})
                    indexdicts[-1]['name'] = current_locus
                    indexdicts[-1]['offset'] = offset
                if line.strip().split()[0] == 'CDS':
                    if 'CDS' not in indexdicts[-1]:
                        indexdicts[-1]['CDS'] = []
                    indexdicts[-1]['CDS'].append({
                        'offset': offset,
                        'location': Location(line.strip().split()[1])
                    })
                if line.strip().split()[0] == 'ORIGIN':
                    indexdicts[-1]['ORIGIN'] = offset + len(line)
                offset += len(line)
        return indexdicts

    def get_locus(self, index):
        locus_index = self.index[index]
        locus_offset = locus_index['offset']
        origin_offset = locus_index['ORIGIN']
        features = {'CDS': []}
        with open(self.filename) as f:
            # Get the CDSs
            if 'CDS' in locus_index:
                for cds in locus_index['CDS']:
                    f.seek(cds['offset'])
                    cds_string = f.readline()
                    line = f.readline()
                    while line[5] == ' ':
                        cds_string += line
                        line = f.readline()
                    features['CDS'].append(GenBankFeature.from_string(cds_string))

            # Get the sequence
            f.seek(origin_offset)
            line = f.readline()
            seq = ''.join(line.strip().split()[1:])
            while line.strip() != '//':
                line = f.readline()
                seq += ''.join(line.strip().split()[1:])

        return GenBankLocus(locus_index['name'], Sequence(seq), features)

    def __iter__(self):
        for i in xrange(len(self)):
            yield self.get_locus(i)

    def __len__(self):
        return len(self.index)
