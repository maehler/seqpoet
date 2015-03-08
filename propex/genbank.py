import collections
import weakref

class GenBankLocus(object):

    def __init__(self, name, seq):
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
                        'location': line.strip().split()[1]
                    })
                if line.strip().split()[0] == 'ORIGIN':
                    indexdicts[-1]['ORIGIN'] = offset + len(line)
                offset += len(line)
        return indexdicts

    def get_locus(self, index):
        locus_offset = self.index[index]['offset']
        origin_offset = self.index[index]['ORIGIN']
        with open(self.filename) as f:
            f.seek(origin_offset)
            line = f.readline()
            seq = ''.join(line.strip().split()[1:])
            while line.strip() != '//':
                line = f.readline()
                seq += ''.join(line.strip().split()[1:])

        return GenBankLocus(self.index[index]['name'], seq)

    def __iter__(self):
        for i in xrange(len(self)):
            yield self.get_locus(i)

    def __len__(self):
        return len(self.index)
