#-*- encoding: utf-8 -*-
"""Classes for dealing with GenBank files.

.. module:: genbank
.. moduleauthor:: Niklas Mähler <niklas.mahler@gmail.com>
"""

import collections
import re

from propex.sequence import Sequence

class Location(object):

    """Represent a GenBank feature location.

    For more information on locations, see
    http://www.insdc.org/files/feature_table.html#3.4

    Location (and GenBank files) are using 1-based positions. To make
    location-bases string handling easier, the Location class represents
    the locations internally as 0-based::

        >>> loc = Location('42..84')
        >>> loc.start
        41
        >>> loc.end
        83
        >>> loc.locstring
        '42..84'

    As of now, not all types of locations are implemented. Those
    that are implemented are single bases, ranges, ranges with
    unknown lower bound, ranges with unknown upper bound and
    locations where the exact position is unknown, but it is one
    of the bases between two positions.

    **Class attributes:**
        * **locstring:** the string representation of the location.
        * **loctype:** the type of the location.
        * **start:** the start position (0-based, including).
        * **end:** the end position (0-based, including).
        * **is_complement:** boolean indicating whether the position represents
          the complement of the sequence.

    :param locstring: a GenBank location string.
    """

    #: Regular expression for finding complement locations.
    loc_complement = re.compile(r'^complement\((.+)\)$')

    #: Regular expression for single base locations.
    loc_single = re.compile(r'^(\d+)$')
    #: Regular expression for range locations,
    loc_range = re.compile(r'^(\d+)\.\.(\d+)$')
    #: Regular expression for locations with unknown lower boundary.
    loc_lower_unknown = re.compile(r'^<(\d+)\.\.(\d+)$')
    #: Regular expression for locations with unknown upper boundary.
    loc_upper_unknown = re.compile(r'^(\d+)\.\.>(\d+)$')
    #: Regular expression for single base locations within a range.
    loc_one_of = re.compile(r'^(\d+)\.(\d+)$')

    def __init__(self, locstring):
        """Location constructor.

        :param locstring: a GenBank location string.
        """
        self.locstring = locstring
        self.loctype, self.start, self.end, self.is_complement = self._parse()

    def _regex_dict(self):
        """Utility function for location regular expressions.

        Returns:
            a dictionary where the keys correspond to the type of
            location and the values are the correpsponding regular
            expression.
        """
        return {
            'single': Location.loc_single,
            'range': Location.loc_range,
            'upper_unknown': Location.loc_lower_unknown,
            'lower_unknown': Location.loc_upper_unknown,
            'one_of': Location.loc_one_of
        }

    def _parse(self):
        """Parse a location string.

        Returns:
            a 4-tuple with the location type, start position, end
            position and a boolean to indicate whether the feature
            is located on the complement strand. Returned positions
            are 0-based.
        Raises:
            ValueError: if the location string is not valid.
        """
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
            start = end = int(regex.match(locstring).group(1))
        else:
            start, end = map(int, regex.match(locstring).groups())

        return re_name, start - 1, end - 1, is_complement

    def overlaps(self, location):
        """Test whether the location overlaps with another location.

        :param location: a Location object.
        :returns: True if the locations overlap with at least one base,
            otherwise False.
        """
        return self.start <= location.end and location.start <= self.end

    def min_distance(self, other):
        if self.overlaps(other):
            return 0
        else:
            return min(abs(self.start - other.end), abs(self.end - other.start))

    @classmethod
    def from_int(cls, start, end=None, strand='+'):
        if end is None:
            locstring = '{0}'.format(start)
        else:
            locstring = '{0}..{1}'.format(start, end)
        if strand == '+':
            return cls(locstring)
        else:
            return cls('complement({0})'.format(locstring))

    def __str__(self):
        return self.locstring

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return '<Location: {0}>'.format(repr(self.locstring))

class GenBankFeature(object):

    """Represent a GenBank feature.

    **Class attributes:**
        * **feature_type**: a string with the feature key.
        * **location**: a Location object representing the location of the feature.
        * **qualifiers**: a dictionary of qualifiers of the feature.

    :param locus: the name of the locus that the feature belongs to.
    :param feature_type: name of the feature.
    :param location: a Location object.
    :param qualifiers: a dictionary of qualifiers with the qualifier names as
                       keys and the qualifier values as values.
    """

    def __init__(self, locus, feature_type, location, qualifiers=None):
        """GenBankFeature constructor.

        Args:
            locus: the locus that the feature belongs to.
            feature_type: the key of the feature, e.g. 'CDS' or 'tRNA'.
            location: a Location object representing the location of the feature
            qualifiers: a dictionary of qualifiers with the qualifier names
                        as keys and the qualifier values as values.
        """
        self.locus = locus
        self.feature_type = feature_type
        self.location = location
        if qualifiers is None:
            self.qualifiers = []
        else:
            self.qualifiers = qualifiers

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self == other

    @classmethod
    def from_string(cls, locus, feature_string):
        """Create a GenBankFeature instance from a string.

        :param feature_string: a string representing a GenBank feature.
        :returns: a GenBankFeature object.
        """
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

        return cls(locus, ftype, Location(location), dict(qualifiers))

    def get_qualifier(self, qualifier_name):
        """Get a feature qualifier.

        :param qualifier_name: a string representing a qualifier.
        :returns: the value of the qualifier.
        :raises: KeyError if the feature does not have a qualifier called
                      ``qualifier_name``.
        """
        if qualifier_name not in self.qualifiers:
            raise KeyError('{0} is not a qualifier for {1}'
                .format(qualifier_name, self))
        return self.qualifiers[qualifier_name]

class GenBankLocus(object):

    """Represent a GenBank locus.

    **Class attributes:**
        * **name:** locus name.
        * **seq:** a Sequence object with the sequence of the locus.
        * **features:** a dictionary containing the features of the locus.

    :param name: the name of the locus.
    :param seq: a Sequence object representing the sequence of the locus.
    :param features: a dictionary containing features of the locus.
    """

    def __init__(self, name, seq, features=None):
        """GenBankLocus constructor.

        Args:
            name: the name of the locus.
            seq: a Sequence object representing the sequence of the locus.
            features: a dictionary containing features of the locus.
        """
        self.name = name
        self.seq = seq
        if features is None:
            self.features = {}
        else:
            self.features = features

    def features_at_location(self, location):
        """Get features at a location.

        :param location: a Location object.
        :returns: a list of GenBankFeature objects. If locus is not None, only
                  features in that locus are returned, otherwise all loci are
                  searched. Returns an empty list if there are no features
                  overlapping the location.
        """
        features = []
        for feat in GenBank.features:
            if feat not in self.features:
                continue
            for feature in self.features[feat]:
                if feature.location.overlaps(location):
                    features.append(feature)
        return features

    def next_upstream(self, feature):
        """Get a neighboring feature upstream of ``feature``.

        The function will find the next feature upstream of
        ``feature``. If the feature is on the forward (coding)
        strand, the function will return the next feature of the
        same type towards the 5' end of the strand. If the
        feature is on the reverse (template) strand, the next
        feature of the same type towards the 3' end will be
        returned.

        If there is no feature upstream of ``feature``,
        ``None`` is returned.

        :param feature: a :py:class:`.GenBankFeature` object.
        :returns: a :py:class:`.GenBankFeature` object or ``None``
                  if no feature is found.
        """
        return self._neighbor(feature, downstream=False)

    def next_downstream(self, feature):
        """Get a neighboring feature downstream of ``feature``.

        The function will find the next feature downstream of
        ``feature``. If the feature is on the forward (coding)
        strand, the function will return the next feature of the
        same type towards the 3' end of the strand. If the
        feature is on the reverse (template) strand, the next
        feature of the same type towards the 5' end will be
        returned.

        If there is no feature downstream of ``feature``,
        ``None`` is returned.

        :param feature: a :py:class:`.GenBankFeature` object.
        :returns: a :py:class:`.GenBankFeature` object or ``None``
                  if no feature is found.
        """
        return self._neighbor(feature, downstream=True)

    def _neighbor(self, feature, downstream=True):
        is_complement = feature.location.is_complement
        ftype = feature.feature_type
        findex = None
        found = False

        if ftype not in self.features:
            return None

        for i, f in enumerate(self.features[ftype]):
            if f == feature:
                # If the feature is on the opposite strand,
                # the direction will be opposite. The same
                # is true if we want to look upstream on the
                # forward strand.
                if (is_complement and downstream) or \
                        (not is_complement and not downstream):
                    findex = i - 1
                else:
                    findex = i + 1
                break


        if findex is None or findex >= len(self.features[ftype]) or findex < 0:
            return None

        # Make sure the feature is on the same strand
        while self.features[ftype][findex] \
                .location.is_complement != is_complement:
            if (is_complement and downstream) or \
                    (not is_complement and not downstream):
                findex -= 1
            else:
                findex += 1

        return self.features[ftype][findex]

class GenBank(object):

    """Represent a GenBank file.

    Class attributes:
        * filename: the filename of the GenBank file.
        * index: a list of dictionaries representing an index of the file.

    :param fname: filename of the GenBank file.
    :raises: ValueError if parsing fails.
    """

    #: List of supported features.
    features = ['CDS']

    def __init__(self, fname):
        """GenBank constructor.

        Args:
            fname: filename of the GenBank file.
        """
        self.filename = fname
        self.index = self._index()

    def _index(self):
        """Create and index of a the GenBank object.

        Returns:
            a list of dictionaries where each element in the list
            represents a locus.
        """
        indexdicts = []
        with open(self.filename) as f:
            offset = 0
            for lineno, line in enumerate(f):
                if lineno == 0 and not line.strip().startswith('LOCUS'):
                    raise ValueError('does not look like a GenBank file: {0}' \
                        .format(self.filename))
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

        # Sort the CDS according to start position
        for s in indexdicts:
            if 'CDS' not in s:
                continue
            s['CDS'] = sorted(s['CDS'], key=lambda x: x['location'].start)

        return indexdicts

    def __getitem__(self, index):
        """Get a specific GenBankLocus object.

        :param index: the index of the wanted locus in the index.
        :returns: a GenBankLocus object.
        """
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
                    features['CDS'].append(
                        GenBankFeature.from_string(locus_index['name'],
                            cds_string))

            # Get the sequence
            f.seek(origin_offset)
            line = f.readline()
            seq = ''.join(line.strip().split()[1:])
            while line.strip() != '//':
                line = f.readline()
                seq += ''.join(line.strip().split()[1:])

        return GenBankLocus(locus_index['name'], Sequence(seq), features)

    def get_locus_from_name(self, name):
        """Get a specific GenBankLocus object from the locus name.

        Since loci in a single GenBank file can have the same name,
        this function returns a list of loci.

        :param name: the name of the locus to retrieve.
        """
        loci = []
        for i, locus in enumerate(self.index):
            if locus['name'] == name:
                loci.append(self[i])
        return loci

    def __iter__(self):
        """Iterate over the loci.
        """
        for i in xrange(len(self)):
            yield self[i]

    def __len__(self):
        return len(self.index)
