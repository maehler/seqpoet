#-*- encoding: utf-8 -*-
"""Classes for dealing with GenBank files.

.. module:: genbank
.. moduleauthor:: Niklas Mähler <niklas.mahler@gmail.com>
"""

import collections
import itertools
import re

from seqpoet.sequence import Sequence

class LocationError(Exception):
    pass

def parse_location(locstring):
    """Parse a location string and return a :py:class:`.Location`
    or :py:class:`.JoinLocation` object.

    :param locstring: a GenBank location string.
    :raises: :py:class:`.LocationError` if parsing fails.
    """
    if locstring.startswith('join') or \
            locstring.startswith('complement(join'):
        return JoinLocation(locstring)
    else:
        return Location(locstring)

class JoinLocation(object):

    """Represent a "join" GenBank feature location.

    For more information on locations, see
    http://www.insdc.org/files/feature_table.html#3.4

    For information on how locations work, see :py:class:`.Location`.

    :param locstring: a GenBank location string.
    :raises: :py:class:`.LocationError` if parsing fails.
    """

    def __init__(self, locstring):
        self.locstring = locstring
        self.loctype = 'join'
        self.locations, self.start, self.end, self.is_complement = \
            self._get_locations()

    def _get_locations(self):
        complement_wrap = False
        if self.locstring.startswith('complement'):
            compmatch = re.match(r'^complement\((join\(.+\))\)$', self.locstring)
            if compmatch is None:
                raise LocationError('invalid join location: {0}' \
                    .format(self.locstring))
            locstring = compmatch.group(1)
            complement_wrap = True
        else:
            locstring = self.locstring

        match = re.match(r'^join\((.+)\)$', locstring)
        if match is None:
            raise LocationError('invalid join location: {0}' \
                .format(self.locstring))
        locstring = match.group(1)
        locations = [Location(x.strip()) for x in locstring.split(',')]
        if len(set(x.is_complement for x in locations)) != 1:
            raise LocationError('joint location is located on both strands')
        start = min(x.start for x in locations)
        end = max(x.end for x in locations)
        if not complement_wrap:
            is_complement = locations[0].is_complement
        else:
            is_complement = True
        return locations, start, end, is_complement

    def overlaps(self, other):
        """Test whether the location overlaps with another location.

        :param other: a :py:class:`.Location` or :py:class:`.JoinLocation`
            object.
        :returns: True if the locations overlap with at least one base,
            otherwise False.
        """
        if isinstance(other, JoinLocation):
            for loc1, loc2 in itertools.product(self.locations,
                    other.locations):
                if loc1.overlaps(loc2):
                    return True
        else:
            for loc in self.locations:
                if loc.overlaps(other):
                    return True
        return False

    def min_distance(self, other):
        """Get the minimum distance to another location.

        :param other: a :py:class:`.Location` or :py:class:`.JoinLocation`
            object.
        :returns: the minimum distance between the locations.
        """
        min_acc = []
        if isinstance(other, JoinLocation):
            for loc1, loc2 in itertools.product(self.locations,
                    other.locations):
                d = loc1.min_distance(loc2)
                if d == 0:
                    return 0
                min_acc.append(d)
        else:
            for loc in self.locations:
                d = loc.min_distance(other)
                if d == 0:
                    return 0
                min_acc.append(d)
        return min(min_acc)

    def __str__(self):
        return self.locstring

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return '<JoinLocation: {0}>'.format(repr(self.locstring))

class Location(object):

    """Represent a GenBank feature location.

    For more information on locations, see
    http://www.insdc.org/files/feature_table.html#3.4

    Location (and GenBank files) are using 1-based positions. To make
    location-based string handling easier, the Location class represents
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

        - **locstring:** the string representation of the location.
        - **loctype:** the type of the location.
        - **start:** the start position (0-based, including).
        - **end:** the end position (0-based, including).
        - **is_complement:** boolean indicating whether the position represents
          the complement of the sequence.

    :param locstring: a GenBank location string.
    :raises: :py:class:`.LocationError` if parsing fails.
    """

    #: Regular expression for finding complement locations.
    _re_complement = re.compile(r'^complement\((.+)\)$')
    #: Regular expression for single base locations.
    _re_single = re.compile(r'^(\d+)$')
    #: Regular expression for range locations,
    _re_range = re.compile(r'^(\d+)\.\.(\d+)$')
    #: Regular expression for locations with unknown lower boundary.
    _re_lower_unknown = re.compile(r'^<(\d+)\.\.(\d+)$')
    #: Regular expression for locations with unknown upper boundary.
    _re_upper_unknown = re.compile(r'^(\d+)\.\.>(\d+)$')
    #: Regular expression for locations with unknown upper and lower boundary.
    _re_lower_upper_unknown = re.compile(r'^<(\d+)\.\.>(\d+)$')
    #: Regular expression for single base locations within a range.
    _re_one_of = re.compile(r'^(\d+)\.(\d+)$')

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
            location and the values are the corresponding regular
            expression.
        """
        return {
            'single': Location._re_single,
            'range': Location._re_range,
            'upper_unknown': Location._re_lower_unknown,
            'lower_unknown': Location._re_upper_unknown,
            'lower_upper_unkown': Location._re_lower_upper_unknown,
            'one_of': Location._re_one_of
        }

    def _parse(self):
        """Parse a location string.

        Returns:
            a 4-tuple with the location type, start position, end
            position and a boolean to indicate whether the feature
            is located on the complement strand. Returned positions
            are 0-based.
        Raises:
            LocationError: if the location string is not valid.
        """
        locstring = self.locstring
        re_name = None
        regex = None
        is_complement = False
        if Location._re_complement.match(locstring):
            is_complement = True
            locstring = Location._re_complement.match(locstring).group(1)
        for name, r in self._regex_dict().iteritems():
            if r.match(locstring) is not None:
                re_name = name
                regex = r
        if re_name is None:
            raise LocationError('unknown location string: {0}' \
                .format(self.locstring))

        if re_name == 'single':
            start = end = int(regex.match(locstring).group(1))
        else:
            start, end = map(int, regex.match(locstring).groups())

        return re_name, start - 1, end - 1, is_complement

    def overlaps(self, other):
        """Test whether the location overlaps with another location.

        :param other: a :py:class:`.Location` or :py:class:`.JoinLocation`
            object.
        :returns: True if the locations overlap with at least one base,
            otherwise False.
        """
        if isinstance(other, JoinLocation):
            return other.overlaps(self)
        return self.start <= other.end and other.start <= self.end

    def min_distance(self, other):
        """Get the minimum distance to another location.

        :param other: a :py:class:`.Location` or :py:class:`.JoinLocation`
            object.
        :returns: the minimum distance between the locations.
        """
        if isinstance(other, JoinLocation):
            return other.min_distance(self)
        if self.overlaps(other):
            return 0
        else:
            return min(abs(self.start - other.end),
                abs(self.end - other.start))

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

        - **feature_type**: a string with the feature key.
        - **location**: a Location object representing the location of
          the feature.
        - **qualifiers**: a dictionary of qualifiers of the feature.

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
            location: a Location object representing the location of the
                      feature
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
        if len(lines) == 1:
            return cls(locus, ftype, parse_location(location), {})
        # Multiline location string
        i = 1
        line = lines[i]
        while not line.startswith('/'):
            i += 1
            location += line
            try:
                line = lines[i]
            except IndexError:
                break

        qualifiers = []

        for line in lines[i:]:
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
                    if isinstance(qualifiers[-1][1], list):
                        qualifiers[-1][1].append(value)
                    else:
                        qualifiers[-1] = (key, [qualifiers[-1][1], value])
                else:
                    qualifiers.append((key, value))
            else:
                # Continuation of qualifier
                key = qualifiers[-1][0]
                if isinstance(qualifiers[-1][1], list):
                    value = qualifiers[-1][1]
                    value[-1] += ' ' + line.strip('"')
                else:
                    value = qualifiers[-1][1] + ' ' + line.strip('"')
                qualifiers[-1] = (key, value)

        return cls(locus, ftype, parse_location(location), dict(qualifiers))

    def get_qualifier(self, qualifier_name):
        """Get a feature qualifier.

        :param qualifier_name: a string representing a qualifier.
        :returns: the value of the qualifier.
        :raises: :py:class:`KeyError` if the feature does not have a qualifier
                 called ``qualifier_name``.
        """
        if qualifier_name not in self.qualifiers:
            raise KeyError('{0} is not a qualifier for {1}'
                .format(qualifier_name, self))
        return self.qualifiers[qualifier_name]

    def __repr__(self):
        return '<GenBankFeature on {0} at {1}>'.format(self.locus, self.location)

class GenBankLocus(object):

    """Represent a GenBank locus.

    **Class attributes:**

        - **name:** locus name.
        - **seq:** a Sequence object with the sequence of the locus.
        - **features:** a dictionary containing the features of the locus.

    :param name: the name of the locus.
    :param seq: a Sequence object representing the sequence of the locus.
    :param features: a dictionary containing features of the locus.
    """

    def __init__(self, name, seq, features=None, header=None):
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
        if header is None:
            self.header = {}
        else:
            self.header = header

    def features_at_location(self, location):
        """Get features at a location.

        :param location: a Location object.
        :returns: a list of GenBankFeature objects. If locus is not None, only
                  features in that locus are returned, otherwise all loci are
                  searched. Returns an empty list if there are no features
                  overlapping the location.
        """
        features = []
        for feat in self.features.iterkeys():
            # Skip the source feature since it always will
            # overlap (assuming that the feature location
            # is within the sequence boundaries).
            if feat == 'source':
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


        if findex is None or findex >= len(self.features[ftype]) or \
                findex < 0:
            return None

        # Make sure the feature is on the same strand
        while self.features[ftype][findex] \
                .location.is_complement != is_complement:
            if (is_complement and downstream) or \
                    (not is_complement and not downstream):
                findex -= 1
            else:
                findex += 1
            if findex >= len(self.features[ftype]) or findex < 0:
                return None

        return self.features[ftype][findex]

class ParsingError(Exception):
    pass

class GenBank(object):

    """Represent a GenBank file.

    **Class attributes:**

        - filename: the filename of the GenBank file.
        - index: a list of dictionaries representing an index of the file.

    :param fname: filename of the GenBank file.
    :raises: :py:exc:`.ParsingError` if parsing fails.
    """

    def __init__(self, fname):
        """GenBank constructor.

        Args:
            fname: filename of the GenBank file.
        """
        self.filename = fname
        self.index, self.features = self._index()

    def _index(self):
        """Create and index of a the GenBank object.

        Returns:
            a list of dictionaries where each element in the list
            represents a locus.
        """
        features = set()
        indexdicts = []
        in_features = False
        with open(self.filename) as f:
            offset = 0
            for lineno, line in enumerate(f):
                if lineno == 0 and not line.strip().startswith('LOCUS'):
                    raise ParsingError(
                        'does not look like a GenBank file: {0}' \
                            .format(self.filename))
                if len(line.strip()) == 0:
                    offset += len(line)
                    continue
                if line.strip().split()[0] == 'LOCUS':
                    current_locus = line.strip().split()[1]
                    indexdicts.append({})
                    indexdicts[-1]['name'] = current_locus
                    indexdicts[-1]['offset'] = offset
                if line.strip().split()[0] == 'ORIGIN':
                    indexdicts[-1]['ORIGIN'] = offset + len(line)
                    in_features = False
                if in_features and line[5] != ' ':
                    feature = line.strip().split()[0]
                    features.add(feature)
                    if feature not in indexdicts[-1]:
                        indexdicts[-1][feature] = []

                    locstring = line.strip().split()[1]

                    nl = f.next()
                    loc_offset = offset + len(line)
                    while len(nl[:21].strip()) == 0 and not nl.strip().startswith('/'):
                        locstring += nl.strip()
                        loc_offset += len(nl)
                        nl = f.next()

                    f.seek(loc_offset)

                    indexdicts[-1][feature].append({
                        'offset': offset,
                        'location': parse_location(locstring)
                    })

                    offset = loc_offset
                    continue
                if line.startswith('FEATURES'):
                    in_features = True
                offset += len(line)

        # Sort the features according to start position
        for f in features:
            for s in indexdicts:
                if f not in s:
                    continue
                s[f] = sorted(s[f], key=lambda x: x['location'].start)

        return indexdicts, features

    def _parse_header(self, hstring):
        """Parse a GenBank header string into a nested dictionary.
        """
        head_data = collections.OrderedDict()

        header_lines = iter(hstring.splitlines(True))

        line = header_lines.next()

        header = line.strip().split()

        name = header[1]
        length = ' '.join(header[2:4])
        molecule = header[4]
        molecule_type = header[5]

        if len(header) == 8:
            division = header[6]
            date = header[7]
        elif len(header) == 7:
            division = ''
            date = header[6]

        head_data['LOCUS'] = {
            'name': name,
            'length': length,
            'molecule': molecule,
            'molecule_type': molecule_type,
            'genbank_division': division,
            'modification_data': date
        }

        last_key = None
        try:
            line = header_lines.next()
        except StopIteration:
            return head_data
        while True:
            if line[0] != ' ':
                key = line[:11].strip()
                last_key = key
                if key in head_data:
                    old_entry = head_data[key]
                    if not isinstance(old_entry, list):
                        head_data[key] = [(old_entry,
                            collections.OrderedDict())]
                    head_data[key].append((line[11:].strip(),
                        collections.OrderedDict()))
                else:
                    head_data[key] = line[11:].strip()
            elif len(line[:11].strip()) != 0:
                sub_key = line[:11].strip()
                old_entry = head_data[last_key]
                if not isinstance(old_entry, list):
                    head_data[last_key] = [(old_entry,
                        collections.OrderedDict())]
                old_entry = head_data[key][-1]
                old_entry[1][sub_key] = line[11:].strip()
                head_data[key][-1] = (old_entry[0], old_entry[1])
            else:
                if isinstance(head_data[last_key], list):
                    sub_key = head_data[last_key][-1][1].keys()[-1]
                    head_data[last_key][-1][1][sub_key] = \
                        '\n'.join([head_data[last_key][-1][1][sub_key],
                            line.strip()])
                else:
                    head_data[last_key] = '\n'.join([head_data[last_key],
                        line.strip()])
            try:
                line = header_lines.next()
            except StopIteration:
                break

        return head_data

    def __getitem__(self, index):
        """Get a specific GenBankLocus object.

        :param index: the index of the wanted locus in the index.
        :returns: a GenBankLocus object.
        """
        locus_index = self.index[index]
        locus_offset = locus_index['offset']
        origin_offset = locus_index['ORIGIN']
        features = collections.defaultdict(list)

        headstring = ''

        with open(self.filename) as f:
            f.seek(locus_offset)
            headstring += f.readline()

            line = f.readline()
            while not line.startswith('FEATURES'):
                headstring += line
                line = f.readline()

            head_data = self._parse_header(headstring)

            for ftype in self.features:
                if ftype not in locus_index:
                    continue
                for feature in locus_index[ftype]:
                    f.seek(feature['offset'])
                    feature_string = f.readline()
                    line = f.readline()
                    while len(line) < 6:
                        line = f.readline()
                    while line[5] == ' ':
                        feature_string += line
                        line = f.readline()
                        while len(line) < 6:
                            line = f.readline()
                    features[ftype].append(
                        GenBankFeature.from_string(locus_index['name'],
                            feature_string))

            # Get the sequence
            f.seek(origin_offset)
            line = f.readline()
            seq = ''.join(line.strip().split()[1:])
            while line.strip() != '//':
                line = f.readline()
                seq += ''.join(line.strip().split()[1:])

        return GenBankLocus(locus_index['name'], Sequence(seq), features,
            head_data)

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
