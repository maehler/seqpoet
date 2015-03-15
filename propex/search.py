#-*- encoding: utf-8 -*-
"""Functions for searching for sequence matches.

.. module:: search
.. moduleauthor:: Niklas Mähler <niklas.mahler@gmail.com>
"""

def search(needle, haystack, mismatches=0):
    """Search for the occurence of ``needle`` in ``haystack``.

    :param needle: string to search for.
    :param haystack: string to search in.
    :param mismatches: the maximum number of mismatches allowed.
    :returns: an integer list with the starting positions of the matches.
    """
    n = len(needle)
    matches = []
    for i in xrange(0, len(haystack) - n):
        hd = hamming_distance(needle, haystack[i:i + n])
        if hd <= mismatches:
            matches.append(i)
    return matches

def hamming_distance(s1, s2):
    """Calculate the Hamming distance between two strings.

    :param s1: string
    :param s2: string
    :returns: the Hamming distance between ``s1`` and ``s2``.
    :raises: ValueError if the strings are not the same length.
    """
    if len(s1) != len(s2):
        raise ValueError('strings are of different length')
    return sum(x != y for x, y in zip(s1, s2))
