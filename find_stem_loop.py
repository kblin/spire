#/usr/bin/env python
#
# SIRE - Search Prokaryote IRE sequences
# Find stem loop structures in a given RNA sequence using a sliding window
# Copyright (C) 2008-2010 Kai Blin <kai.blin@biotech.uni-tuebingen.de>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""find IRE stem loops using a sliding window"""

import sys

class SeqTooShortException(Exception):
    """exception thrown if sequence is too short to contain a stem loop"""
    pass

def can_pair(base_a, base_b):
    """check if two RNA bases can pair"""
    first = base_a.lower()
    second = base_b.lower()
    if first == 'a':
        return second == 'u'
    if first == 'c':
        return second == 'g'
    if first == 'g':
        if second in ('c', 'u'):
            return True
        return False
    if first == 'u':
        if second in ('a', 'g'):
            return True
        return False

def find_stem_loop(seq, stem_length=5, loop_length=5, loop_fuzzy=1, max_miss=0):
    """(seq, stem_length, loop_length, loop_fuzzy, max_miss)->loop list

    look for stem loops in a sequence.
    seq             is the sequence string to search
    stem_length     is the length of the stem in base pairs
    loop_length     is the length of the loop in bases
    loop_fuzzy      is the possible additional loop length
    max_miss        is the maximum number of allowed mismatches in the stem

    The function returns a list of stem loops, where a stem loop is a triple of
    (stem_start_position, stem_stop_position, stem_loop_sequence_string)

    If the sequence is too short to fold to a stem loop with the given
    parameters, a SeqTooShortException will be raised.

    """
    #pylint: disable-msg=R0914
    seq_len = len(seq)
    seq_offset = (stem_length * 2) + loop_length + loop_fuzzy
    if seq_len < seq_offset:
        #print >> sys.stderr, "seq_len is %s, seq_offset is %s" % (seq_len,
        #                                                          seq_offset)
        raise SeqTooShortException
    stem_loops = []

    i = 0
    #print "==== starting ===="
    while i < seq_len - seq_offset:
        left_start = i
        left_end = i + stem_length
        right_start = i + stem_length + loop_length
        right_end = i + stem_length + loop_length + loop_fuzzy + stem_length
        left_stem = seq[left_start:left_end][::-1]
        right_stem = seq[right_start:right_end]
        #print "%s: matching '%s' to '%s'" % (i, left_stem, right_stem)
        j = 0
        offset = 0
        mismatches = 0
        while j < stem_length:
            #print "%s = %s, ofs %s, miss %s" % (left_stem[j],
            #                                    right_stem[j + offset],
            #                                    offset, mismatches)
            if not can_pair(left_stem[j], right_stem[j + offset]):
                #print "mismatch"
                mismatches += 1
                if mismatches > max_miss and offset < loop_fuzzy:
                    #print "increasing offset"
                    offset += 1
                    j = 0
                    mismatches = 0
                    continue
            if mismatches > max_miss:
                #print "maximum number of mismatches reached," \
                #      " does not form a stem loop"
                break
            j += 1
        if mismatches <= max_miss:
            #print "found stem loop (%s, %s, '%s')" % (left_start, right_end,
            #                                       seq[left_start:right_end])
            stem_loops.append((left_start, right_end,
                               seq[left_start:right_end]))
        i += 1

    #print "==== done ===="
    return stem_loops

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print "Usage: %s <sequence> <stem_length> <loop_length> <loop_fuzzy>" \
              "<mismatches>" % sys.argv[0]
        sys.exit(1)
    print >> sys.stderr, sys.argv
    # results is _not_ a constant, no matter what pylint likes to believe
    #pylint: disable-msg=C0103
    results = find_stem_loop(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]),
                             int(sys.argv[4]), int(sys.argv[5]))
    for r in results:
        print "Stem loop at %s" % r
