#/usr/bin/env python
#
# SIRE - Search Prokaryote IRE sequences
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

import sys
sys.path.append('.')
sys.path.append('..')
from search_ire import *

successes = 0
failures = 0
total = 0

class FakeReMatch:
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def span(self):
        return (self.start, self.end)

def test_reverse_complement(name):
    seq =   "AGGCU"
    compl = "AGCCU"
    if reverse_complement(seq) != compl:
        print "%s: got '%s', expected '%s'" % (name, reverse_complement(seq), compl)
        return False
    if reverse_complement(compl) != seq:
        print "%s: got '%s', expected '%s'" % (name, reverse_complement(compl), seq)
        return False
    return True

def test_match(name):
    start = 11
    end = 21
    seq = "XXXXXCNNNNNCAGUGYMMMMZXXXXX"
    before = "NNNNN"
    after =  "YMMMMZ"
    dir = FORWARD
    m = Match(FakeReMatch(start, end), seq, dir)

    if m.start() != start:
        print "%s: m.start() is %d, expected %d" % (name, m.start(), start)
        return False
    if m.end() != end:
        print "%s: m.end() is %d, expected %d" % (name, m.end(), end)
        return False
    if m.get_before() != before:
        print "%s: m.get_before() is '%s', not '%s'" % (name, m.get_before(), before)
        return False
    if m.get_after() != after:
        print "%s: m.get_after() is '%s', not '%s'" % (name, m.get_after(), after)
        return False
    return True

def test_find_matches(name):
    #             Match without extra base            Match with extra base
    #                1         2         3         4         5         6         7
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNNCAGUGYMMMMZ                  CNNNNNCAGUGYMMMMZ
    seq = "AAAAAAAAAACAAGUUCAGUGAGCUUUUUUUUUUUGGGGGGGGGGCAGGCUCAGUGCAGCUUCCCCCCCC"
    #                                    1         1         1         1         1
    #      7         8         9         0         1         2         3         4
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNQCAGUGYMMMMZ                  ZMMMMYCACUGNNNNNG
    seq +="AAAAAAAAAACAAGUCCAGUGAGCUUUUUUUUUUUGGGGGGGGGGAGCAACCACUGUUGCUGCCCCCCCC"
    #                  should not fold               reverse match should fold
    #      1         1         1         1         1         1         2         2
    #      4         5         6         7         8         9         0         1
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNQCAGUGYMMMMZ                  ZMMMMYCACUGNNNNNG
    seq +="AAAGUUCUUGCUUCAACAGUGUUUGAACGGAACAAGGGGGGGGGGAGCAACCACUGUUGCUGCCCCCCCC"
    #             hferritin should fold              reverse match should fold
    #      2         2         2         2         2         2         2         2
    #      1         2         3         4         5         6         7         8
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNQGAGAGYMMMMZ                    CNNNNNUAGUAYMMMMZ
    seq +="AAAGUUCUUGCUUCAAGAGAGUUUGAACGGAACAAGGGGGGGGGGGGCAACUAUAGUACUGGUUCCCCCC"
    #             GAGAG pattern, should fold          UAGUA pattern, should fold


    forward_pattern = re.compile(FORWARD_PATTERN)
    reverse_pattern = re.compile(REVERSE_PATTERN)

    matches = find_matches(seq, forward_pattern, FORWARD)
    matches.extend(find_matches(seq, reverse_pattern, REVERSE))

    if not matches:
        print "%s: failed to find matches"
        return False
    if len(matches) != 8:
        print "%s: found %d matches, expected 8" % (name, len(matches))
        for match in matches:
            print match
        return False

    m1 = matches[0]
    m2 = matches[1]
    m5 = matches[3] #forward matches are found first
    m4 = matches[6]

    if m1.start() != 10:
        print "%s: First match reported at %d, not 10" % (name, m1.start())
        return False
    if m1.end() != 21:
        print "%s: First match ends at %d, not 21" % (name, m1.end())
        return False
    if m1.sequence != "AAAAACAAGUUCAGUGAGCUUUUUUU":
        print "%s: First match sequence is '%s', not 'AAAAACAAGUUCAGUGAGCUUUUUUU'" % (name, m1.sequence)
        return False
    if m2.start() != 45:
        print "%s: Second match reported at %d, not 45" % (name, m2.start())
        return False
    if m2.end() != 56:
        print "%s: Second match ends at %d, not 56" % (name, m2.end())
        return False
    if m2.sequence != "GGGGGCAGGCUCAGUGCAGCUUCCCC":
        print "%s: Second match sequence is '%s', not 'GGGGGCAGGCUCAGUGCAGCUUCCCC'" % (name, m2.sequence)
        return False
    if m4.direction != REVERSE:
        print "%s: Reverse match not detected as revese" % name
        return False
    if m4.start() != 121:
        print "%s: Reverse match reported at %d, not 121" % (name, m4.start())
        return False
    if m4.end() != 132:
        print "%s: Reverse match ends at %d, not 132" % (name, m4.end())
        return False
    if m4.sequence != reverse_complement("GGGGAGCAACCACUGUUGCUGCCCCC"):
        print "%s: Reverse match sequence is '%s', not '%s'" % (name, m4.sequence, reverse_complement("GGGGAGCAACCACUGUUGCUGCCCCC"))
        return False
    if len(m4.sequence) != len("GGGGAGCGACCACUGUUGCUGCCCCC"):
        print "%s: Reverse complement length is %s, not %s" % (name, len(m4.sequence), len("GGGGAGCGACCACUGUUGCUGCCCCC"))
        return False
    return True

def test_find_cds_features(name):
    print "%s: implement me!" % name
    return False

def test_is_upstream(name):
    print "%s: implement me!" % name
    return False

def test_is_downstream(name):
    print "%s: implement me!" % name
    return False

def test_is_within(name):
    print "%s: implement me!" % name
    return False

def test_can_pair(name):
    test_seq = "AcGGu"
    correct_seq = "UGcUA"
    wrong_seq = "CAaAC"

    for i in range(0,5):
        t = test_seq[i]
        c = correct_seq[i]
        w = wrong_seq[i]

        if not can_pair(t, c):
            print "%s: failed to pair %s=%s" % (name, t, c)
            return False
        if can_pair(t, w):
            print "%s: erroneously paired %s=%s" % (name, t, w)
            return False
    return True

def test_filter_fold(name):
    #             Match without extra base            Match with extra base
    #                1         2         3         4         5         6
    #      0123456789012345678901234567890123456789012345678901234567890123456789
    #                CNNNNNCAGUGYMMMMZ                  CNNNNNCAGUGYMMMMZ
    seq = "AAAAAAAAAACAAGUUCAGUGAGCUUUUUUUUUUUGGGGGGGGGGCAGGCUCAGUGCAGCUUCCCCCCCC"
    #                  should not fold               reverse match should fold
    #                                    1         1         1         1
    #      7         8         9         0         1         2         3
    #      0123456789012345678901234567890123456789012345678901234567890123456789
    #                CNNNNQCAGUGYMMMMZ                  ZMMMMYCACUGNNNNNG
    seq +="AAAAAAAAAACCCGUCCAGUGAGCUUUUUUUUUUUGGGGGGGGGGAGCAACCACUGUUGCUGCCCCCCCC"
    #             hferritin should fold              reverse GAGAG should fold
    #      1         1         1         1         1         1         2         2
    #      4         5         6         7         8         9         0         1
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNQCAGUGYMMMMZ                  ZMMMMYCUCUCNNNNNG
    seq +="AAAGUUCUUGCUUCAACAGUGUUUGAACGGAACAAGGGGGGGGGGAGCAACCUCUCUUGCUGCCCCCCCC"
    #             GAGAG pattern, should fold          UAGUA pattern, should fold
    #      2         2         2         2         2         2         2         2
    #      1         2         3         4         5         6         7         8
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNQGAGAGYMMMMZ                    CNNNNNUAGUAYMMMMZ
    seq +="AAAGUUCUUGCUUCAAGAGAGUUUGAACGGAACAAGGGGGGGGGGGGCAACUAUAGUACUGGUUCCCCCC"
    #                  should not fold            reverse UAGUA pattern, should fold
    #      2         2         3         3         3         3         3         3
    #      8         9         0         1         2         3         4         5
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNQCAGUGYMMMMZ                    ZMMMMYUACUANNNNNG
    seq +="AAAUUAAGCGCCCCCCCAGUGGUGAGGUCUUCGGAACAAGGGGGGGGAACCAGUACUAUGGUUGCCCCCC"
    #                FtsZ, should fold          Mb.tub. ideR, folds with one mismatch
    #      3         3         3         3         3         4         4         4
    #      5         6         7         8         9         0         1         2
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNQGAGAGYMMMMZ                    CNNNNNCAGCGYMMMMZ
    seq +="CCACGUAACUCGAGGCGAGAGGCCUUCGACGUGGCAAAAAAUGCCCGCCGCGCCAGCGGCGGGCAUACCG"

    positions = [(10, 21), (45,56), (121, 132), (150, 161), (191, 202), (220, 231), (257,268),
             (333, 344), (360, 371), (397, 408)]

    forward_pattern = re.compile(FORWARD_PATTERN)
    reverse_pattern = re.compile(REVERSE_PATTERN)

    f_matches = find_matches(seq, forward_pattern, FORWARD)
    r_matches = find_matches(seq, reverse_pattern, REVERSE)

    matches = f_matches
    matches.extend(r_matches)

    if len(matches) != 12:
        print "%s: Found %d of 12 possible matches" % (name, len(matches))
        for match in matches:
            print match
        return False

    stem_loops = filter_fold(matches)

    if len(stem_loops) != 10:
        print "%s: Found %d stemloops, 10 expected" % (name, len(stem_loops))
        for match in stem_loops:
            print match
        return False

    for match in stem_loops:
        if (match.start(), match.end()) not in positions:
            print "%s: Unknown match at position (%s:%s)" %(name, match.start(), match.end())
            return False

    return True

def test_annotate_with_cds(name):
    print "%s: implement me!" % name
    return False

def test_filter_UTR(name):
    print "%s: implement me!" % name
    return False

def test_create_pretty_fold_graph(name):
    print "%s: implement me!" % name
    return False

###############################################################################
def testit(test_fun, name):
    global successes
    global failures
    global total

    print "\nTesting %s" % name
    total += 1
    if test_fun(name):
        print "SUCCESS: %s" % name
        successes += 1
    else:
        print "FAILED: %s" % name
        failures += 1

def knownfail(test_fun, name):
    global successes
    global failures
    global total

    print "\nTesting %s" % name
    total += 1
    if test_fun(name):
        print "FAILED: %s unexpected success" % name
        failures += 1
    else:
        print "SUCCES: %s failed as expected" % name
        successes += 1

if __name__ == "__main__":
    testit(test_reverse_complement, "Reverse complement")
    testit(test_match, "Match object member functions")
    testit(test_find_matches, "Match finder")
    knownfail(test_find_cds_features, "CDS filter")
    knownfail(test_is_upstream, "is_upstream")
    knownfail(test_is_downstream, "is_downstream")
    knownfail(test_is_within, "is_within")
    testit(test_can_pair, "Pairing function")
    testit(test_filter_fold, "Folding filter")
    knownfail(test_annotate_with_cds, "Feature annotation")
    knownfail(test_filter_UTR, "UTR region filter")
    knownfail(test_create_pretty_fold_graph, "RNAfold call")

    print "\n%d tests run, %d suceeded, %d failed" % (total, successes, failures)

