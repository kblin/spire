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
from nose.tools import *

sys.path.append('.')
sys.path.append('..')
from search_ire import *

class FakeReMatch:
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def span(self):
        return (self.start, self.end)

def test_reverse_complement():
    seq =   "AGGCU"
    compl = "AGCCU"
    assert_equal(reverse_complement(seq), compl)
    assert_equal(reverse_complement(compl), seq)

def test_match():
    start = 11
    end = 21
    seq = "XXXXXCNNNNNCAGUGYMMMMZXXXXX"
    before = "NNNNN"
    after =  "YMMMMZ"
    loop =   "CAGUG"
    direction = FORWARD
    m = Match(FakeReMatch(start, end), seq, direction)

    assert_equal(m.start(), start)
    assert_equal(m.end(), end)
    assert_equal(m.get_before(), before)
    assert_equal(m.get_loop(), loop)
    assert_equal(m.get_after(), after)

def test_find_matches():
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

    ok_(matches, "failed to find matches")
    assert_equal(len(matches), 8)

    m1 = matches[0]
    m2 = matches[1]
    m5 = matches[3] #forward matches are found first
    m4 = matches[6]

    assert_equal(m1.start(), 10)
    assert_equal(m1.end(), 21)
    assert_equal(m1.sequence, "AAAAACAAGUUCAGUGAGCUUUUUUU")

    assert_equal(m2.start(), 45)
    assert_equal(m2.end(), 56)
    assert_equal(m2.sequence, "GGGGGCAGGCUCAGUGCAGCUUCCCC")

    assert_equal(m4.direction, REVERSE)
    assert_equal(m4.start(), 121)
    assert_equal(m4.end(), 132)
    assert_equal(m4.sequence, reverse_complement("GGGGAGCAACCACUGUUGCUGCCCCC"))
    assert_equal(len(m4.sequence), len("GGGGAGCGACCACUGUUGCUGCCCCC"))

def test_find_cds_features():
    # FIXME: Implement this test
    pass


def test_is_same_position():
    # FIXME: Implement this test
    pass

def test_can_pair():
    test_seq = "AcGGu"
    correct_seq = "UGcUA"
    wrong_seq = "CAaAC"

    for i in range(0,5):
        t = test_seq[i]
        c = correct_seq[i]
        w = wrong_seq[i]

        assert_true(can_pair(t, c))
        assert_false(can_pair(t, w))

def test_filter_fold():
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

    assert_equal(len(matches), 12)

    stem_loops = filter_fold(matches)

    assert_equal(len(stem_loops), 10)

    for match in stem_loops:
        ok_((match.start(), match.end()) in positions, "Unknown match found")

def test_annotate_with_cds():
    # FIXME: Implement this test
    pass

def test_filter_UTR():
    # FIXME: Implement this test
    pass

def test_create_pretty_fold_graph():
    forward_pattern = re.compile(FORWARD_PATTERN)
    seq = "CCACGUAACUCGAGGCGAGAGGCCUUCGACGUGGCA"
    exp =           ".(((((.....)))))"
    matches = find_matches(seq, forward_pattern, FORWARD)
    assert_equal(len(matches), 1)
    assert_equal(matches[0].fold_graph, None)
    matches = filter_fold(matches)
    ok_(matches[0].fold_graph != None)
    assert_equal(matches[0].fold_graph, exp)
    pretty_fold_expected = """
     G
  A     A
 G       G
   C - G
   G - C
   G - C
   A - U
   G - U
 C
   U   C
"""
    matches = create_pretty_fold_grap(matches)
    assert_equal(matches[0].two_d_fold_graph, pretty_fold_expected)
