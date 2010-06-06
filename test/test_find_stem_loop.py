#/usr/bin/env python
#
# SIRE - Search Prokaryote IRE sequences
# Test sliding stem loop detection
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
from find_stem_loop import *

def test_can_pair():
    assert_true(can_pair("A", "U"), "Failed to pair A/U")
    assert_true(can_pair("U", "A"), "Failed to pair U/A")
    assert_true(can_pair("C", "G"), "Failed to pair C/G")
    assert_true(can_pair("G", "C"), "Failed to pair G/C")
    assert_true(can_pair("G", "U"), "Failed to pair G/U")
    assert_true(can_pair("U", "G"), "Failed to pair U/G")
    assert_false(can_pair("A", "A"), "Erroneously paired A/A")
    assert_false(can_pair("C", "C"), "Erroneously paired C/C")
    assert_false(can_pair("G", "G"), "Erroneously paired G/G")
    assert_false(can_pair("U", "U"), "Erroneously paired U/U")
    assert_false(can_pair("A", "C"), "Erroneously paired A/C")
    assert_false(can_pair("C", "A"), "Erroneously paired C/A")
    assert_false(can_pair("C", "U"), "Erroneously paired C/U")
    assert_false(can_pair("U", "C"), "Erroneously paired U/C")

def test_find_stem_loop():
    ###
    ### Simple tests
    ###

    seq = "ATGC"
    assert_raises(SeqTooShortException, find_stem_loop, seq, 5, 5, 1, 0)

    seq = "CCCAACUUCCC"

    # find 3-base "(.)" loop ACU
    results = find_stem_loop(seq, 1, 1, 0, 0)
    ok_(len(results) == 1, "looking for stem loop ACU failed")

    res = results[0]
    assert_equal(res, (4,7, 'ACU'), "stem loop ACU at wrong position.")

    # find 5-base "(...)" loop AACUU
    results = find_stem_loop(seq, 1, 3, 0, 0)
    ok_(len(results) == 1, "looking for stem loop AACUU failed")

    res = results[0]
    assert_equal(res, (3,8, 'AACUU'), "stem loop AACUU at wrong position.")

    # find 5-base "((.))" loop AACUU
    results = find_stem_loop(seq, 2, 1, 0, 0)
    ok_(len(results) == 1, "looking for stem loop AACUU failed")

    res = results[0]
    assert_equal(res, (3, 8, 'AACUU'), "stem loop AACUU at wrong position")

    # find 6-base "((..))" loop CCAUGG, where U is a loop fuzz
    seq = "AAACCAAGGAAA"
    results = find_stem_loop(seq, 2, 1, 1, 0)

    ok_(len(results) == 1, "looking for stem loop CCAAGG failed")

    res = results[0]
    assert_equal(res, (3, 9, 'CCAAGG'),
                 "stem loop 'CCAAGG' at wrong position.")

    # find 16-base "(((((......)))))" loop UUCAACAGUGUUUGAA
    seq = "GUUCUUGCUUCAACAGUGUUUGAACGGAAC"
    results = find_stem_loop(seq, 5, 5, 1, 0)

    ok_(len(results) == 1, "looking for stem loop 'UUCAACAGUGUUUGAA' failed")

    res = results[0]
    assert_equal(res, (8, 24, 'UUCAACAGUGUUUGAA'),
                 "stem loop 'UUCAACAGUGUUUGAA' at wrong position.")

    ###
    ### More complex test cases
    ###

    #             Match without extra base            Match with extra base
    #                1         2         3         4         5         6         7
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNNCAGUGYMMMMZ                  CNNNNNCAGUGYMMMMZ
    seq = "XXXXXXXXXXCAAGUUCAGUGAGCUUUXXXXXXXXXXXXXXXXXXCAGGCUCAGUGCAGCUUXXXXXXXX"
    #                  should not fold               no reverse matches....
    #                                    1         1         1         1         1
    #      7         8         9         0         1         2         3         4
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNNCAGUGYMMMMZ
    seq +="XXXXXXXXXXCCCGUCCAGUGAGCUUUXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    #             hferritin should fold              no reverse matches....
    #      1         1         1         1         1         1         2         2
    #      4         5         6         7         8         9         0         1
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNNCAGUGYMMMMZ
    seq +="XXXXXXXXXXCUUCAACAGUGUUUGAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    #             GAGAG pattern, should fold          UAGUA pattern, should fold
    #      2         2         2         2         2         2         2         2
    #      1         2         3         4         5         6         7         8
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNNGAGAGYMMMMZ                    CNNNNNUAGUAYMMMMZ
    seq +="XXXXXXXXXXCUUCAGGAGAGCCUGAAXXXXXXXXXXXXXXXXXXXXCAACUAUAGUACUGGUUXXXXXX"
    #                  should not fold             no reverse matches...
    #      2         2         3         3         3         3         3         3
    #      8         9         0         1         2         3         4         5
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNNCAGUGYMMMMZ
    seq +="XXXXXXXXXXCCCCCCCAGUGGUGAGGXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    #                FtsZ, should fold          Mb.tub. ideR, folds with one mismatch
    #      3         3         3         3         3         4         4         4
    #      5         6         7         8         9         0         1         2
    #      01234567890123456789012345678901234567890123456789012345678901234567890
    #                CNNNNNGAGAGYMMMMZ                    CNNNNNCAGCGYMMMMZ
    seq +="XXXXXXXXXXCGAGGCGAGAGGCCUUCXXXXXXXXXXXXXXXXXXXXCCGCGCCAGCGGCGGGCXXXXXX"

    positions = [(11,  27,  'AAGUUCAGUGAGCUUU'),
                 (46,  62,  'AGGCUCAGUGCAGCUU'),
                 (151, 167, 'UUCAACAGUGUUUGAA'),
                 (221, 237, 'UUCAGGAGAGCCUGAA'),
                 (258, 274, 'AACUAUAGUACUGGUU'),
                 (361, 377, 'GAGGCGAGAGGCCUUC'),
                 (397, 413, 'CCGCGCCAGCGGCGGG'),
                 (398, 414, 'CGCGCCAGCGGCGGGC')]

    # find 16-base loops without mismatches
    results = find_stem_loop(seq, 5, 5, 1, 0)

    ok_(len(results) == 6, "Did not find 6 stem loops")

    for r in results:
        ok_(r in positions, "found unknown match")


    # find 16-base loops with one or less mismatches
    results = find_stem_loop(seq, 5, 5, 1, 1)

    ok_(len(results) == 8, "Did not find 8 stem loops with < 2 mismatches")

    for r in results:
        ok_(r in positions, "found unknown match")
