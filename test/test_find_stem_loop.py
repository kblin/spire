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

from find_stem_loop import *

successes = 0
failures = 0
total = 0

def test_can_pair(name):
    if not can_pair("A", "U"):
        return False
    if not can_pair("U", "A"):
        return False
    if not can_pair("C", "G"):
        return False
    if not can_pair("G", "C"):
        return False
    if not can_pair("G", "U"):
        return False
    if not can_pair("U", "G"):
        return False
    if can_pair("A", "C"):
        return False
    if can_pair("C", "A"):
        return False
    if can_pair("C", "U"):
        return False
    if can_pair("U", "C"):
        return False
    return True

def test_find_stem_loop(name):
    seq = "ATGC"
    got_exception = False
    try:
        results = find_stem_loop(seq, 5, 5, 1, 0)
    except SeqTooShortException:
        got_exception = True

    if not got_exception:
        print "%s: Failed to trigger exception on too short sequence" % name
        return False

    seq = "CCCAACUUCCC"
    results = find_stem_loop(seq, 1, 1, 0, 0)

    if len(results) != 1:
        print "%s: looking for stem loop ACU failed" % name
        return False

    res = results[0]
    if res != (4,7, 'ACU'):
        print "%s: stem loop ACU at wrong position. Expected (4,7, 'ACU'), got %s" % (name, res)
        return False

    results = find_stem_loop(seq, 1, 3, 0, 0)

    if len(results) != 1:
        print "%s: looking for stem loop AACUU failed" % name
        return False

    res = results[0]
    if res != (3,8, 'AACUU'):
        print "%s: stem loop AACUU at wrong position. Expected (3, 8, 'AACUU'), got %s" % (name, res)
        return False

    results = find_stem_loop(seq, 2, 1, 0, 0)

    if len(results) != 1:
        print "%s: looking for stem loop AACUU failed" % name
        return False

    res = results[0]
    if res != (3, 8, 'AACUU'):
        print "%s: stem loop AACUU at wrong position. Expected (3, 8, 'AACUU'), got %s" % (name, res)
        return False

    seq = "AAACCAAGGAAA"
    results = find_stem_loop(seq, 2, 1, 1, 0)

    if len(results) != 1:
        print "%s: looking for stem loop CCAUGG failed, found %s loops" % (name, len(results))
        for res in results:
            print res
        return False

    res = results[0]
    if res != (3, 9, 'CCAAGG'):
        print "%s: stem loop 'CCAAGG' at wrong position. Expected (3, 9, 'CCAAGG'), got %s" % (name, res)
        return False

    seq = "GUUCUUGCUUCAACAGUGUUUGAACGGAAC"
    results = find_stem_loop(seq, 5, 5, 1, 0)

    if len(results) != 1:
        print "%s: looking for stem loop 'UUCAACAGUGUUUGAA' failed: found %s loops" % (name, len(results))
        return False

    res = results[0]
    if res != (8, 24, 'UUCAACAGUGUUUGAA'):
        print "%s: stem loop 'UUCAACAGUGUUUGAA' at wrong position. Expected (8, 24, 'UUCAACAGUGUUUGAA'), got %s" % (name, res)
        return False

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

    results = find_stem_loop(seq, 5, 5, 1, 0)

    if len(results) != 6:
        print "%s: Found %s stem loops with 0 mismatches, expected 6" % (name, len(results))
        for r in results:
            print r
        return False

    for r in results:
        if not r in positions:
            print "%s: found unknown match %s" % (name, r)
            return False

    results = find_stem_loop(seq, 5, 5, 1, 1)

    if len(results) != 8:
        print "%s: Found %s stem loops with 1 mismatch, expected 8" % (name, len(results))
        for r in results:
            print r
        return False

    for r in results:
        if not r in positions:
            print "%s: found unknown match %s" % (name, r)
            return False

    return True

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
    testit(test_can_pair, "RNA pairing")
    testit(test_find_stem_loop, "Stem loop detection")
