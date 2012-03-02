#!/usr/bin/env python
#
# SPIRE - Search Prokaryote IRE sequences
# Copyright (C) 2008-2012 Kai Blin <kai.blin@biotech.uni-tuebingen.de>
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
"""spire - Search Prokaryote IRE sequences
Usage: spire <genome file> <blast database name>

<genome file> is a prokaryote genome file (GenBank format) to scan
<blast database name> is the blast database to use

Output will go to stdout, with some diagnostics printed to stderr.

"""

import sys
import logging
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from optparse import OptionParser
from helperlibs.bio.ire import (
        find_cds_features,
        annotate_with_cds,
        filter_fold,
        blast,
        filter_align,
        filter_utr,
        filter_real_protein,
        filter_for_self_matches,
        calculate_stats,
        process_sequence,
    )


def filter_matches(matches, cds_list, opts):
    """Run all the requested filters"""
    matches = filter_fold(matches)
    logging.info("Found %s stem loops" % len(matches))
    annotate_with_cds(matches, cds_list, opts.hits_within_genes)
    print "Found %s stem loops" % len(matches)
    logging.info("Running UTR filter")
    matches = filter_utr(matches)
    print "Found %s matches in UTRs" % len(matches)
    logging.info("%s UTR matches" % len(matches))
    if opts.blastdb != "none":
        blast(matches, opts.blastdb)
        matches = filter_align(matches)
        print "Found %s matches with alignments with %s" % \
              (len(matches), opts.blastdb)
    if opts.filter_hypothetical:
        matches = filter_real_protein(matches)
        print "Found %s matches after filtering hypothetical proteins" % \
              len(matches)

    if opts.self_hits:
        matches = filter_for_self_matches(matches)
        print "Found %s matches after filtering for self-matches" % \
              len(matches)

    return matches


def main(argv):
    """search prokaryote IRE sequences on a genome"""
    usage = "Usage: %prog [options] <genebank file>"
    parser = OptionParser(usage=usage, version="%prog 0.2.0")
    parser.add_option("-n", "--nucleotide", dest="nucleotide",
                      action="store_true", default=False,
                      help="print nucleotide-based match information")
    parser.add_option("-p", "--protein", dest="protein",
                      action="store_true", default=False,
                      help="print protein fasta sequence of match")
    parser.add_option("-g", "--genbank", dest="genbank",
                      action="store_true", default=False,
                      help="write results into an embl file")
    parser.add_option("-d", "--database", dest="blastdb",
                      metavar="DB", default="none",
                      help="blast database to use, use 'none' to skip blast")
    parser.add_option("-F", "--fold-rna", dest="fold_rna",
                      action="store_true", default=False,
                      help="create RNA fold graphs using RNAFold")
    parser.add_option("-s", "--stats", dest="stats",
                      action="store_true", default=False,
                      help="print statistics on the found IRE sequences")
    parser.add_option("-H", "--filter-hypothetical", dest="filter_hypothetical",
                      action="store_true", default=False,
                      help="Ignore IRE hits in hypothetical proteins")
    parser.add_option("-t", "--text", dest="text",
                      action="store_true", default=False,
                      help="Print textual information")
    parser.add_option("-w", "--wighin-genes", dest="hits_within_genes",
                      action="store_true", default=False,
                      help="allow hits within genes, not only in UTRs")
    parser.add_option("-S", "--self-hits", dest="self_hits",
                      action="store_true", default=False,
                      help="Display only IRE hits that have similar IREs "\
                           "on the same genome")
    parser.add_option("-v", "--verbose", dest="verbose",
                      action="store_true", default=False,
                      help="Print diagnostic messages while running")

    opts, args = parser.parse_args()
    if len(args) < 1:
        parser.error("Please specify a GenBank file")


    # Set up logging
    if opts.verbose:
        logging.basicConfig(format='%(levelname)s:%(message)s',
                            level=logging.INFO)

    # Load sequence
    handle = open(args[0])
    seq_iterator = SeqIO.parse(handle, "genbank")
    seq_i = seq_iterator.next()
    seq = seq_i.seq.transcribe()

    matches  = process_sequence(seq)

    cds_list = find_cds_features(seq_i, seq)

    matches = filter_matches(matches, cds_list, opts)

    if opts.nucleotide:
        for hit in matches:
            print "%s" % hit.feature_fasta()
    if opts.protein:
        for hit in matches:
            print "%s" % hit.protein_fasta()
    if opts.genbank:
        for hit in matches:
            note = "IRE motif, loop sequence %s" % hit.sequence[11:16]
            feature = SeqFeature(FeatureLocation(hit.start(), hit.end()),
                                 strand=hit.direction, type="misc_feature",
                                 qualifiers={'note':note})
            seq_i.features.append(feature)

        out_handle = open("spire_%s" % args[0], 'w')
        SeqIO.write([seq_i], out_handle, "genbank")
        out_handle.close()
    if opts.stats:
        calculate_stats(matches)
    if opts.text:
        for hit in matches:
            print hit
            if opts.fold_rna:
                print hit.two_d_fold_graph

    handle.close()

if __name__ == "__main__":
    main(sys.argv)
