#!/usr/bin/env python
#
# SPIRE - Search Prokaryote IRE sequences
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
"""spire - Search Prokaryote IRE sequences
Usage: spire <genome file> <blast database name>

<genome file> is a prokaryote genome file (GenBank format) to scan
<blast database name> is the blast database to use

Output will go to stdout, with some diagnostics printed to stderr.

"""

import re
import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from optparse import OptionParser

FORWARD_PATTERN = '(C|A).....(CAGUG|CAAUG|GAGAG|UAGUA|CAGCG|CUGUG)'
REVERSE_PATTERN = '(CACAG|CACUG|CAUUG|CUCUC|UACUA|CGCUG).....(G|U)'
FORWARD = 1
REVERSE = -1
UTR_LEN = 200
E_VAL_THRESH = 0.04
MAX_ALN_RESULTS = 3

def reverse_complement(sequence):
    """Get the reverse complement of an RNA seqence string"""
    seq = Seq(sequence, IUPAC.unambiguous_rna)
    return seq.reverse_complement().tostring()

def get_best_hsp(alignment):
    """Return the best HSP of an alignment"""
    max_score = 0
    best_hsp = None
    for hsp in alignment.hsps:
        if hsp.score > max_score:
            max_score = hsp.score
            best_hsp = hsp

    return best_hsp

class CDS:
    """coding sequence datatype"""
    def __init__(self, feature, seq_str):
        """create CDS"""
        self.qualifiers = feature.qualifiers
        self.location = feature.location
        self.strand = feature.strand
        self.mrna = seq_str[feature.location.nofuzzy_start-UTR_LEN:\
                            feature.location.nofuzzy_end+UTR_LEN]
        if self.strand == REVERSE:
            self.mrna = reverse_complement(self.mrna)
    def start(self):
        """get the start of the CDS"""
        return self.location.nofuzzy_start
    def end(self):
        """get the end of the CDS"""
        return self.location.nofuzzy_end

class FeatureMatch:
    """The main datatype of SPIRE"""
    def __init__(self, feature, dist, position=""):
        """Create a FeatureMatch"""
        self.feature = feature
        self.dist = dist
        self.position = position
        self.qualifiers = feature.qualifiers
        self.n_align = []
        self.p_align = []
    def start(self):
        """Start of the FeatureMatch"""
        return self.feature.start()
    def end(self):
        """End of the FeatureMatch"""
        return self.feature.end()
    def __str__(self):
        if self.qualifiers.has_key('locus_tag'):
            return self.qualifiers['locus_tag'][0]
        elif self.qualifiers.has_key('gene'):
            return self.qualifiers['gene'][0]
        else:
            return self.qualifiers['note'][0]
    def feature_fasta(self):
        """Print a FASTA version of the feature's RNA sequence"""
        ret = ""
        ret = ">"
        ret += "%s" % self.__str__()
        if self.qualifiers.has_key('protein_id'):
            ret += "|%s" % self.qualifiers['protein_id'][0]
        else:
            ret += "None"
        ret += "|%s:%s|" % (self.start(), self.end())
        if self.qualifiers.has_key('product'):
            ret += "%s %s %s" % (self.dist, self.position,
                                 self.qualifiers['product'][0])
        else:
            ret += "not near any known protein"
        ret += "\n"
        ret += "%s" % self.feature.mrna
        ret += "\n"
        return ret
    def protein_fasta(self):
        """Print a FASTA version of the feature's protein sequence"""
        ret = ""
        ret = ">"
        ret += "%s" % self.__str__()
        if self.qualifiers.has_key('protein_id'):
            ret += "|%s" % self.qualifiers['protein_id'][0]
        else:
            ret += "None"
        ret += "|%s:%s|" % (self.start(), self.end())
        if self.qualifiers.has_key('product'):
            ret += "%s %s %s" % (self.dist, self.position,
                                 self.qualifiers['product'][0])
        else:
            ret += "not near any known protein"
        ret += "\n"
        ret += "%s" % self.qualifiers['translation'][0]
        ret += "\n"
        return ret

def is_same_position(feature_match, alignment):
    """Check if a feature is up- or downstream in the alignments as well"""
    # a simple match on the feature's position doesn't work, as that'll
    # break the full genome blast.
    if feature_match.position == "upstream of":
        position = "downstream"
    elif feature_match.position == "downstream of":
        position = "upstream"
    else:
        #no position? better pretend it's the same position
        return True
    pattern = re.compile(position)
    match = pattern.search(alignment.title)
    if match is None:
        return True
    return False

class Match:
    """datatype for an IRE match"""
    def __init__(self, re_match, sequence, direction):
        """initialize IRE match"""
        self.re_match = re_match
        if direction == REVERSE:
            self.sequence = reverse_complement(sequence)
        else:
            self.sequence = sequence
        self.match_graph = None
        self.direction = direction
        self.features = []
        self.position = 0
    def start(self):
        """start of the match"""
        return self.re_match.span()[0]
    def end(self):
        """end of the match"""
        return self.re_match.span()[1]
    def get_loop(self):
        """get the loop match"""
        return self.sequence[11:16]
    def get_before(self):
        """get stem before the loop match"""
        return self.sequence[6:11]
    def get_after(self):
        """get stem after the loop match"""
        return self.sequence[16:22]
    def __str__(self):
        ret = "Match at (%s:%s)" % (self.start(), self.end())
        ret += "\n\tDirection: "
        if self.direction == FORWARD:
            ret += "forward"
        else:
            ret += "reverse"
        ret += "\n\tSequence: %s" % self.sequence
        if self.match_graph is not None:
            ret += "\n\tFolds to: %s" % self.match_graph
        for feature in self.features:
            ret += "\n\tFeature: %s " % feature.position
            product = "unknown product"
            if feature.qualifiers.has_key('product'):
                product = feature.qualifiers["product"][0]
            protein_id = "no id"
            if feature.qualifiers.has_key('protein_id'):
                protein_id = feature.qualifiers["protein_id"][0]
            ret += "%s (%s)" % (product, protein_id)
            ret += "\n\tDistance: %s" % feature.dist
            ret += "\n\tPosition: %s..%s" % (feature.start(), feature.end())
            if feature.qualifiers.has_key('locus_tag'):
                ret += "\n\tURL: http://www.ncbi.nlm.nih.gov/sites/entrez?"
                ret += "db=gene&cmd=search&term="
                ret += feature.qualifiers["locus_tag"][0]
            ret += "\n\tAlignments: %s protein-based, %s rna-based" % \
                    (len(feature.p_align), len(feature.n_align))
            count = 0
            for align in feature.p_align:
                if not is_same_position(feature, align):
                    continue
                if not count < MAX_ALN_RESULTS:
                    ret += "\n\t\t..."
                    break
                ret += "\n\t\tp-a: %s" % align.title
                length = len(feature.qualifiers['translation'][0])
                hsp = get_best_hsp(align)
                len_pc = (hsp.align_length / float(length)*100)
                id_pc = (hsp.identities / float(length)*100)
                ret += " (%0.2f%%, id: %0.2f%%, gaps: %s)" % (len_pc, id_pc,
                                                              hsp.gaps)
                count += 1
            count = 0
            for align in feature.n_align:
                if not is_same_position(feature, align):
                    continue
                if not count < MAX_ALN_RESULTS:
                    ret += "\n\t\t..."
                    break
                ret += "\n\t\tn-a: %s" % align.title
                length = feature.end() - feature.start()
                hsp = get_best_hsp(align)
                len_pc = (hsp.align_length / float(length)*100)
                id_pc = (hsp.identities / float(length)*100)
                ret += " (%0.2f%%, id: %0.2f%%, gaps: %s)" % (len_pc, id_pc,
                                                              hsp.gaps)
                count += 1
        return ret
    def fasta(self):
        """create FASTA output for match"""
        ret = ">%s:%s" % (self.start(), self.end())
        for feature in self.features:
            if feature.qualifiers.has_key('protein_id'):
                ret += "|%s" % feature.qualifiers['protein_id'][0]
            else:
                ret += "|None"
            if feature.qualifiers.has_key('product'):
                ret += "%s %s" % (self.position,
                                  feature.qualifiers['product'][0])
            else:
                ret += "not near any known protein"
        ret += "\n"
        ret += self.sequence
        return ret
    def feature_fasta(self):
        """create nucleotide FASTA output for all features"""
        ret = ""
        for feature in self.features:
            ret += feature.feature_fasta()
        return ret
    def protein_fasta(self):
        """create protein FASTA output for all features"""
        ret = ""
        for feature in self.features:
            ret += feature.protein_fasta()
        return ret

def find_matches(sequence, search_pattern, direction):
    """find IRE regex matches on a sequence"""
    match_list = []
    itr = search_pattern.finditer(sequence)
    for i in itr:
        offset_d = 5
        offset_u = 10
        if direction == REVERSE:
            offset_d = 10
            offset_u = 5

        match_list.append(Match(i,
                                sequence[i.span()[0]-offset_d:\
                                         i.span()[1]+offset_u],
                                direction))

    return match_list

def find_cds_features(seq_item, sequence):
    """find CDS feature annotations"""
    cds_list = []
    for feature in seq_item.features:
        if feature.type == "CDS":
            cds_list.append(CDS(feature, sequence))
    return cds_list

def set_position(match, feature, direction):
    """create FeatureMatch instance for a match close to a feature"""
    if match.end() < feature.start() and \
       feature.start() - match.end() < UTR_LEN:
        f_match = FeatureMatch(feature, feature.start() - match.end())
        if direction == FORWARD:
            f_match.position = "upstream of"
        else:
            f_match.position = "downstream of"
        match.features.append(f_match)
    elif match.start() > feature.end() and \
         match.start() - feature.end() < UTR_LEN:
        f_match = FeatureMatch(feature, match.start() - feature.end())
        if direction == FORWARD:
            f_match.position = "downstream of"
        else:
            f_match.position = "upstream of"
        match.features.append(f_match)
    elif match.start() > feature.start() and match.end() < feature.end():
        match.features.append(FeatureMatch(feature, 0, "within"))
        return
    else:
        return

def can_pair(base_a, base_b):
    """Check if two RNA bases can pair"""
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

def filter_fold(match_list):
    """Filter regex matches so only stem loops remain"""
    stem_loops = []
    for match in match_list:
        #pylint: disable-msg=W0511
        # Now comes the tricky part, a seq looks like this:
        # XXXXXCNNNNNCAGUGMMMMMXXXXX
        # Now, the NNNNN before the loop should pair to MMMMM after the loop.
        # Unfortunately, U=G pairs can happen as well, so simply
        # checking if the reverse complement of NNNNN is equal to MMMMM
        # does not work. Also, we might have one more base than just
        # CAGUG in the loop, as in hferritin.
        # We need to check base by base, taking the extra base into
        # account
        stem_before = match.get_before()
        stem_after = match.get_after()

        does_pair = False
        mismatches = 0
        # offset to skip non-pairing first base
        offset = 0
        i = 0
        while(i<5):
            i_b = (-1 * i) - 1
            i_a = i + offset

            try:
                does_pair = can_pair(stem_before[i_b], stem_after[i_a])
            except IndexError:
                does_pair = False
                break
            # We allow one mismatch total
            if not does_pair:
                mismatches += 1
            # We allow a mismatch at the first base after the loop
            if offset ==  0 and mismatches > 1:
                offset = 1
                i = 0
                mismatches = 0
                continue

            i += 1
            if mismatches > 1:
                break
        if mismatches < 2:
            stem_loops.append(match)

    return stem_loops

def find_close_features(match, cds_list):
    """get a list of features close to the match"""
    len_cds_list = len(cds_list)
    left = 0
    right = len_cds_list
    while left < right :
        i = (left + right) / 2
        feature = cds_list[i]
        if feature.end() > match.start():
            right = i
        elif match.start() - feature.end() > UTR_LEN:
            left = i + 1
        else:
            break
    start = i > 5 and i - 5 or 0

    return cds_list[start:i+5]

def annotate_with_cds(match_list, cds_list):
    """Annotate features with nearby CDS"""
    for match in match_list:
        close_features = find_close_features(match, cds_list)
        for feature in close_features:
            # do a binary search for close features
            if match.direction != feature.strand:
                continue
            set_position(match, feature, feature.strand)

def filter_utr(match_list):
    """Filter out matches not in the 5' or 3' UTR"""
    filtered_matches = []
    for match in match_list:
        if match.features != []:
            filtered_matches.append(match)
    return filtered_matches

def create_pretty_fold_grap(match_list):
    """Create a fold graph by running RNAfold"""
    for match in match_list:
        pipe = subprocess.Popen("RNAfold -noPS", shell=True,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE)
        try:
            pipe.stdin.write(match.sequence)
            pipe.stdin.close()
            lines = pipe.stdout.readlines()
            match.match_graph = lines[1].rstrip()
        finally:
            pipe.stdin.close()
            pipe.stdout.close()

def run_blastn(match, blastdb):
    """run blastn"""
    from Bio.Blast.Applications import NcbiblastnCommandline
    for feature in match.features:
        print >> sys.stderr, ".",
        rec = None
        try:
            cline = NcbiblastnCommandline(db=blastdb, outfmt=5, num_threads=4)
            pipe = subprocess.Popen(str(cline), shell=True,
                                    stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            pipe.stdin.write(feature.feature_fasta())
            pipe.stdin.close()
            recs = NCBIXML.parse(pipe.stdout)
            rec = recs.next()
            pipe.stdout.close()
            pipe.stderr.close()
        except OSError, err:
            print >> sys.stderr, "Failed to run blastp: %s" % err
            continue
        except ValueError, err:
            print >> sys.stderr, "Parsing blast output failed: %s" % err
            continue
        if not rec:
            continue
        for aln in rec.alignments:
            for hsp in aln.hsps:
                if hsp.expect < E_VAL_THRESH:
                    feature.n_align.append(aln)
                    break

def run_blastp(match, blastdb):
    """run blastp"""
    from Bio.Blast.Applications import NcbiblastpCommandline
    for feature in match.features:
        print >> sys.stderr, ".",
        rec = None
        fasta = feature.protein_fasta()
        if fasta == "":
            continue
        try:
            cline = NcbiblastpCommandline(db=blastdb, outfmt=5, num_threads=4)
            pipe = subprocess.Popen(str(cline), shell=True,
                                    stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            pipe.stdin.write(fasta)
            pipe.stdin.close()
            recs = NCBIXML.parse(pipe.stdout)
            rec = recs.next()
            pipe.stdout.close()
            pipe.stderr.close()
        except OSError, err:
            print >> sys.stderr, "Failed to run blastp: %s" % err
            continue
        except ValueError, err:
            print >> sys.stderr, "Parsing blast output failed: %s" % err
            continue
        if not rec:
            continue
        for aln in rec.alignments:
            for hsp in aln.hsps:
                if hsp.expect < E_VAL_THRESH:
                    feature.p_align.append(aln)
                    break


def blast(match_list, blastdb):
    """run blast searches"""
    for match in match_list:
        run_blastn(match, blastdb)
        run_blastp(match, blastdb)

def filter_align(match_list):
    """filter out matches that don't have an alignment to other genomes"""
    matches = []
    for match in match_list:
        appended = False
        for feature in match.features:
            if feature.n_align != [] or feature.p_align != []:
                if not appended:
                    matches.append(match)
                    appended = True
            else:
                match.features.remove(feature)

    return matches

def filter_real_protein(match_list):
    """filter out hypothetical proteins"""
    pat = re.compile("hypothetical protein")
    matches = []
    for match in match_list:
        remove_list = []
        for feat in match.features:
            if not feat.qualifiers.has_key('product'):
                continue
            product = feat.qualifiers['product'][0]
            mat = pat.search(product)
            if mat is not None:
                remove_list.append(feat)
        for feat in remove_list:
            match.features.remove(feat)
        if match.features != []:
            matches.append(match)
    return matches


def calculate_stats(matches):
    """Calculate stats for the different loop sequences"""
    stats = {}
    for hit in matches:
        loop = hit.get_loop()
        if stats.has_key(loop):
            stats[loop] += 1
        else:
            stats[loop] = 1

    print >> sys.stderr, "Statistics:\nLoop seq\tcount\n--------\t-----"
    for key in sorted(stats.keys()):
        print "%s\t\t%s" % (key, stats[key])


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

    opts, args = parser.parse_args()
    if len(args) < 1:
        parser.error("Please specify a GenBank file")


    handle = open(args[0])
    seq_iterator = SeqIO.parse(handle, "genbank")

    forward_pattern = re.compile(FORWARD_PATTERN)
    reverse_pattern = re.compile(REVERSE_PATTERN)

    seq_str = ""
    seq_i = seq_iterator.next()

    # Find our pattern in the genome
    seq_str = seq_i.seq.transcribe().tostring()
    forward_matches = find_matches(seq_str, forward_pattern, FORWARD)
    reverse_matches = find_matches(seq_str, reverse_pattern, REVERSE)

    cds_list = find_cds_features(seq_i, seq_str)

    matches = forward_matches
    matches.extend(reverse_matches)

    print "Found %s forward, %s reverse matches" % (len(forward_matches),
                                                    len(reverse_matches))

    matches = filter_fold(matches)
    print >> sys.stderr, "Found %s stem loops" % len(matches)
    annotate_with_cds(matches, cds_list)
    print "Found %s stem loops" % len(matches)
    print >> sys.stderr, "Running UTR filter"
    matches = filter_utr(matches)
    print "Found %s matches in UTRs" % len(matches)
    print >> sys.stderr, "%s UTR matches" % len(matches)
    if opts.blastdb != "none":
        blast(matches, opts.blastdb)
        matches = filter_align(matches)
        print "Found %s matches with alignments with %s" % \
              (len(matches), blastdb)
    if opts.filter_hypothetical:
        matches = filter_real_protein(matches)
        print "Found %s matches after filtering hypothetical proteins" % \
              len(matches)
    if opts.fold_rna:
        print >> sys.stderr, "Creating RNAFold graphs"
        create_pretty_fold_grap(matches)

    if opts.nucleotide:
        for hit in matches:
            print "%s" % hit
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

        out_handle = open("spire_%s" % genome_file, 'w')
        SeqIO.write([seq_i], out_handle, "genbank")
        out_handle.close()
    if opts.stats:
        calculate_stats(matches)

    handle.close()

if __name__ == "__main__":
    main(sys.argv)
