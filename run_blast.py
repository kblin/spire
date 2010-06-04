#!/usr/bin/env python
#
# SIRE - Search Prokaryote IRE sequences
# Run blastp or blastn against a database
#
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
"""helper to run blastp or blastn without running out of file descirptors"""

import os
import sys
import tempfile
from Bio.Blast import NCBIStandalone
from Bio.Blast import NCBIXML

DB_PATH = "/home/kai/uni/dipl/materials/dbs/"
BLAST_EXE = "/usr/bin/blastall"


def run_blast(fasta_seq, blastdb, program):
    """run the specified blast program"""
    blast_db = DB_PATH + blastdb

    fdesc, file_name = tempfile.mkstemp()
    fhandle = open(file_name, 'w')
    fhandle.write(fasta_seq)
    fhandle.close()

    result_handle, error_handle = NCBIStandalone.blastall(BLAST_EXE, program,
                                                          blast_db, file_name)
    blast_records = NCBIXML.parse(result_handle)

    try:
        # The parser does have a next() function, no matter what pylint thinks
        #pylint: disable-msg=E1101
        blast_record = blast_records.next()
    except ValueError, error:
        print "file not ready yet"
        raise error
    finally:
        error_handle.close()
        result_handle.close()
        os.close(fdesc)
        os.unlink(file_name)

    return blast_record

if __name__ == "__main__":
    #pylint: disable-msg=C0103
    if len(sys.argv) < 4:
        print "Too few arguments"
        sys.exit(2)

    handle = open(sys.argv[1], 'r')
    fasta_lines = handle.readlines()
    handle.close()
    fasta = "".join(fasta_lines)

    print "%s" % run_blast(fasta, sys.argv[2], sys.argv[3])
