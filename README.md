SPIRE
=====

**S**earch for **P**rokaryote **I**ron **R**esponse **E**lements


Running SIPRE
-------------

SPIRE depends on the excellent [BioPython](http://biopython.org) library, and
the [bioinf-helperlibs](https://github.com/kblin/bioinf-helperlibs) helper
library.

If you install the two dependencies, you should be able to run spire from the
directory it is located in.

> ./spire example.gbk

Run without any additional arguments, spire will print out diagnostics on how
many IRE sequences were found. Calling spire with `--help` will show the
additional output options. `--text` will show a human-readable overview of
the results, `--stats` will show some statistics. Multiple output options can
be used together.


License
-------

SPIRE is licensed under the GNU General Public License version 3 or later, see
the `LICENSE` file for details.

> Copyright 2008-2012 Kai Blin


Citing SPIRE
------------

If you found SPIRE useful, please cite
The bifunctional role of aconitase in *Streptomyces viridochromogenes* Tü494
Ewelina Michta, Klaus Schad, Kai Blin, Regina Ort-Winklbauer, Marc Röttig, Oliver Kohlbacher, Wolfgang Wohlleben, Eva Schinko, and Yvonne Mast
*Environmental Microbiology* (2012)
doi: [10.1111/1462-2920.12006](http://dx.doi.org/10.1111/1462-2920.12006)
