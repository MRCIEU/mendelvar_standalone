#!/usr/bin/env python
from sys import argv
from GeneMap import GeneMap
#Script to parse OMIM genemap to produce a clean list of gene-disease associations.
"""Rules:
apping key must be 3. Exclude entries in "[ ]", indicating "nondiseases," mainly genetic variations that lead to apparently abnormal laboratory test values (e.g., dysalbuminemic euthyroidal hyperthyroxinemia).
"""

script, genemap_f = argv

omim = GeneMap(genemap_f)
omim.read_in_genemap()
omim.filter_genemap()
omim.print_genemap()