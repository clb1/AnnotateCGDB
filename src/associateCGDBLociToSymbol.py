#!/usr/bin/env python

from collections import defaultdict
import cPickle
import bz2
import sys

import pdb


if (__name__ == "__main__"):
    symbol_to_uniprot_tsv, isoform_groups_bz2 = sys.argv[1:]

    uniprot_to_symbol = defaultdict(set)
    cgdb_to_symbol = defaultdict(set)
    all_cgdb_loci = set()
    
    with open(symbol_to_uniprot_tsv, 'r') as ip:
        for line in ip:
            fields = line.split("\t")
            symbol = fields[0]
            try:
                uniprot = fields[5].strip()
            except IndexError:
                print line
                pdb.set_trace()
                
            if (uniprot != '' and fields[0] != "Approved Symbol"):
                uniprot_to_symbol[uniprot].add(symbol)

    with bz2.BZ2File(isoform_groups_bz2,'r') as ip:
        for line in ip:
            fields = line.split("\t")
            CGDB_locus = fields[0].split('.')[0]
            all_cgdb_loci.add(CGDB_locus)
            for uniprot in fields[1].split(','):
                if (uniprot_to_symbol.has_key(uniprot)):
                    for symbol in uniprot_to_symbol[uniprot]:
                        cgdb_to_symbol[CGDB_locus].add(symbol)

    cPickle.dump(cgdb_to_symbol,open("cgdb_to_symbol.pkl",'wb'))

    # Which CGDB loci are associated with multiple symbols?
    op_debug = open("multi_symbols.txt", 'w')
    for CGDB_locus, symbols in cgdb_to_symbol.items():
        if (len(symbols)>1):
            op_debug.write("%s\t%s\n" % (CGDB_locus, " ".join(symbols)))
    op_debug.close()

    # Which CGDB loci are associated with no symbols?
    op_debug = open("no_symbols.txt", 'w')
    for CGDB_locus in all_cgdb_loci:
        if (not cgdb_to_symbol.has_key(CGDB_locus)):
            op_debug.write("%s\n" % CGDB_locus)
    op_debug.close()
    
    sys.exit(1)
