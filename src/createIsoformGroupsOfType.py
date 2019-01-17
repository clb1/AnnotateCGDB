#!/usr/bin/env python

import bz2
import gzip
import sys
from collections import defaultdict

# 0) Get IGs for all pORFs that are of type 'type_selector'.
# 1) Each isoform is associated with one or more IG. Get the IGs for each isoform.
# 2) Group the isoforms that have the same the same set of IGs.
# 3) Isoforms with the same set of IGs all encode the same pORF(s). So these isoforms should be primered together and so form the primer groups.


def collectEvidencePerPORF(type_selector, porf_annotations):
    isoform_to_uniprot = defaultdict(set)

    assert (type_selector in ["PM", "GPI", "PM_and_GPI"])

    if (type_selector in ["PM", "PM_and_GPI"]):
        ip = gzip.open(porf_annotations, 'rb')
        for line in ip:
            fields = line.strip().split("\t")
            if (fields[1] != '-'):
                isoform_to_uniprot[fields[0]].add( (fields[1],fields[3]) )
        ip.close()

    if (type_selector in ["GPI", "PM_and_GPI"]):
        ip = gzip.open(porf_annotations, 'rb')
        for line in ip:
            fields = line.strip().split("\t")
            if (fields[2] != '-'):
                isoform_to_uniprot[fields[0]].add( (fields[2],fields[3]) )
        ip.close()

    return isoform_to_uniprot


def selectGroups(isoform_to_uniprot, same_porf_isoform_groups, output_tsv):
    ip = bz2.BZ2File(same_porf_isoform_groups, 'rb')
    line = ip.readline() # header
    
    isoforms_of_type = set(isoform_to_uniprot.keys())
    isoform_to_IGs_map = defaultdict(set)
    gene_symbols_per_IG = {}
    
    #GroupID IsoformIDs      GeneSymbol(s)   pORFs   ProteinSequence
    for line in ip:
        fields = line.split("\t")
        isoformIDs = set(fields[1].split('|'))
        if (isoformIDs.issubset(isoforms_of_type)):
            gene_symbols_per_IG[fields[0]] = set(filter(lambda x: x!='None', fields[2].split(',')))
            for isoformID in isoformIDs:
                isoform_to_IGs_map[isoformID].add(fields[0])
    ip.close()


    isoforms_per_IG_tuple = defaultdict(list)
    for isoformID, IGs_set in isoform_to_IGs_map.items():
        IGs_tuple = tuple(sorted(IGs_set))
        isoforms_per_IG_tuple[IGs_tuple].append(isoformID)

    op = bz2.BZ2File(output_tsv, 'wb')
    op.write("IsoformIDs\tUniProtIDsAssoc\tGeneSymbols\tIsoformGroups\n")
    for IGs_tuple, isoformIDs in isoforms_per_IG_tuple.items():
        uniprot_assoc_tuples = set.union(*map(lambda x: isoform_to_uniprot[x], isoformIDs))
        uniprot_assoc_tuples_str = map(lambda x: "%s,%s" % x, uniprot_assoc_tuples)
        associated_uniprot_str = " ".join(uniprot_assoc_tuples_str) if (len(uniprot_assoc_tuples_str)>0) else "None"
        associated_gene_symbols = sorted(set.union(*map(lambda y: gene_symbols_per_IG[y], IGs_tuple)))
        associated_gene_symbols_str = ",".join(associated_gene_symbols)
        IGs_tuple_str = ",".join(IGs_tuple)
        op.write("%s\t%s\t%s\t%s\n" % (" ".join(isoformIDs), associated_uniprot_str, associated_gene_symbols_str, IGs_tuple_str))
    op.close()

    
if (__name__ == "__main__"):
    type_selector, same_porf_isoform_groups, porf_annotations, output_tsv = sys.argv[1:]

    isoform_to_uniprot = collectEvidencePerPORF(type_selector, porf_annotations)
    selectGroups(isoform_to_uniprot, same_porf_isoform_groups, output_tsv)

    sys.exit(0)
    
