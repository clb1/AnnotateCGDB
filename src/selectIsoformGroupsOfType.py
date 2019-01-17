#!/usr/bin/env python

import bz2
import gzip
import sys
from collections import defaultdict

import pdb

# 0) Get IGs for all pORFs that are of type 'type_selector'.
# 1) Each isoform is associated with one or more IG. Get the IGs for each isoform.
# 2) Group the isoforms that have the same the same set of IGs.
# 3) Isoforms with the same set of IGs all encode the same pORF(s). So these isoforms should be primered together and so form the primer groups.


def collectEvidencePerPORF(type_selector, porf_annotations):
    IGs_of_type = set()

    assert (type_selector in ["PM", "GPI"])
    index = 1 if (type_selector == "PM") else 2
    
    ip = gzip.open(porf_annotations, 'rb')
    for line in ip:
        fields = line.strip().split("\t")
        if (fields[index] != '-'):
            IGs_of_type.add( fields[0] )
    ip.close()

    return IGs_of_type


def selectGroups(IGs_of_type, same_porf_isoform_groups, output_tsv):
    ip = bz2.BZ2File(same_porf_isoform_groups, 'rb')

    isoform_to_IGs_map = defaultdict(set)
    
    #GroupID IsoformIDs      GeneSymbol(s)   pORFs   ProteinSequence
    for line in ip:
        fields = line.split("\t")
        if (fields[0] in IGs_of_type):
            pORFs = fields[3].split("|")
            for isoformID in map(lambda x: x.split(',')[0], pORFs):
                isoform_to_IGs_map[isoformID].add(fields[0])
                
    isoforms_per_IG_tuple = defaultdict(list)
    for isoformID, IGs_set in isoform_to_IGs_map.items():
        IGs_tuple = tuple(sorted(IGs_set))
        isoforms_per_IG_tuple[IGs_tuple].append(isoformID)

    ip.close()

    op = bz2.BZ2File(output_tsv, 'wb')
    for IGs_tuple, isoformIDs in isoforms_per_IG_tuple.items():
        op.write("%s\t%s\n" % (" ".join(isoformIDs), str(IGs_tuple)))
    op.close()

    
if (__name__ == "__main__"):
    type_selector, same_porf_isoform_groups, porf_annotations, output_tsv = sys.argv[1:]

    IGs_of_type = collectEvidencePerPORF(type_selector, porf_annotations)
    selectGroups(IGs_of_type, same_porf_isoform_groups, output_tsv)

    sys.exit(0)
    
