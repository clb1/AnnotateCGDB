#!/usr/bin/env python

from collections import defaultdict
import pandas as pd
import sys


def agglomByIsoformGroups(seq_and_porf_stream):
    pORFs_w_same_seq = defaultdict(set)
    ip = open(seq_and_porf_stream, 'r')
    for line in ip:
        seq, pORF = line.strip().split("\t")
        pORFs_w_same_seq[seq].add(pORF)
    ip.close()
    print >> sys.stderr, "INFO: read all pORF sequences"
    
    unfilt_pORFs_and_seqs_by_isoform_groups = defaultdict(list)
    for seq, pORFs in pORFs_w_same_seq.items():
        isoformIDs = set(map(lambda x:x.split(',')[0], pORFs))
        sorted_isoformIDs = sorted(isoformIDs)
        unfilt_pORFs_and_seqs_by_isoform_groups[tuple(sorted_isoformIDs)].append( (pORFs, seq) )
    print >> sys.stderr, "INFO: grouped pORFs and sequences by isoform IDs"
    
    # Eliminate "contained" (pORFs,seq) tuples
    pORFs_and_seqs_by_isoform_groups = {}
    for isoformIDs_tuple, unfilt_pORFs_and_seqs_list in unfilt_pORFs_and_seqs_by_isoform_groups.items():
        # Remove (non-self) pORF groups whose common sequence is a subsequence of another sequence
        pORFs_and_seqs_list = filter(lambda x: not any(map(lambda y:len(x[1]) < len(y[1]) and x[1] in y[1], unfilt_pORFs_and_seqs_list)), unfilt_pORFs_and_seqs_list)
        pORFs_and_seqs_by_isoform_groups[isoformIDs_tuple] = pORFs_and_seqs_list
    print >> sys.stderr, "INFO: removed contained (pORFs,seq) tuples"

    return pORFs_and_seqs_by_isoform_groups

    
def nameAndWriteIsoformGroups(isoform_gene_symbols, pORFs_and_seqs_by_isoform_groups, output_tsv):
    curr_group_num = 1
    records_for_table = []
    for isoformIDs_tuple, pORFs_and_seqs_list in pORFs_and_seqs_by_isoform_groups.items():
        isoformIDs_tuple_str = "|".join(isoformIDs_tuple)
        gene_symbols = list(set(map(lambda x:isoform_gene_symbols[x], isoformIDs_tuple)))
        gene_symbols.sort()
        gene_symbols_str = ",".join(gene_symbols)
        for pORFs, seq in pORFs_and_seqs_list:
            groupID = "IG%d" % curr_group_num
            curr_group_num += 1
            records_for_table.append( (groupID, isoformIDs_tuple_str, gene_symbols_str, "|".join(pORFs), seq) )

    df = pd.DataFrame.from_records(records_for_table, columns=["GroupID", "IsoformIDs", "GeneSymbol(s)", "pORFs", "ProteinSequence"])
    df.to_csv(output_tsv, sep="\t", index=False)

    
def mapIsoformsToGeneSymbols(CG_named_loci):
    named_loci_df = pd.read_csv(CG_named_loci, sep="\t")
    return named_loci_df.set_index("IsoformID")["GeneSymbol"].to_dict()


if (__name__ == "__main__"):
    seq_and_porf_stream, CG_named_loci, output_tsv = sys.argv[1:]

    isoform_gene_symbols = mapIsoformsToGeneSymbols(CG_named_loci)
    pORFs_and_seqs_by_isoform_groups = agglomByIsoformGroups(seq_and_porf_stream)
    nameAndWriteIsoformGroups(isoform_gene_symbols, pORFs_and_seqs_by_isoform_groups, output_tsv)
    
    sys.exit(0)
    
