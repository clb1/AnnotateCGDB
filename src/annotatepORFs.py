#!/usr/bin/env python

from collections import defaultdict
import gzip
import os
import pandas as pd
import sys

import pdb
    


def collectDataForMatchesToPMProteins(annotated_PM_proteins, blast_directory):
    pORF_to_PM_protein_matches = defaultdict(set)

    blast_file_columns = ["query", "subject", "perc_ident", "align_len", "num_mismatches", "num_gaps", \
                          "query_start", "query_end", "subject_start", "subject_end", "E_value", "bit_score"]

    ip = gzip.open(annotated_PM_proteins, 'rb')
    header = ip.readline()
    for line in ip:
        uniprot_id, aa_seq, tm_segments = line.strip().split("\t")
        print >> sys.stderr, "\rINFO: processing Blast results for %s" % uniprot_id,
        tm_segments_positions = []
        for tm_segment in tm_segments.split(','):
            start, stop = tm_segment.split('-')
            tm_segments_positions.append( set(range(int(start),int(stop)+1)) )

        blast_file = "%s/%s.blastp.bz2" % (blast_directory, uniprot_id)
        assert(os.path.exists(blast_file))

        min_TM_match_len = 15
        blast_df = pd.read_csv(blast_file, sep="\t", header=None, skiprows=5, names=blast_file_columns, compression="bz2", skipfooter=1)
        for tup in blast_df.itertuples(index=False):
            match_positions = set(range(int(tup[6]), int(tup[7]+1)))
            if (any(map(lambda x: len(match_positions & x) >= min_TM_match_len, tm_segments_positions))):
                pORF_to_PM_protein_matches[tup[1]].add(uniprot_id)

    ip.close()
    print >> sys.stderr, "INFO: done\n"
    return pORF_to_PM_protein_matches


def writeOutput(pORF_to_PM_protein_matches, output_tsv):
    print >> sys.stderr, "INFO: writing output to %s" % output_tsv
    op = gzip.open(output_tsv, 'wb')
    op.write("pORF_ID\tPM_PROTEIN_ID\tGPI_PROTEIN_ID\n")
    for pORF_ID, PM_protein_IDs in pORF_to_PM_protein_matches.items():
        for protein_ID in PM_protein_IDs:
            op.write("%s\t%s\t-\n" % (pORF_ID, protein_ID))
    op.close()


if (__name__ == "__main__"):
    annotated_PM_proteins, blast_directory, output_tsv = sys.argv[1:]

    # Record which regions of each HMPAS protein aligned to which pORFs
    pORF_to_PM_protein_matches = collectDataForMatchesToPMProteins(annotated_PM_proteins, blast_directory)

    writeOutput(pORF_to_PM_protein_matches, output_tsv)

    sys.exit(0)
