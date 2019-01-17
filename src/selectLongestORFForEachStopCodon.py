#!/usr/bin/env python

import re
import sys
from Bio import SeqIO

def findLongestReprProteins(ip):
    longest_per_isoform_and_stop = {}
    reID = re.compile("(CGDB\d+.\d+.\d+.\d+),(\S\S\S)(\d+)-(\S\S\S)(\d+)")
    
    #ID format: CG23111.2.1,aug15-uga240
    for seq_record in SeqIO.parse(ip, "fasta"):
        orfID = seq_record.id
        mo = reID.match(orfID)
        isoformID, start_codon, start_codon_pos, stop_codon, stop_codon_pos = mo.groups()
        seq_len = int(stop_codon_pos) - int(start_codon_pos)
        stop_codon_group = "%s_%s%s" % (isoformID, stop_codon, stop_codon_pos)
        if (not longest_per_isoform_and_stop.has_key(stop_codon_group) or
            longest_per_isoform_and_stop[stop_codon_group][1] < seq_len):
            longest_per_isoform_and_stop[stop_codon_group] = (seq_record, seq_len)

    return longest_per_isoform_and_stop


def writeReprProteins(longest_per_isoform_and_stop, output_fasta):
    op = open(output_fasta, 'w')
    for seq_record, seq_len in longest_per_isoform_and_stop.values():
        op.write(">%s\n%s\n" % (seq_record.id, str(seq_record.seq)))
    op.close()
    

if (__name__ == "__main__"):
    output_fasta = sys.argv[1]

    longest_per_isoform_and_stop = findLongestReprProteins(sys.stdin)
    sys.stderr.write("INFO: identified %d covering ORFs.\n" % len(longest_per_isoform_and_stop))
    writeReprProteins(longest_per_isoform_and_stop, output_fasta)

    sys.exit(0)
