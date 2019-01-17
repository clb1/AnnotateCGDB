#!/usr/bin/env python

import bz2
from collections import defaultdict
from itertools import product
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re
import sys


def identifyORFs(nuc_seq, start_codons, stop_codons, min_ORF_len):
    """Expects 'nuc_seq' to be in 5'->3' orientation.
       Identifies and returns: 

          1) all ORFs containing standard start and stop codons, 
          2) partial ORFs for stop codons with no start codon,
          3) partial ORFs for start codons with no stop codon,
          4) for frames with no complete or partial ORF, the partial ORF
          5) the longest ORF or partial ORF for each stop codon (or last complete codon in the case of partial ORFs)
          subject to an ORF being >= min_ORF_len
    """
    orfs = []
    starts_w_stop = set()
    
    # Record start & stop codons with tuples of the form (codon triplet, start index, end index, frame)
    starts, stops = [], []
    for sc in start_codons:
        starts.extend( [(sc, a.start(), a.end(), a.start()%3) for a in re.finditer(sc, nuc_seq)] )
    starts = sorted(starts, key=itemgetter(1))
    
    for sc in stop_codons:
        stops.extend( [(sc, a.start(), a.end(), a.start()%3) for a in re.finditer(sc, nuc_seq)] )
    stops = sorted(stops, key=itemgetter(1))
    
    for i, stop_tup in enumerate(stops):
        frame = stop_tup[3]
        upstream_stops = filter(lambda y:y[3]==frame, stops[0:i])
        if (len(upstream_stops) > 0):
            upstream_starts = filter(lambda x:x[1]<=stop_tup[1] and x[3]==frame and x[1]>upstream_stops[-1][1], starts)
        else:
            upstream_starts = filter(lambda x:x[1]<=stop_tup[1] and x[3]==frame, starts)
            
        if (len(upstream_starts) > 0):
            starts_w_stop.update(upstream_starts)
            most_5prime_start = upstream_starts[0]
            if (stop_tup[1]-most_5prime_start[1] >= min_ORF_len):
                coding_rna = Seq(nuc_seq[most_5prime_start[1]:stop_tup[2]], IUPAC.ambiguous_rna)
                orfs.append( (most_5prime_start, stop_tup, coding_rna, True) )
                
            if (len(upstream_starts) > 1):                
                for upstream_start in filter(lambda y: stop_tup[1]-y[1] >= min_ORF_len, upstream_starts[1:]):
                    coding_rna = Seq(nuc_seq[upstream_start[1]:stop_tup[2]], IUPAC.ambiguous_rna)
                    orfs.append( (upstream_start, stop_tup, coding_rna, False) )

        else:
            # Identify partial ORFs with a stop codon but no start codon
            partial_pORF_start_tup = ('xxx', frame, frame+3, frame)
            if (stop_tup[1]-frame >= min_ORF_len and not any(map(lambda y:y[3]==stop_tup[3], stops[0:i]))):
                coding_rna = Seq(nuc_seq[frame:stop_tup[2]], IUPAC.ambiguous_rna)
                orfs.append( (partial_pORF_start_tup, stop_tup, coding_rna, True) ) 

    # Identify partial ORFs with a start codon but no stop codon
    partial_orfs_for_codon = defaultdict(list)
    for start_tup in filter(lambda x: not x in starts_w_stop, starts):
        frame = start_tup[3]
        start_of_last_complete_codon = filter(lambda y:(y - start_tup[1])%3==0, range(len(nuc_seq)-5, len(nuc_seq)-2))[0]
        if (start_of_last_complete_codon - frame >= min_ORF_len):
            partial_pORF_stop_tup = ("x%dx" % frame, start_of_last_complete_codon, start_of_last_complete_codon+3, frame)
            coding_rna = Seq(nuc_seq[start_tup[1]:start_of_last_complete_codon+3], IUPAC.ambiguous_rna)
            partial_orfs_for_codon[partial_pORF_stop_tup].append( (start_tup, partial_pORF_stop_tup, coding_rna) )
    
    for assoc_start_tups in partial_orfs_for_codon.values():
        assoc_start_tups = sorted(assoc_start_tups, key=itemgetter(1))
        orfs.append( assoc_start_tups[0] + (True,) )
        if (len(assoc_start_tups) > 1):
            for tup in assoc_start_tups[1:]:
                orfs.append( tup + (False,) )
    
    # Check that all three frames have an associated start or stop. If not, then it's a partial, mid-ORF piece of sequence.
    frames_wo_orf = set([0,1,2]) - set(map(lambda y:y[3], starts)) - set(map(lambda y:y[3], stops))
    #print frames_wo_orf
    for frame in frames_wo_orf:
        partial_pORF_start_tup = ('xxx', frame, frame+3, frame)
        start_of_last_complete_codon = filter(lambda y:(y - frame)%3==0, range(len(nuc_seq)-5, len(nuc_seq)-2))[0]
        if (start_of_last_complete_codon - frame >= min_ORF_len):
            partial_pORF_stop_tup = ("x%dx" % frame, start_of_last_complete_codon, start_of_last_complete_codon+3, frame)                        
            coding_rna = Seq(nuc_seq[frame:start_of_last_complete_codon+3], IUPAC.ambiguous_rna)
            orfs.append( (partial_pORF_start_tup, partial_pORF_stop_tup, coding_rna, True) ) 
            
    return orfs


def identifyAndTranslateAllORFs(mRNA_fasta_file, start_codons, stop_codons, min_protein_len, short_protein_len):
    min_ORF_len = 3 * min_protein_len

    short_protein_fastas = {}
    long_protein_fastas = {}
    for start_codon in start_codons:
        short_protein_fastas[start_codon] = bz2.BZ2File("CG-%s-pORFs-lt%d.fa.bz2" % (start_codon, short_protein_len), 'wb')
        long_protein_fastas[start_codon] = bz2.BZ2File("CG-%s-pORFs-ge%d.fa.bz2" % (start_codon, short_protein_len), 'wb')

    longest_long_protein_fasta = bz2.BZ2File("CG-longest-pORFs-ge%d.fa.bz2" % short_protein_len, 'wb')
    longest_short_protein_fasta = bz2.BZ2File("CG-longest-pORFs-lt%d.fa.bz2" % short_protein_len, 'wb')
    
    num_short, num_long = (0,0)
    incomplete_orf_stops = set(["x0x","x1x","x2x"])
    for counter, seq_record in enumerate(SeqIO.parse(mRNA_fasta_file, "fasta")):
        isoform_ID = seq_record.id
        if (counter%1000==0):
            print >> sys.stderr, "INFO: Processed %d" % counter

        nuc_seq = str(seq_record.seq).lower()
        orfs = identifyORFs(nuc_seq, start_codons, stop_codons, min_ORF_len)
        for start_codon_tup, stop_codon_tup, coding_rna, is_longest in orfs:
            protein = coding_rna.translate(table=1)
            protein_str = str(protein)
            if (start_codon_tup[0] != "xxx" and protein_str[0] != 'M'):
                protein_str = "M" + protein_str[1:]

            if (protein_str.count('*') > 1 or (protein_str[-1] != '*' and not stop_codon_tup[0] in incomplete_orf_stops)):
                print >> sys.stderr, "ERROR: protein translation problem for %s: %s" % (isoform_ID, protein_str), start_codon_tup, stop_codon_tup
                pdb.set_trace()
                sys.exit(1)

            if (protein_str[-1] == "*"):
                protein_str = protein_str[0:-1]
                assert (protein_str[-1] != "*")

            protein_id_and_seq = ">%s,%s%d-%s%d\n%s\n" % (isoform_ID, start_codon_tup[0], start_codon_tup[1], stop_codon_tup[0], stop_codon_tup[1], protein_str)
            if (len(protein_str) < short_protein_len):
                num_short += 1
                short_protein_fastas[start_codon_tup[0]].write(protein_id_and_seq)
                if (len(protein_str) >= min_protein_len and is_longest):
                    longest_short_protein_fasta.write(protein_id_and_seq)
            else:
                num_long += 1
                long_protein_fastas[start_codon_tup[0]].write(protein_id_and_seq)
                if (is_longest):
                    longest_long_protein_fasta.write(protein_id_and_seq)

    for op in short_protein_fastas.values():
        op.close()
    for op in long_protein_fastas.values():
        op.close()
    longest_long_protein_fasta.close()
    longest_short_protein_fasta.close()
    
    #print >> sys.stderr, "INFO: wrote %s and %s with %d and %d pORFs, respectively for start codon %s." % \
    #    (short_protein_fasta_file, long_protein_fasta_file, num_short, num_long, start_codon)

    
if (__name__ == "__main__"):
    short_protein_len, mRNA_fasta_file = sys.argv[1:]
    short_protein_len = int(short_protein_len)
    
    start_codons = ["atg", "ctg", "gtg", "ttg", "acg", "xxx"] # As cDNA
    stop_codons = ["taa","tag","tga", "x0x", "x1x", "x2x"] # The 'x#x' are for non-stop codons in incomplete ORF of reading frames 0,1,2
    min_protein_len = 10
    
    identifyAndTranslateAllORFs(mRNA_fasta_file, start_codons, stop_codons, min_protein_len, short_protein_len)

    sys.exit(0)
