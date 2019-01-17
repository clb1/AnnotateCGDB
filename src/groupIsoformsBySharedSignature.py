#!/usr/bin/env python

import gzip
import sys
from collections import defaultdict
from commonCGDB import associateIsoformsToGenes, buildIsoformPartsList, annotateAndOrderIsoformParts
from IsoformSignatureClasses import *
from Bio import SeqIO


def readIsoformSeqs(isoform_fasta):
    isoform_mRNAs = {}
    if (isoform_fasta.endswith(".gz")):
        ip = gzip.open(isoform_fasta, 'rb')
    else:
        ip = open(isoform_fasta, 'r')

    for seq_record in SeqIO.parse(ip, "fasta"):
        isoform_mRNAs[seq_record.id] = seq_record.seq.tostring().lower()

    ip.close()
    return isoform_mRNAs


def groupBySharedSignature(isoforms_per_gene, isoform_parts_list, isoform_seqs):
    isoforms_per_signature = defaultdict(set)
    num_amplicons_per_signature = {}
    for locus_isoforms in isoforms_per_gene.values():
        sequences_per_signature = defaultdict(set)
        isoforms_per_signature_per_locus = defaultdict(set)
        for isoform in locus_isoforms:
            if (not isoform_seqs.has_key(isoform)):
                print >> sys.stderr, "WARNING: no sequence information for %s." % isoform
            else:
                isg = IsoformSignatureGenerator(isoform, isoform_seqs[isoform], isoform_parts_list[isoform]["otp"], isoform_parts_list[isoform]["strand"])

                signatures = isg.getAllSignatures()
                for signature in signatures:
                    signature_tuple = signature.getSignatureTuple()
                    signature_seq = signature.getmRNASignatureSequence()
                    isoforms_per_signature_per_locus[signature_tuple].add(isoform)
                    sequences_per_signature[signature_tuple].add(signature_seq)

        for k,v in sequences_per_signature.items():
            num_amplicons_per_signature[k] = len(v)

        isoforms_per_signature.update(isoforms_per_signature_per_locus)

    return isoforms_per_signature, num_amplicons_per_signature


def writeOutput(isoforms_per_signature, num_amplicons_per_signature, output_file):
    op = open(output_file, 'w')
    op.write("#Signature\tNumber_of_Amplicons\tIsoforms_With_Signature\n")
    for signature, isoforms in isoforms_per_signature.items():
        signature_as_string = "%s_%d_%d_%s" % (signature[0], signature[1], signature[-2], signature[-1])
        num_amplicons = num_amplicons_per_signature[signature]
        isoforms_as_string = ",".join(isoforms)
        op.write("%s\t%d\t%s\n" % (signature_as_string, num_amplicons, isoforms_as_string))
    op.close()


if (__name__ == "__main__"):
    coding_isoform_fasta, noncoding_isoform_fasta, input_gtf, all_parts_bed, all_exonic_parts_tbl, output_file = sys.argv[1:]

    print >> sys.stderr, "INFO: associating isoforms to genes...",
    isoforms_per_gene = associateIsoformsToGenes(input_gtf)
    print >> sys.stderr, "done"

    print >> sys.stderr, "INFO: building isoforms' parts lists...",
    isoform_parts_list = buildIsoformPartsList(all_parts_bed)
    print >> sys.stderr, "done"

    print >> sys.stderr, "INFO: annotating and ordering isoform parts...",
    annotateAndOrderIsoformParts(all_exonic_parts_tbl, isoform_parts_list, True)
    print >> sys.stderr, "done"

    print >> sys.stderr, "INFO: reading isoform mRNA sequences...",
    isoform_seqs = readIsoformSeqs(coding_isoform_fasta)
    noncoding_isoform_seqs = readIsoformSeqs(noncoding_isoform_fasta)
    isoform_seqs.update(noncoding_isoform_seqs)
    print >> sys.stderr, "done"

    print >> sys.stderr, "INFO: grouping isoforms by shared signatures...",
    isoforms_per_signature, num_amplicons_per_signature = groupBySharedSignature(isoforms_per_gene, isoform_parts_list, isoform_seqs)
    print >> sys.stderr, "done"

    print >> sys.stderr, "INFO: ... writing output",
    writeOutput(isoforms_per_signature, num_amplicons_per_signature, output_file)
    print >> sys.stderr, "done"    

    sys.exit(0)
