#!/usr/bin/env python

import sys
from numpy import abs, argmin, array, mean, std
from collections import defaultdict, Counter
from itertools import chain
from operator import itemgetter
from classMFE import MFE, MFEOptions
from Bio import SeqIO
from IsoformSignatureClasses import Signature, IsoformSignatureGenerator
from commonCGDB import buildIsoformPartsList, annotateAndOrderIsoformParts

import designParams
from designSetsClasses import OligoThermodynamics


def findIsoformGroupSignatures(target_isoform_groups, isoform_group_defs, isoforms_per_signature_file):
    signatures_per_isoform_group = None

    # Get the isoform group definitions for the target isoform groups
    target_IGs = []
    with open(target_isoform_groups,'r') as ip:
        for line in ip:
            IG_ID = line.strip()
            target_IGs.append(IG_ID)

    isoforms_to_IG_map = {}
    IG_to_isoforms_map = {}
    IG_ID_to_gene_symbol = {}
    with open(isoform_group_defs,'r') as ip:
        header = ip.readline()
        for line in ip:
            IG_ID, combined_isoformIDs, gene_symbol = line.strip().split("\t")[0:3]
            if (IG_ID in target_IGs):
                isoformIDs = tuple(sorted(combined_isoformIDs.split("|")))
                isoforms_to_IG_map[isoformIDs] = IG_ID
                IG_to_isoforms_map[IG_ID] = set(isoformIDs) # isoformIDs
                IG_ID_to_gene_symbol[IG_ID] = gene_symbol

    # Find the signatures that have the exact set of isoforms for some target isoform group
    IG_signatures = defaultdict(set)
    ip = open(isoforms_per_signature_file, 'r')
    header = ip.readline()
    for line in ip:
        signature, num_amplicons, isoforms_csv = line.strip().split("\t")
        if (int(num_amplicons) == 1):
            isoformIDs = tuple(sorted(isoforms_csv.split(",")))
            if (isoforms_to_IG_map.has_key(isoformIDs)):
                IG_signatures[ isoforms_to_IG_map[isoformIDs] ].add(signature)
    ip.close()

    return (IG_signatures, IG_to_isoforms_map, IG_ID_to_gene_symbol, target_IGs)


def readTargetIsoformMRNAs(target_isoform_IDs, isoform_fasta):
    target_isoforms = {}

    for seq_record in SeqIO.parse(isoform_fasta, "fasta"):
        if (seq_record.id in target_isoform_IDs):
            target_isoforms[seq_record.id] = str(seq_record.seq)

    assert (len(target_isoform_IDs) == len(target_isoforms.keys()))
    return target_isoforms


def designPrimers(IG_signatures, IG_to_isoforms_map, isoform_parts_list, isoform_seqs, transcriptomic_MFE, genomic_MFE, oligo_thermo):
    """Finds one primer pair per isoform group."""
    IG_signature_primers = {}

    short_enough_len = 90
    
    IG_with_primered_signature = set()
    mfe_options = MFEOptions()

    for IG_ID, candidate_signatures in IG_signatures.items():
        print >> sys.stderr, "\nINFO: finding signature with primers for %s" % IG_ID
        for candidate_signature in candidate_signatures:
            isoforms_of_signature = IG_to_isoforms_map[IG_ID]
            isoform = list(isoforms_of_signature)[0] # Any isoform will suffice, since primers must target all (and only) isoforms of the signature TODO: choose shortest
            isg = IsoformSignatureGenerator(isoform, isoform_seqs[isoform], isoform_parts_list[isoform]["otp"], isoform_parts_list[isoform]["strand"])
            chromosome, s, e, part_type = candidate_signature.split("_")
            signature_candidate_tuple = (chromosome, int(s), int(e), part_type)
            S = isg.lookupSignature(signature_candidate_tuple)

            # TODO: send findPrimerPair() penalty and amplicon length of best solution so far and only set a primer pair if it can find one equal/better on both criteria
            msg = S.findPrimerPair(transcriptomic_MFE, genomic_MFE, mfe_options, isoforms_of_signature, oligo_thermo)
            if (S.hasPrimerPair()):
                if (IG_signature_primers.has_key(IG_ID)):
                    if (IG_signature_primers[IG_ID].getAmpliconLength() > S.getAmpliconLength()):
                        print >> sys.stderr, msg, "Is new best because amplicon is shorter."                        
                        IG_signature_primers[IG_ID] = S

                    if (S.getAmpliconLength() <= short_enough_len):
                        break
                else:
                    print >> sys.stderr, msg
                    IG_signature_primers[IG_ID] = S
                
        if (not IG_signature_primers.has_key(IG_ID)):
            print >> sys.stderr, "WARNING: could not find suitable primer pair for any of the %d signature(s) of %s" % (len(candidate_signatures), IG_ID)

    return IG_signature_primers


def writeSignatures(ordered_target_isoform_groups, IG_signature_primers, IG_ID_to_gene_symbol, output_file, output_gtf):
    op = open(output_file, 'w')
    op.write("GroupID\tGeneSymbol\tStrand\tGenomicSignature\tAmpliconLength\tFP_Tm\tRP_Tm\tFP\tRP\tAmplicon\n")
    for IG_ID in ordered_target_isoform_groups:
        if (IG_signature_primers.has_key(IG_ID)):
            S = IG_signature_primers[IG_ID]
            fp, fp_tm, rp, rp_tm, strand = S.getPrimerPair()
            gene_symbol = IG_ID_to_gene_symbol[IG_ID]
            signature_genomic = S.getSignatureTuple()
            genomic_signature_str = "%s:%d-%d_%s" % tuple(signature_genomic[0:2] + signature_genomic[-2:])
            amplicon = S.getAmplicon()
            op.write("%s\t%s\t%s\t%s\t%d\t%4.2f\t%4.2f\t%s\t%s\t%s\n" % \
                     (IG_ID, gene_symbol, strand, genomic_signature_str, len(amplicon), fp_tm, rp_tm, fp, rp, amplicon))
    op.close()

    gtf = open(output_gtf, 'w')
    gtf.write("track type=gtf name=Amplicons description=Amplicons visibility=full useScore=1\n")
    for S in IG_signature_primers.values():
        gtf_lines = S.getAmpliconGenomicCoords("gtf")
        gtf.write("%s\n" % "\n".join(gtf_lines))
    gtf.close()

    print >> sys.stderr, "\nINFO: wrote files for primers found for %d of %d target isoform groups." % (len(IG_signature_primers.keys()), len(ordered_target_isoform_groups))


if (__name__ == "__main__"):
    target_isoform_groups, isoform_group_defs, isoforms_per_signature_file, hg_fasta, mRNA_isoform_fasta, input_gtf, all_parts_bed, all_exonic_parts_bed, output_file, output_gtf = sys.argv[1:]

    print >> sys.stderr, "INFO: finding signatures for isoform groups...",
    IG_signatures, IG_to_isoforms_map, IG_ID_to_gene_symbol, ordered_target_isoform_groups = findIsoformGroupSignatures(target_isoform_groups, isoform_group_defs, isoforms_per_signature_file)
    print >> sys.stderr, "done"
    print >> sys.stderr, "INFO: found at least one signature for %d of the %d isoform groups" % (len(IG_signatures.keys()), len(ordered_target_isoform_groups))
    
    print >> sys.stderr, "INFO: building isoforms' parts lists...",
    isoform_parts_list = buildIsoformPartsList(all_parts_bed)
    print >> sys.stderr, "done"

    print >> sys.stderr, "INFO: annotating and ordering isoform parts...",
    annotateAndOrderIsoformParts(all_exonic_parts_bed, isoform_parts_list, True)
    print >> sys.stderr, "done"

    print >> sys.stderr, "INFO: reading isoform mRNA sequences...",
    target_isoform_IDs = set([x for x in chain( * IG_to_isoforms_map.values() )])
    isoform_seqs = readTargetIsoformMRNAs(target_isoform_IDs, mRNA_isoform_fasta)
    print >> sys.stderr, "done"

    transcriptomic_MFE = MFE(mRNA_isoform_fasta)
    genomic_MFE = MFE(hg_fasta)
    
    designParams.setParameters(param_settings_descriptor)
    oligo_thermo = OligoThermodynamics()

    print >> sys.stderr, ""
    IG_signature_primers = designPrimers(IG_signatures, IG_to_isoforms_map, isoform_parts_list, isoform_seqs, transcriptomic_MFE, genomic_MFE, oligo_thermo)
    
    transcriptomic_MFE.finalize()
    genomic_MFE.finalize()
    
    print >> sys.stderr, "INFO: writing reference candidates to %s..." % output_file,
    writeSignatures(ordered_target_isoform_groups, IG_signature_primers, IG_ID_to_gene_symbol, output_file, output_gtf)
    print >> sys.stderr, "done"
    
    sys.exit(0)
