#!/usr/bin/env python

import sys
from numpy import ceil, log2
from operator import itemgetter
from classMFE import MFE, MFEOptions
from IsoformSignatureClasses import Signature, IsoformSignatureGenerator
from commonCGDB import associateIsoformsToTargetGenesPlusOldID, buildTargetIsoformPartsList, annotateAndOrderTargetIsoformParts
from Bio import SeqIO
from IsoformSignatureClasses import *


def readTargetIsoformMRNAs(target_isoform_IDs, transcriptome_ref):
    target_IDs, target_isoforms = set(), {}

    with open(target_isoform_IDs, 'r') as ip:
        for line in ip:
            line = line.strip()
            target_IDs.add(line)
            
    for seq_record in SeqIO.parse(transcriptome_ref, "fasta"):
        if (seq_record.id in target_IDs):
            target_isoforms[seq_record.id] = str(seq_record.seq)

    assert (len(target_IDs) == len(target_isoforms.keys()))
    return target_isoforms


def writeSignatures(target_isoform_signature_with_primer_pair, isoformOID, output_file, output_gtf):
    op = open(output_file, 'w')
    op.write("IsoformID\tGeneSymbol\tStrand\tGenomicSignature\tTranscriptomicSignature\tAmpliconLength\tFP_Tm\tRP_Tm\tFP\tRP\tAmplicon\n")
    for isoformID, potential_signature in target_isoform_signature_with_primer_pair.items():
        fp, fp_tm, rp, rp_tm, strand = potential_signature.getPrimerPair()
        transcriptID = potential_signature.getParentTranscript()
        oID = isoformOID[transcriptID]
        assert (isoformID == transcriptID)
        signature_genomic = potential_signature.getSignatureTuple()
        genomic_signature_str = "%s:%d-%d_%s" % tuple(signature_genomic[0:2] + signature_genomic[-2:])
        signature_transcriptomic = potential_signature.getSignatureTupleInmRNACoords()
        transcriptomic_signature_str = "%s:%d-%d,%s" % tuple(signature_transcriptomic[0:2] + signature_transcriptomic[-2:])
        amplicon = potential_signature.getAmplicon()
        op.write("%s\t%s\t%s\t%s\t%s\t%d\t%4.2f\t%4.2f\t%s\t%s\t%s\n" % \
                 (transcriptID, oID, strand, genomic_signature_str, transcriptomic_signature_str, len(amplicon), fp_tm, rp_tm, fp, rp, amplicon))
    op.close()

    print >> sys.stderr, "INFO: transcripts for which primer pairs were found:"
    gtf = open(output_gtf, 'w')
    gtf.write("track type=gtf name=Amplicons description=Amplicons visibility=full useScore=1\n")
    for potential_signature in target_isoform_signature_with_primer_pair.values():
        print >> sys.stderr, "\t%s" % potential_signature.parent_transcript
        gtf_lines = potential_signature.getAmpliconGenomicCoords("gtf")
        gtf.write("%s\n" % "\n".join(gtf_lines))
    gtf.close()


def designPrimers(target_isoforms, isoforms_per_gene, isoform_parts_list, genome_ref, transcriptome_ref, fwd_tag, rev_tag, do_in_parallel=False):
    primer3_param_sets = primer32ParamSets(1)
    target_isoform_IDs = set(target_isoforms.keys())
    target_isoform_signature_with_primer_pair = {}

    transcriptomic_MFE = MFE(transcriptome_ref)
    genomic_MFE = MFE(genome_ref)
    mfe_options = MFEOptions()
    num_isoforms_no_signature, num_isoforms_signature_but_no_primers = (0,0)

    # Process each locus/gene in turn, but only need to find a signature for the isoforms in 'target_isoform_IDs'
    for gene, locus_isoforms in isoforms_per_gene.items():
        all_signature_generators = {}
        for isoform in locus_isoforms:
            isoform_mRNA_seq = None if not target_isoforms.has_key(isoform) else target_isoforms[isoform]
            all_signature_generators[isoform] = IsoformSignatureGenerator(isoform, isoform_mRNA_seq, isoform_parts_list[isoform]["otp"], isoform_parts_list[isoform]["strand"])

        this_locus_target_isoform_IDs = locus_isoforms.intersection(target_isoform_IDs)

        # For each isoform in turn, find a signature, see check for best primers, then check for genomic specificity.  Add to
        # class IsoformSignatureGenerator some clean way to say "exhausted signature options" so that I will know that what
        # I've got to that point is the best I'll get to work with.
        for target_isoform_ID in this_locus_target_isoform_IDs:
            #if (target_isoform_ID not in ["ERCC-00134"]):
            #    continue
            num_valid_signatures_for_isoform = 0
            print >> sys.stderr, "\nINFO: finding signature with primers for %s" % target_isoform_ID
            the_other_isoforms = locus_isoforms - set([target_isoform_ID])

            if (do_in_parallel):
                q = Queue()
                all_unique_signatures, jobs = [], []
                for signature in all_signature_generators[target_isoform_ID].getAllSignatures():
                    if (not any([all_signature_generators[x].signatureIndicatesMe(signature) for x in the_other_isoforms])):
                        num_valid_signatures_for_isoform += 1
                        all_unique_signatures.append(signature)
                        job = Process(target=signature.findPrimerPair, args=(primer3_param_sets, transcriptomic_MFE, genomic_MFE, mfe_options, set([target_isoform_ID]), False, q))
                        job.start()
                        print >> sys.stderr, "Started job"
                        jobs.append(job)

                if (len(all_unique_signatures)>0):
                    for job in jobs:
                        job.join()
                    primer_results = dict([q.get() for job in jobs])
                    signature_results = []
                    for signature in all_unique_signatures:
                        signature_tuple = signature.getSignatureTuple()
                        if (primer_results[signature_tuple] != None):
                            signature.setPrimerPair(  * primer_results[signature_tuple] )
                            signature_results.append( (signature.getPrimerPairStratum(), signature.getAmpliconLength(), signature) )
                    if (len(signature_results) > 0):
                        signature_results = sorted(signature_results, key=itemgetter(0,1))
                        target_isoform_ID_signature_with_primer_pair[target_isoform_ID] = signature_results[0][-1]
            else:
                potential_signature = all_signature_generators[target_isoform_ID].nextPotentialSignature()
                best_primer_pair_stratum = 1e10
                while (potential_signature != None): # and best_primer_pair_stratum > 0): #  (not potential_signature.hasPrimerPair() or
                    if (not any([all_signature_generators[x].signatureIndicatesMe(potential_signature) for x in the_other_isoforms])):
                        num_valid_signatures_for_isoform += 1
                        print >> sys.stderr, "\tEvaluating signature", potential_signature.signature_tuple
                        msg = potential_signature.findPrimerPair(primer3_param_sets, transcriptomic_MFE, genomic_MFE, mfe_options, target_isoform_ID,
                                                                 fwd_tag, rev_tag, True, True)
                        if (potential_signature.hasPrimerPair()):
                            if (target_isoform_signature_with_primer_pair.has_key(target_isoform_ID)):
                                if (target_isoform_signature_with_primer_pair[target_isoform_ID].getPrimerPairStratum() > potential_signature.getPrimerPairStratum()):
                                    print >> sys.stderr, msg, "Is new best because primers from lower parameter stratum."
                                    target_isoform_signature_with_primer_pair[target_isoform_ID] = potential_signature
                                    best_primer_pair_stratum = potential_signature.getPrimerPairStratum()
                                elif (target_isoform_signature_with_primer_pair[target_isoform_ID].getAmpliconLength() > potential_signature.getAmpliconLength() and
                                      target_isoform_signature_with_primer_pair[target_isoform_ID].getPenalty() > potential_signature.getPenalty()):
                                    print >> sys.stderr, msg, "Is new best because amplicon is shorter and has smaller Primer3 penalty."
                                    target_isoform_signature_with_primer_pair[target_isoform_ID] = potential_signature
                                    best_primer_pair_stratum = potential_signature.getPrimerPairStratum()
                            else:
                                print >> sys.stderr, msg
                                target_isoform_signature_with_primer_pair[target_isoform_ID] = potential_signature
                                best_primer_pair_stratum = potential_signature.getPrimerPairStratum()
                    potential_signature = all_signature_generators[target_isoform_ID].nextPotentialSignature()

            if (not target_isoform_signature_with_primer_pair.has_key(target_isoform_ID)):
                if (num_valid_signatures_for_isoform > 0):
                    num_isoforms_signature_but_no_primers += 1
                    print >> sys.stderr, "WARNING: could not find suitable primer pair for any of the %d signature(s) of %s" % \
                      (num_valid_signatures_for_isoform, target_isoform_ID)
                else:
                    num_isoforms_no_signature += 1
                    print >> sys.stderr, "WARNING: could not find distinguishing signature for %s" % target_isoform_ID

    transcriptomic_MFE.finalize()
    genomic_MFE.finalize()

    print >> sys.stderr, "\nINFO: of %d isoforms, %d had no signatures and %d had no signature primers" % \
        (len(target_isoforms), num_isoforms_no_signature, num_isoforms_signature_but_no_primers)

    return target_isoform_signature_with_primer_pair


def evalPrimerPool(target_isoform_signature_with_primer_pair):
    all_primers = []
    for isoformID, potential_signature in target_isoform_signature_with_primer_pair.items():
        fp, fp_tm, rp, rp_tm, strand = potential_signature.getPrimerPair()
        all_primers.extend( [(isoformID, "fwd", fp), (isoformID, "rev", rp)] )

    OligosEvaluator = DNAStructureEvaluator(56.0)

    # Evaluate all pairs of primers (excluding ones from same isoformID)
    all_primer_pairs = []
    for i in xrange(len(all_primers)-1):
        for j in xrange(i+1, len(all_primers)):
            if (all_primers[i][0] != all_primers[j][0]):
                pair_label = "%s,%s+%s,%s" % (all_primers[i][0], all_primers[i][1], all_primers[j][0], all_primers[j][1])
                all_primer_pairs.append( (pair_label, (all_primers[i][2], all_primers[j][2])) )

    all_primer_pair_dG = OligosEvaluator.checkForDimer(all_primer_pairs, True)

    

def primer32ParamSets(num_ranges=5):
    common = [("PRIMER_TASK","generic"), ("PRIMER_EXPLAIN_FLAG",1), ("PRIMER_FIRST_BASE_INDEX",1), ("PRIMER_OPT_SIZE",18), 
              ("PRIMER_MIN_SIZE",18), ("PRIMER_MAX_SIZE",24), ("PRIMER_PRODUCT_OPT_SIZE",150), ("PRIMER_PRODUCT_SIZE_RANGE","90-210"),
              ("PRIMER_PAIR_MAX_DIFF_TM",6), ("PRIMER_THERMODYNAMIC_ALIGNMENT",1)] # , ("SEQUENCE_INCLUDED_REGION","1,149")

    param_sets = [ dict( common + [("PRIMER_MIN_TM",58), ("PRIMER_MAX_TM",63), ("PRIMER_OPT_TM",60), ("PRIMER_SALT_DIVALENT",2.5), ("PRIMER_DNTP_CONC",0.8)] ),
                   dict( common + [("PRIMER_MIN_TM",63), ("PRIMER_MAX_TM",66), ("PRIMER_OPT_TM",63), ("PRIMER_SALT_DIVALENT",2.5), ("PRIMER_DNTP_CONC",0.8)] ),
                   dict( common + [("PRIMER_MIN_TM",54), ("PRIMER_MAX_TM",57), ("PRIMER_OPT_TM",57), ("PRIMER_SALT_DIVALENT",2.5), ("PRIMER_DNTP_CONC",0.8)] ),
                   dict( common + [("PRIMER_MIN_TM",66), ("PRIMER_MAX_TM",69), ("PRIMER_OPT_TM",66), ("PRIMER_SALT_DIVALENT",2.5), ("PRIMER_DNTP_CONC",0.8)] ),
                   dict( common + [("PRIMER_MIN_TM",51), ("PRIMER_MAX_TM",54), ("PRIMER_OPT_TM",54), ("PRIMER_SALT_DIVALENT",2.5), ("PRIMER_DNTP_CONC",0.8)] ) ]

    return param_sets[0:num_ranges]

                    
if (__name__ == "__main__"):
    fwd_tag, rev_tag, target_isoform_IDs, genome_ref, transcriptome_ref, input_gtf, all_parts_bed, all_exonic_parts_tbl, output_file, output_gtf = sys.argv[1:]
    

    target_isoforms = readTargetIsoformMRNAs(target_isoform_IDs, transcriptome_ref) # target_isoform_mRNAs
    isoforms_per_gene, all_target_loci_isoforms, isoformOID = associateIsoformsToTargetGenesPlusOldID(input_gtf, target_isoforms)
    isoform_parts_list = buildTargetIsoformPartsList(all_parts_bed, all_target_loci_isoforms)
    annotateAndOrderTargetIsoformParts(all_exonic_parts_tbl, isoform_parts_list, all_target_loci_isoforms, True)
    target_isoform_signature_with_primer_pair = designPrimers(target_isoforms, isoforms_per_gene, isoform_parts_list, genome_ref, transcriptome_ref, fwd_tag, rev_tag)

    #cPickle.dump(target_isoform_signature_with_primer_pair, open("signature_primers.pkl","wb"))
    #target_isoform_signature_with_primer_pair = cPickle.load(open("signature_primers.pkl","rb"))

    evalPrimerPool(target_isoform_signature_with_primer_pair)

    writeSignatures(target_isoform_signature_with_primer_pair, isoformOID, output_file, output_gtf)
    
    sys.exit(0)
    

