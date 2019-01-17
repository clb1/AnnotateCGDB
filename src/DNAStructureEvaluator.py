import re
import RNA
from subprocess import Popen, PIPE
import sys


class DNAStructureEvaluator:
    '''A class evaluating the hairpin, and homo-/heterodimer stability of DNA primers.'''

    def __init__(self, temperature):
        self.temperature = temperature
        self.RNAfold_pattern = re.compile("(\S+)\s+.\s+(-?\d+\.\d\d)\s*")
        #self.RNAcofold_pattern = re.compile("(\S+)&(\S+)\s+.\s+(-?\d+\.\d\d)")
        self.RNAcofold_pattern = re.compile("(\S+)&(\S+)\s+\[\s*(-?\d+\.\d\d)\]")

    def setTemperature(self, new_temperature):
        self.temperature = new_temperature


    def checkForHairpins(self, seqID_and_seq_tuples, verbose=False):
        '''Failure by hairpin occurs if there is a stem > 3bp or deltaG < -2.0 kcal/mol.
           seqID_and_seq_tuples is [(seqID,seq),...]
        '''
        hairpin_deltaG = {}

        cmd = "RNAfold -p -T %f --noPS --noconv" % self.temperature
        RNAfold = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)        

        input_string = ""
        for seqID, seq in seqID_and_seq_tuples:
            input_string += ">%s\n%s\n" % (seqID, seq)

        stdout_results, stderr_results = RNAfold.communicate(input=input_string)

        for seqID, RNAfold_lines in map(lambda y: (y[0],y[1:]), map(lambda x: x.split("\n")[0:-1], stdout_results.split(">")[1:])):
            seqID = int(seqID)
            hairpin_deltaG[seqID] = None
            #nuc_seq = RNAfold_lines[0]
            if (verbose):
                print >> sys.stderr, "\n###\n", "\n".join(RNAfold_lines)
            seq_has_hairpin = False
            for struct_line in RNAfold_lines[1:-1]:
                mo = self.RNAfold_pattern.match(struct_line)
                struct, dG = mo.groups()
                dG = float(dG)
                struct_line_segments = re.split("\.|,", struct)
                seq_has_hairpin = any(map(lambda z: len(z)>3, struct_line_segments)) or dG < -2.0
                if (verbose):
                    print >> sys.stderr, struct_line, seq_has_hairpin
                if (seq_has_hairpin):
                    hairpin_deltaG[seqID] = None
                    break
                elif (hairpin_deltaG[seqID] == None):
                    hairpin_deltaG[seqID] = dG
                else:
                    hairpin_deltaG[seqID] = min(hairpin_deltaG[seqID], dG)

        if (verbose):
            print >> sys.stderr, "%d of %d passed hairpin check" % (len(filter(lambda x: x!=None, hairpin_deltaG.values())), len(hairpin_deltaG))

        return hairpin_deltaG


    def checkForDimer(self, seqID_and_seqs_tuples, verbose=False):
        '''seqID_and_seqs_tuples is [(seqID, (seq1, seq2)),...]'''
        primer_deltaG = {}

        input_string = ""
        for counter, (seqID, seqs) in enumerate(seqID_and_seqs_tuples):
            input_string += ">%s\n%s&%s\n" % (counter, seqs[0], seqs[1]) # seqID

        cmd = "RNAcofold -p -T %f --noPS --noconv" % self.temperature
        RNAcofold = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)        

        stdout_results, stderr_results = RNAcofold.communicate(input=input_string)
        
        counter = 0
        for seqID, RNAcofold_lines in map(lambda y: (y[0],y[1:]), map(lambda x: x.split("\n")[0:-1], stdout_results.split(">")[1:])):
            seqID = int(seqID)
            primer_deltaG[seqID] = None

            assert (len(RNAcofold_lines) == 4)
            seq1, seq2 = RNAcofold_lines[0].strip().split("&")
            mo = self.RNAcofold_pattern.match(RNAcofold_lines[2])
            struct1, struct2, ensemble_dG = mo.groups()
            dG = float(ensemble_dG)
            if (dG >= -6.0):
                struct1_segments = re.split("\.|,", struct1)
                struct2_segments = re.split("\.|,", struct2)
                len_struct1_3p_dimer = len(struct1_segments[-1])
                len_struct2_3p_dimer = len(struct2_segments[-1])

                seq1_has_3p_dimer = (len_struct1_3p_dimer>3 or
                                     (len_struct1_3p_dimer>1 and (seq1[-len_struct1_3p_dimer:].count("C") + seq1[-len_struct1_3p_dimer:].count("G") > 1 and dG <= -5.0)))

                seq2_has_3p_dimer = (len_struct2_3p_dimer>1 or
                                     (len_struct2_3p_dimer>1 and (seq2[-len_struct2_3p_dimer:].count("C") + seq2[-len_struct2_3p_dimer:].count("G") > 1 and dG <= -5.0)))

                seq1_has_internal_dimer = any(map(lambda x: x>5, map(len, struct1_segments[0:-1]))) #and dG <= -6.0
                seq2_has_internal_dimer = any(map(lambda x: x>5, map(len, struct2_segments[0:-1]))) #and dG <= -6.0

                seq_is_disqualified = (seq1_has_3p_dimer or seq2_has_3p_dimer or seq1_has_internal_dimer or seq2_has_internal_dimer)
                if (verbose):
                    print >> sys.stderr, counter, ": ", struct1, struct2, dG, "No Pass" if (seq_is_disqualified) else "Pass"

                if (seq_is_disqualified):
                    primer_deltaG[seqID] = None
                else:
                    primer_deltaG[seqID] = dG
            elif (verbose):
                    print >> sys.stderr, counter, ": Skipped, deltaG = ", dG
            counter += 1

        if (verbose):
            print >> sys.stderr, "%d of %d passed dimer check" % (len(filter(lambda x: x!=None, primer_deltaG.values())), len(primer_deltaG))

        return primer_deltaG
