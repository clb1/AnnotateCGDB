#!/usr/bin/env python
from __future__ import division

import sys
from collections import defaultdict
from operator import itemgetter
import math
import sqlite3
import twobitreader
sys.path.append("/usr/local/src/MFEprimer")
sys.path.append("/usr/local/src/MFEprimer/chilli")

from MFEprimer import get_pos_data, get_pos_range, get_align_seq, cal_PPC, draw_graphical_alignment_primer

from chilli import Seq
from chilli import TmDeltaG
from chilli import chilli

from string import maketrans, translate
DNA_complement_table = maketrans("ACGTNacgtn","TGCANtgcan")

import pdb

class MFEOptions:
    def __init__(self, size_start=50, size_stop=800, k=8):
        #self.amplicon=False
        #self.database=""
        self.dg_start=-9223372036854775807
        self.dg_stop=0
        self.diva_conc=1.5
        self.dntp_conc=0.25
        self.k_value=k
        self.mono_conc=50
        self.oligo_conc=50
        self.ppc=30
        self.size_start=size_start
        self.size_stop=size_stop
        self.tm_start=40.0
        self.tm_stop=70.0
    

class MFE:
    def __init__(self, db, options):
        self.options = options

        fcdict_cache = db + '.uni'
        self.fcdict = chilli.get_cache(fcdict_cache)

        dbname = db + '.sqlite3.db'
        self.conn = sqlite3.connect(dbname)
        self.cur = self.conn.cursor()

        self.db_two_bit = twobitreader.TwoBitFile(db + '.2bit')
        
        # These are set with call to processPrimer()
        self.fwd_primer_seq = None
        self.rev_primer_seq = None
        self.amp_list = None

        
    def clear(self):
        self.fwd_primer_seq = None
        self.rev_primer_seq = None
        self.amp_list = None


    def finalize(self):
        self.cur.close()
        self.conn.close()
        

    def getAmpliconData(self):
        # See MFEprimer.py:tab_out() for more details
        return self.amp_list

    
    def getPrimerProducts(self, apply_directionality_logic=True, max_report_length=100000):
        '''Directionality here is to be understood in the context of an RT-qSeq experiment.

           In a non-directional experiment, fwd and rev primers are used together in a 1-step 1st strand + 2nd strand cDNA synthesis reaction.
           After generating the sequencing library from the cDNA and sequencing, it is not possible to determine whether the fwd or the rev
           primer generated the 1st cDNA strand. So regardless of the chromosome strand of the target isoform, the sequencing results are the
           same. For non-directional experiments, this method returns each product in a standardized orientation such that the product begins
           with the fwd primer.
           
           In a directional experiment, rev primers are used to generate the 1st cDNA strand and the fwd primers are used to generate the 2nd
           cDNA strand in a 2-step cDNA synthesis reaction. Thus any product returned by MFEprimer that begins with the rev primer is not a
           product that could have been generated in the directional experiment. These MFEprimer products are ignored and not returned if 
           apply_directionality_logic is 'True'.

           mrna_hyb_start_pos and mrna_hyb_stop_pos are both in 0-based coordinates (on the mRNA template). I think...
        '''
        assert(isinstance(apply_directionality_logic, bool) and isinstance(max_report_length,int)), "Parameter types mixed up in classMFE:getPrimerProducts()"

        primer_products = defaultdict(list)
        for amp in map(itemgetter(3), self.amp_list):
            mrna_hyb_start_pos = amp['f3_pos'] - len(amp['p_sseq'])
            mrna_hyb_stop_pos = amp['r3_pos'] + len(amp['m_sseq'])

            product = amp['product']
            if (len(product) > max_report_length):
                continue

            if (apply_directionality_logic):
                if (amp["fwd_primer_seq"] == self.fwd_primer_seq and amp["rev_primer_seq"] == self.rev_primer_seq):
                    primer_products[ amp['real_hid'] ].append( (mrna_hyb_start_pos, mrna_hyb_stop_pos, product, True, amp) )
                elif (amp["fwd_primer_seq"] == amp["rev_primer_seq"] == self.rev_primer_seq):
                    # Enable detection of products that are formed in an RT reaction that are formed
                    # by the reverse primer acting to prime a 1st strand and then a 2nd strand synthesis.
                    primer_products[ amp['real_hid'] ].append( (mrna_hyb_start_pos, mrna_hyb_stop_pos, product, False, amp) )
                else:
                    try:
                        assert (amp["rev_primer_seq"] == self.fwd_primer_seq), "Unexpected situation in classMFE:getPrimerProducts()"
                    except AssertionError, ae:
                        pdb.set_trace()
                        print >> sys.stderr, ae.message
            else:
                if (amp["fwd_primer_seq"] == self.rev_primer_seq and amp["rev_primer_seq"] == self.fwd_primer_seq):
                    product = product[::-1].translate(DNA_complement_table)
                else:
                    assert (amp["fwd_primer_seq"] == self.fwd_primer_seq and amp["rev_primer_seq"] == self.rev_primer_seq)

                if (amp["fwd_primer_seq"] != amp["rev_primer_seq"]):
                    primer_products[ amp['real_hid'] ].append( (mrna_hyb_start_pos, mrna_hyb_stop_pos, product, True, amp) )
                else:
                    primer_products[ amp['real_hid'] ].append( (mrna_hyb_start_pos, mrna_hyb_stop_pos, product, False, amp) )

        return primer_products


    def checkForStableHeteroduplex(self, seq1, seq2, max_Tm):
        dont_exceed_max_Tm = False
        # use logic of get_align_seq(seq_list, options, product):
        return dont_exceed_max_Tm


    def verifyPrimerSpecificity(self, intended_target_IDs=set()):
        sys.exit(1) # 'mid_seq' is not the primer product, just the sequence between the primers. Use my computed 'product' instead.
        # See MFEprimer.py:tab_out() for more details
        is_okay, message = False, ""
        
        primer_product_IDs = set(map(lambda x:x[3]['real_hid'], self.amp_list))
        primer_product_seqs = set(map(lambda x:x[3]['mid_seq'], self.amp_list)) 
        primer_product_tuples = map(lambda x:(x[3]['real_hid'], x[3]['mid_seq']), self.amp_list)

        unintended_target_IDs = primer_product_IDs - intended_target_IDs 
        
        if (len(intended_target_IDs)==0):
            is_okay = len(primer_product_IDs)==0 and len(primer_product_seqs)==0
            if (not is_okay):
                message = "Product amplified, but none expected."
        else:
            is_okay = primer_product_IDs==intended_target_IDs and len(primer_product_seqs)==1
            if (is_okay):
                is_okay = all(map(lambda x: len(x[3]['m_tail'])==0 and len(x[3]['p_tail'])==0 and 
                              x[3]['pseq']==x[3]['p_qseq'] and x[3]['mseq']==x[3]['m_qseq'], self.amp_list))
            else:
                if (len(primer_product_seqs) == 0):
                    message = "No primer product. "
                elif (len(primer_product_seqs) != 1):
                    message = "Primer products of differing lengths. "
                if (primer_product_IDs != intended_target_IDs):
                    message += " Primer product IDs (%d) did not match intended targets (%d)." % (len(primer_product_IDs), len(intended_target_IDs))

        return is_okay, message, primer_product_tuples


    def processPrimer(self, fwd_primer_seq, rev_primer_seq):
        self.curr_oligos = [{'size': len(fwd_primer_seq), 'id': 'fwd', 'seq': fwd_primer_seq, 'desc': ''},
                            {'size': len(rev_primer_seq), 'id': 'rev', 'seq': rev_primer_seq, 'desc': ''}]

        self.fwd_primer_seq = fwd_primer_seq
        self.rev_primer_seq = rev_primer_seq
        
        products, binding_range = self.my_primer_process(self.options)
        
        seq_list = []
        for hid, s, e in binding_range:
            seq_list.append( self.db_two_bit[hid][s:e] )
            
        filtered_products = get_align_seq(seq_list, self.options, products)
        self.amp_list = self.my_primer_analysis(filtered_products, self.options)

        return len(self.amp_list) > 0

    
    def my_primer_process(self, options):
        oligo_pos = []
        oligo_id_list = []
        for oligo in self.curr_oligos:
            primer_seq = oligo['seq']
            oligo_id_list.append(oligo['id'])

            mer = primer_seq[-options.k_value:]
            mer_id = chilli.DNA2int(mer)

            # p for plus strand, m for minus strand
            p_pos_list, m_pos_list = self.my_get_position(options, mer_id)

            oligo_pos.append({'p_list' : p_pos_list, 'm_list' : m_pos_list})

        product = []
        binding_range = []
        binding_primer = []
        for i in xrange(len(self.curr_oligos)):
            p_list = oligo_pos[i]['p_list']
            p_oligo_length = self.curr_oligos[i]['size']
            for k in xrange(len(self.curr_oligos)):
                m_list = oligo_pos[k]['m_list']
                m_oligo_length = self.curr_oligos[k]['size']

                for j in p_list.iterkeys():
                    hid = str(j) # Because the database has been reformatted
                    try:
                        p_pos = p_list[j]
                        m_pos = m_list[j]
                    except:
                        continue   # TODO: replace expensive exception clause with simple if statement

                    for p in p_pos:
                        left = get_pos_range(p, m_pos)
                        for pos_index in xrange(left, len(m_pos)):
                            f3_pos = p+1
                            r3_pos = m_pos[pos_index]

                            # !!!!!! TODO !!!!!!
                            # Reconsider the logic here. Doesn't the condition f3_pos >= r3_pos, which I
                            # assume to be the primers' 3' end positions, natively create something indistinguishable
                            # from a primer dimer (as long as neither is beyond the primer 5' position). If so, I do
                            # not want to ignore this case, as is done here. Also, a simple heteroduplex analysis
                            # performed beforehand wouldn't flag a problem if the overlap was only a few positions.
                            # The amplicon size <= p.len + m.len
                            if f3_pos >= r3_pos:
                                continue

                            product_size = p_oligo_length + m_pos[pos_index] - p + m_oligo_length - 1

                            if product_size < options.size_start:
                                continue
                            if product_size > options.size_stop:
                                break

                            # TODO: Should these p_start & m_stop overruns raise an exception instead?
                            p_start = p - p_oligo_length + 1
                            if p_start < 0:
                                p_start = 0

                            m_stop = r3_pos + m_oligo_length
                            if m_stop > self.fcdict[hid]['size']:
                                m_stop = self.fcdict[hid]['size'] 

                            binding_range.append( (hid, p_start, p + 1) )
                            # Reverse: Correction for reverse
                            binding_range.append( (hid, r3_pos, m_stop) )

                            amp = {
                                'hid' : hid,
                                'pid' : self.curr_oligos[i]['id'],
                                'mid' : self.curr_oligos[k]['id'],
                                'plen' : p_oligo_length,
                                'mlen' : m_oligo_length,
                                'pseq' : self.curr_oligos[i]['seq'],
                                'mseq' : Seq.rev_com(self.curr_oligos[k]['seq']),
                                'size' : product_size,
                                'f3_pos' : f3_pos,
                                'r3_pos' : r3_pos,
                                }

                            product.append(amp)

        return product, binding_range


    def my_get_position(self, options, mer_id):
        '''Get position of the mer from SQLite3 database'''
        plus_pos = {}
        minus_pos = {}
        query = "select plus, minus from pos where mer_id=%s" %  mer_id
        self.cur.execute(query)
        try:
            (plus, minus) = self.cur.fetchone()
        except:
            print "Error found when retrieving position values from indexed database"
            print "Is the k-value right?"
            exit(1)

        if plus:
            plus_pos = get_pos_data(plus)
        if minus:
            minus_pos = get_pos_data(minus)

        return plus_pos, minus_pos


    def my_primer_analysis(self, product, options):
        '''Analysis the candidate forward and reverse primer and check whether they can amplify an amplicon'''
        id_and_midseq_coord_list = []
        tmp_list = []
        amp_list = []
        filter_product = []

        # Forward primer
        for i in xrange(len(product)):
            # Reverse primer
            #print i
            amp = product[i]
            hid = amp['hid']
            pid = amp['pid']
            mid = amp['mid']
            f_len = amp['plen']
            r_len = amp['mlen']
            pseq = amp['pseq']
            mseq = amp['mseq']
            size = amp['size']
            f_3_pos = amp['f3_pos'] 
            r_3_pos = amp['r3_pos']

            plen = amp['plen']
            mlen = amp['mlen']

            p_qseq = amp['p_qseq']
            p_aseq = amp['p_aseq']
            p_sseq = amp['p_sseq']
            p_tail = amp['p_tail']

            m_qseq = amp['m_qseq']
            m_aseq = amp['m_aseq']
            m_sseq = amp['m_sseq']
            m_tail = amp['m_tail']

            p_Tm = amp['p_Tm']
            p_DeltaG = amp['p_DeltaG']
            m_Tm = amp['m_Tm']
            m_DeltaG = amp['m_DeltaG']
            p_3_DeltaG = TmDeltaG.calDeltaG(p_qseq[-5:], Seq.complement(p_sseq[-5:]), mono_conc=options.mono_conc, diva_conc=options.diva_conc, dntp_conc=options.dntp_conc)
            m_3_DeltaG = TmDeltaG.calDeltaG(m_qseq[:5], Seq.complement(m_sseq[:5]), mono_conc=options.mono_conc, diva_conc=options.diva_conc, dntp_conc=options.dntp_conc)

            # Filter DeltaG
            if p_3_DeltaG < float(options.dg_start) or p_3_DeltaG > float(options.dg_stop):
                continue
            if m_3_DeltaG < float(options.dg_start) or m_3_DeltaG > float(options.dg_stop):
                continue

            ppc = cal_PPC(len(p_qseq), f_len, len(m_qseq), r_len)
            # Filter by PPC
            if ppc < options.ppc:
                continue

            id_and_midseq_coord_list.append( (hid, f_3_pos, r_3_pos) )

            ave_Tm = (p_Tm + m_Tm) / 2 # For sort
            to_be_added = (ave_Tm, ppc, p_3_DeltaG, m_3_DeltaG)
            tmp_list.append(to_be_added)
            filter_product.append(amp)

        for i in xrange(len(filter_product)):
            (ave_Tm, ppc, p_3_DeltaG, m_3_DeltaG) = tmp_list[i]
            amp = filter_product[i]
            pid = amp['pid']
            mid = amp['mid']

            hid, s, e = id_and_midseq_coord_list[i]
            midseq = self.db_two_bit[hid][s:e]

            assert (hid == amp['hid'])
            real_hid = self.fcdict[hid]['id']
            hdesc = self.fcdict[hid]['desc']
            #amp['amp_graphic'] = draw_graphical_alignment_primer(amp, self.curr_oligos, options, mid_seq)
            size = amp['size']
            amp['p_3_DeltaG'] = p_3_DeltaG
            amp['m_3_DeltaG'] = m_3_DeltaG
            amp['real_hid'] = real_hid
            amp['hdesc'] = hdesc
            amp_list.append([ave_Tm, ppc, size, amp])

            fwd_primer_seq = amp['pseq'].upper()
            rev_primer_seq = amp['mseq'][::-1].translate(DNA_complement_table).upper()

            amp['fwd_primer_seq'] = fwd_primer_seq
            amp['rev_primer_seq'] = rev_primer_seq

            amplicon = "%s%s%s" % (amp['pseq'], midseq, amp['mseq'])
            assert (len(amplicon) == size), "Constructed amplification product length differs from MFEprimer size value"
            amp['product'] = amplicon

        return amp_list


if (__name__ == "__main__"):
    #genomicMFE = MFE()

    CG = "/raid2/projects/CGDB/models/CG.fa"
    transcriptomicMFE = MFE(CG)
    options = MFEOptions()
    transcriptomicMFE.processPrimer("agcgtgcgattcaaaggatt", "gactggtgccgacgatgact", options)
    transcriptomicMFE.finalize()
    
    sys.exit(0)
