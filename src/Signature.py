class Signature:
    def __init__(self, parent_transcript_id, parent_mRNA_seq, chrom, strand, signature_tuple, signature_tuple_in_mRNA_coords, coord_correspondence):
        self.parent_transcript = parent_transcript_id
        self.mRNA_seq = parent_mRNA_seq
        self.chrom = chrom
        self.strand = strand
        self.signature_tuple = signature_tuple
        self.mRNA_coords = signature_tuple_in_mRNA_coords
        self.coord_correspondence = coord_correspondence
        self.primer_targets = None  # set of (isoformID, primer product nucleotide sequence)
        self.unintended_targets = None
        
        # Primer pair data
        self.pp_stratum = -1 # This is the primer3 parameter set by which the primer pair was found (i.e. "param_set_rank")
        self.fp_name = None
        self.fp_seq = None
        self.fp_tm = None
        self.fp_gc = None
        self.rp_name = None
        self.rp_seq = None
        self.rp_tm = None
        self.rp_gc = None
        self.penalty = None

        
    def setPrimerPair(self, stratum, fp_name, forward_seq, forward_tm, forward_gc, rp_name, reverse_seq, reverse_tm, reverse_gc, penalty, targeted_isoforms, primer_targets):
        self.pp_stratum = stratum
        self.fp_name = fp_name
        self.fp_seq = forward_seq
        self.fp_tm = forward_tm
        self.fp_gc = forward_gc
        self.rp_name = rp_name 
        self.rp_seq = reverse_seq
        self.rp_tm = reverse_tm
        self.rp_gc = reverse_gc
        self.penalty = penalty
        self.targeted_isoforms = targeted_isoforms
        self.primer_targets = primer_targets
        self.num_unintended_w_target_sequence = -1
        
        if (isinstance(targeted_isoforms,set)):
            self.unintended_targets = filter(lambda x: x[0] not in targeted_isoforms, primer_targets)
        elif (isinstance(targeted_isoforms,str)):
            self.unintended_targets = filter(lambda x: x[0] != targeted_isoforms, primer_targets)
        else:
            print >> sys.stderr, "ERROR: shouldn't reach this point. Exiting."
            sys.exit(1)

        # Even though other isoforms may be amplified, when products are sequenced will the intended isoform's
        # sequence be distinguishable from the rest?
        try:
            target_sequence = filter(lambda x:x[0]==self.targeted_isoforms, self.primer_targets)[0][1]
        except IndexError:
            import pdb
            pdb.set_trace()
            
        self.num_unintended_w_target_sequence = len(filter(lambda y:y[1]==target_sequence, self.unintended_targets))


    def distinguishesIntendedTarget(self):
        assert (self.num_unintended_w_target_sequence != -1)
        return (self.num_unintended_w_target_sequence == 0)


    def numUnintendedTargets(self):
        return len(self.unintended_targets)

    
    def getPrimerPair(self):
        return (self.fp_seq, self.fp_tm, self.rp_seq, self.rp_tm, self.strand)

            
    def hasPrimerPair(self):
        return (self.fp_name != None)


    def getPenalty(self):
        return self.penalty

    
    def getAmplicon(self):
        li = int(self.fp_name.split("_")[0])
        ri = int(self.rp_name.split("_")[0])
        return self.mRNA_seq[li-1:ri]


    def getAmpliconLength(self):
        li = int(self.fp_name.split("_")[0])
        ri = int(self.rp_name.split("_")[0])
        return len(self.mRNA_seq[li-1:ri])


    def getAmpliconGenomicCoords(self, op_format, gtf_score=500):
        # fp_name and rp_name are defined in terms of 1-based coordinates into the mRNA sequence
        li = int(self.fp_name.split("_")[0]) - 1
        ri = int(self.rp_name.split("_")[0])

        lines = []
        exon_number = 1
        gene_id = self.parent_transcript.split(".")[0]
        targeted_isoforms_str = ",".join(self.targeted_isoforms) if isinstance(self.targeted_isoforms,set) else self.targeted_isoforms
        if (op_format == "gtf"):
            cc = self.coord_correspondence # A 0-based mRNA index to 1-based genomic index tuple list
            curr_faux_exon_start = cc[li][1]
            last_in_genome = cc[li][1]
            for i in xrange(li+1,ri):
                curr_in_genome = cc[i][1]
                if (abs(curr_in_genome-last_in_genome)>1 or i==ri-1): # If adjacent positions in the mRNA correspond to non-adjacent in genome or at end of the transcript...
                    tuple_descr = "%s:%d-%d_%s" % (self.signature_tuple[0], self.signature_tuple[1], self.signature_tuple[-2], self.signature_tuple[-1])
                    if (self.strand == "+"):
                        lines.append( "%s\tAmplicon\texon\t%d\t%d\t%d\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number %d; signature_description \"%s\"" % \
                                      (self.chrom, curr_faux_exon_start, last_in_genome, gtf_score, self.strand, gene_id, targeted_isoforms_str, exon_number, tuple_descr) )
                    else:
                        lines.append( "%s\tAmplicon\texon\t%d\t%d\t%d\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number %d; signature_description \"%s\"" % \
                                      (self.chrom, last_in_genome, curr_faux_exon_start, gtf_score, self.strand, gene_id, targeted_isoforms_str, exon_number, tuple_descr) )
                    curr_faux_exon_start = curr_in_genome
                    exon_number += 1
                last_in_genome = curr_in_genome
                    
        else:
            print >> sys.stderr, "ERROR: have not yet handled case of %s output format for amplicons" % op_format

        return lines


    def getPrimerPairStratum(self):
        assert (self.pp_stratum > -1)
        return self.pp_stratum


    def getParentTranscript(self):
        return self.parent_transcript
    

    def getSignatureTuple(self):
        return self.signature_tuple


    def getSignatureTupleInmRNACoords(self):
        return self.mRNA_coords


    def getmRNASignatureSequence(self):
        li, ri = self.mRNA_coords[1:3]
        return self.mRNA_seq[li:ri+1]
    

    def isSameAs(self, aSignature):
        return self.signature_tuple == aSignature.getSignatureTuple()


    def composePrimer3RegionSpec(self, allow_ue_part_type):
        """mRNA information is all in left-to-right 5'->3' orientation, regardless of DNA strand on which it resides.
        Currently hardcodes minimum/maximum primer lengths as 18/23 nucleotides and min/max amplicon length as 60/450 nucleotides."""
        assert (self.mRNA_seq != None)
        region_specs = []
        part_type = self.mRNA_coords[-1]

        if (part_type == "ue" and allow_ue_part_type):
            mRNA_part = self.mRNA_coords
            assert(len(mRNA_part) == 4)
            region_spec = Primer3RegionSpecification(part_type, mRNA_part[1]+1, -1, -1, mRNA_part[2], len(self.mRNA_seq), False, False)
            if (region_spec.isFeasible()):
                region_specs.append(region_spec)
        elif (part_type != "ue"):
            mRNA_part = self.mRNA_coords
            if (len(mRNA_part) == 4):
                part_left_index, part_right_index = mRNA_part[1]+1, mRNA_part[2]
                part_right_index_internal, part_left_index_internal = -1, -1
            else:
                part_left_index, part_right_index_internal, part_left_index_internal, part_right_index = mRNA_part[1]+1, mRNA_part[2], mRNA_part[3]+1, mRNA_part[4]
                
            # Return region specs ordered by 1) both primers overlap must overlap, 2) only forward primer overlaps,
            # 3) only reverse primer overlaps, 4) neither primer overlaps
            for fp_ol, rp_ol in [(True,True), (True,False), (False,True), (False,False)]:
                region_spec = Primer3RegionSpecification(part_type, part_left_index, part_right_index_internal, part_left_index_internal, part_right_index, len(self.mRNA_seq),
                                                         fp_ol, rp_ol)
                if (region_spec.isFeasible()):
                    region_specs.append(region_spec)
        else:
            print >> sys.stderr, "WARNING: cannot handle part type \'%s\'. Skipping this signature." % self.signature_tuple[3]

        return region_specs


    def taggedPrimersHairpinCheck(self, param_set, fwd_seq, rev_seq, fwd_Tm, rev_Tm, fwd_tag, rev_tag):
        """Returns True if the (approximated) hairpin propensity of both tagged primers with random UMIs is sufficiently low.
        The Hairpin Tm is approximated as the max Tm between the tag and the primer sequence, and the UMI is excluded.
        Done this way since complete tag+umi+primer exceeds Primer3's sequence limit.
        """
        is_okay_count = 0
        
        common_primer3_input_string = "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/usr/local/src/primer3-2.3.5/src/primer3_config/\n"
        common_primer3_input_string += "PRIMER_TASK=check_primers\nPRIMER_THERMODYNAMIC_ALIGNMENT=1\nPRIMER_EXPLAIN_FLAG=1\n"
        common_primer3_input_string += "PRIMER_SALT_DIVALENT=%f\nPRIMER_DNTP_CONC=%f\n" % (param_set["PRIMER_SALT_DIVALENT"], param_set["PRIMER_DNTP_CONC"])
        common_primer3_input_string += "PRIMER_MIN_TM=40\nPRIMER_MAX_TM=92\nPRIMER_OPT_TM=75\n" # Set at these values to allow the hairpin calculations to go through

        for (seq, Tm, tag) in [(fwd_seq, fwd_Tm, fwd_tag), (rev_seq, rev_Tm, rev_tag)]:
            max_hairpin_Th = Tm - 10.0

            primer3_input_string = common_primer3_input_string
            primer3_input_string += "SEQUENCE_PRIMER=%s\nSEQUENCE_PRIMER_REVCOMP=%s\n" % (seq, tag)
            primer3_input_string += "PRIMER_PAIR_MAX_COMPL_ANY_TH=%f\n" %  max_hairpin_Th
            primer3_input_string += "=\n"
                
            primer3 = Popen("/usr/local/src/primer3-2.3.5/src/primer3_core", shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            stdout_results, stderr_results = primer3.communicate(input=primer3_input_string)
            if (primer3.returncode != 0):
                print >> sys.stderr, "WARNING: primer3_core return code was %d during hairpin checking." % primer3.returncode

            results_dict = dict( map(lambda x: x.strip().split("="), stdout_results.split("\n")[0:-1]) )                

            try:
                if (results_dict["PRIMER_PAIR_EXPLAIN"] == "considered 1, ok 1"):
                    is_okay_count += 1
                elif (results_dict["PRIMER_PAIR_EXPLAIN"] != "considered 1, high any compl 1, ok 1"):
                    print >> sys.stderr, "ERROR: unexpected EXPLAIN in hairpin checking -> %s" % results_dict["PRIMER_PAIR_EXPLAIN"]
            except KeyError:
                import pdb
                pdb.set_trace()
                #sys.exit(1)

            if (is_okay_count == 0):
                break

        return (is_okay_count == 2)


    def checkDNAPolEfficiency(self, fwd_seq, rev_seq, product, min_GandC_count, use_ideal):
        pol_efficiency_okay = 0
        fwd_runway = product[len(fwd_seq):len(fwd_seq)+2]
        rev_runway_rc = product[-len(rev_seq)-2:-len(rev_seq)]

        for seq, runway in [(fwd_seq, fwd_runway), (rev_seq,rev_runway_rc)]:
            GandC_count = seq[-6:].count("C") + seq[-6:].count("G")
            runway_count = runway.count("C") + runway.count("G")

            if (GandC_count >= min_GandC_count and seq[-3:] not in ["ACG", "AGC"]): # and runway_count > 0
                if (use_ideal):
                    if (seq[-6:-4] in ["GC", "CG"]):
                        pol_efficiency_okay += 1
                else:
                    pol_efficiency_okay += 1

        return (pol_efficiency_okay==2)


    def setPrimer3Args(self, region_specification):
        # Set the primer3-py Global arguments
        # "PRIMER_FIRST_BASE_INDEX" : 1
        # "PRIMER_PAIR_MAX_DIFF_TM" : 3
        primer3_global_args = {"PRIMER_TASK" : "generic",
                               "PRIMER_PICK_LEFT_PRIMER" : 1,
                               "PRIMER_PICK_INTERNAL_OLIGO" : 0,
                               "PRIMER_PICK_RIGHT_PRIMER" : 1,
                               "PRIMER_OPT_SIZE" : designParams.opt_primer_len,
                               "PRIMER_MIN_SIZE" : designParams.min_primer_len,
                               "PRIMER_MAX_SIZE" : designParams.max_primer_len,
                               "PRIMER_EXPLAIN_FLAG" : 1,
                               "PRIMER_NUM_RETURN" : 300,
                               "PRIMER_SALT_MONOVALENT" : designParams.monovalent_salt_conc,
                               "PRIMER_SALT_DIVALENT" : designParams.divalent_salt_conc,
                               "PRIMER_DNTP_CONC" : designParams.dntp_conc,
                               "PRIMER_DNA_CONC" : designParams.input_primer_conc_nM,
                               "PRIMER_MIN_TM" : designParams.min_primer_Tm,
                               "PRIMER_MAX_TM" : designParams.max_primer_Tm,
                               "PRIMER_OPT_TM" : designParams.opt_primer_Tm,
                               "PRIMER_TM_FORMULA" : 1,
                               "PRIMER_SALT_CORRECTIONS" : 1,
                               "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT" : 1,
                               "PRIMER_PRODUCT_OPT_SIZE" : designParams.amplicon_opt_len,
                               "PRIMER_PRODUCT_SIZE_RANGE" : (designParams.amplicon_min_len, designParams.amplicon_max_len)}
                        
        if (len(self.nuc_seq) >= 10000):
            primer3_global_args["PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT"] = 0
            primer3_global_args["PRIMER_MAX_TEMPLATE_MISPRIMING"] = 12
        else:
            primer3_global_args["PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT"] = 1
            primer3_global_args["PRIMER_MAX_TEMPLATE_MISPRIMING_TH"] = 50

        template_seq = self.mRNA_seq.upper()
        template_seq = template_seq.replace('Y','T')

        primer3_seq_args = {"SEQUENCE_ID" : self.parent_transcript,
                            "SEQUENCE_TEMPLATE" : template_seq}

        return (primer3_global_args, primer3_seq_args)


    def getCandidatePrimers(self, oligo_thermo, region_specification, fwd_tag, rev_tag, verbose=False):
    """Runs Primer3"""
        min_num_unique_primers = 25
        candidate_pairs = []
        primer_pairs = []

        assert ((fwd_tag != '' and rev_tag != '') or (fwd_tag=='', rev_tag==''))
        primer3_global_args, primer3_seq_args = self.setPrimer3Args(region_specification)
        design_results = designPrimers(primer3_seq_args, primer3_global_args)

        if (design_results.has_key("PRIMER_ERROR")):
            print >> sys.stderr, "Primer3 Error: %s" % design_results["PRIMER_ERROR"]
            pdb.set_trace()
        elif (design_results.has_key("PRIMER_WARNING")):
            print >> sys.stderr, "Primer3 Warning: %s" % design_results["PRIMER_WARNING"]
            pdb.set_trace()
        else:
            num_pps_returned = design_results["PRIMER_PAIR_NUM_RETURNED"]
            for i in xrange(num_pps_returned):
                fwd_primer_5p_pos, fwd_primer_len = design_results["PRIMER_LEFT_%d" % i]
                rev_primer_5p_pos, rev_primer_len = design_results["PRIMER_RIGHT_%d" % i]

                if (not fwd_primer_cache.has_key( (fwd_primer_5p_pos, fwd_primer_len) ) ):
                    fwd_seq = design_results["PRIMER_LEFT_%d_SEQUENCE" % i]
                    fwd_Tm, fwd_frac_duplexed, fwd_has_hairpin, fwd_penalty = oligo_thermo.calcPrimerThermoDetails("fwd", fwd_seq)
                    if (fwd_has_hairpin or fwd_penalty > 0.01):
                        fwd_primer_cache[(fwd_primer_5p_pos, fwd_primer_len)] = None
                        continue
                    else:
                        nuc_imbalance = abs(3 - fwd_seq[-6:].count('C') - fwd_seq[-6:].count('G'))
                        AT_3p = 1 if (fwd_seq[-1] == 'A' or fwd_seq[-1] == 'T') else 0
                        aux_3p_penalty = nuc_imbalance + AT_3p
                        fwd_primer = PrimerSingleton(fwd_primer_5p_pos, fwd_seq, fwd_primer_len, fwd_Tm, fwd_frac_duplexed, fwd_penalty, aux_3p_penalty)
                        fwd_primer_cache[(fwd_primer_5p_pos, fwd_primer_len)] = fwd_primer

                if (not rev_primer_cache.has_key( (rev_primer_5p_pos, rev_primer_len) ) ):
                    rev_seq = design_results["PRIMER_RIGHT_%d_SEQUENCE" % i] 
                    rev_Tm, rev_frac_duplexed, rev_has_hairpin, rev_penalty = oligo_thermo.calcPrimerThermoDetails("rev", rev_seq)
                    if (rev_has_hairpin or rev_penalty > 0.01):
                        rev_primer_cache[(rev_primer_5p_pos, rev_primer_len)] = None
                        continue
                    else:
                        nuc_imbalance = abs(3 - rev_seq[-6:].count('C') - rev_seq[-6:].count('G'))
                        AT_3p = 1 if (rev_seq[-1] == 'A' or rev_seq[-1] == 'T') else 0
                        aux_3p_penalty = nuc_imbalance + AT_3p
                        rev_primer = PrimerSingleton(rev_primer_5p_pos, rev_seq, rev_primer_len, rev_Tm, rev_frac_duplexed, rev_penalty, aux_3p_penalty)
                        rev_primer_cache[(rev_primer_5p_pos, rev_primer_len)] = rev_primer

                fwd_primer = fwd_primer_cache[(fwd_primer_5p_pos, fwd_primer_len)]
                rev_primer = rev_primer_cache[(rev_primer_5p_pos, rev_primer_len)]

                if (fwd_primer != None and rev_primer != None):
                    sum_individ_penalties = fwd_primer.thermo_penalty + rev_primer.thermo_penalty
                    sum_aux_penalties = fwd_primer.aux_3p_penalty + rev_primer.aux_3p_penalty
                    rounded_amplicon_len = round(rev_primer.template_5p_pos - fwd_primer.template_5p_pos + 1, -1) # round to nearest tens
                    candidate_pairs.append( (fwd_primer, rev_primer, sum_individ_penalties, round(sum_individ_penalties,3), sum_aux_penalties, rounded_amplicon_len) )


        # Select primer pairs with a diversity of fwd and rev primer sequences
        candidate_pairs.sort(key=itemgetter(3,4,5))
        unique_fwd_seq = set()
        unique_rev_seq = set()
        i = 0
        while (i < len(candidate_pairs) and len(unique_fwd_seq) < min_num_unique_primers and len(unique_rev_seq) < min_num_unique_primers):
            fwd_primer, rev_primer, sum_individ_penalties, rounded_sum_individ_penalties, sum_aux_penalties, amplicon_len = candidate_pairs[i]
            pair_penalty_term = oligo_thermo.calcPrimerPairThermoDetails(fwd_primer[1], rev_primer[1])
            primer_pair_penalty = sum_individ_penalties + pair_penalty_term
            if (primer_pair_penalty < 0.01):
                primer_pairs.append( (fwd_primer, rev_primer, round(primer_pair_penalty,3), sum_aux_penalties, amplicon_len) )
                unique_fwd_seq.add(fwd_primer[1])
                unique_rev_seq.add(rev_primer[1])
            i += 1

        primer_pairs.sort(key=itemgetter(2,3,4))
        print >> sys.stderr, "Selected %d candidate pairs" % len(primer_pairs)
        return primer_pairs


    def savitzky_golay(self, y, window_size, order, deriv=0, rate=1):
        r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
        The Savitzky-Golay filter removes high frequency noise from data.
        It has the advantage of preserving the original shape and
        features of the signal better than other types of filtering
        approaches, such as moving averages techniques.
        Parameters
        ----------
        y : array_like, shape (N,)
        the values of the time history of the signal.
        window_size : int
        the length of the window. Must be an odd integer number.
        order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
        deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
        Returns
        -------
        ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
        Notes
        -----
        The Savitzky-Golay is a type of low-pass filter, particularly
        suited for smoothing noisy data. The main idea behind this
        approach is to make for each point a least-square fit with a
        polynomial of high order over a odd-sized window centered at
        the point.
        Examples
        --------
        t = np.linspace(-4, 4, 500)
        y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
        ysg = savitzky_golay(y, window_size=31, order=4)
        import matplotlib.pyplot as plt
        plt.plot(t, y, label='Noisy signal')
        plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
        plt.plot(t, ysg, 'r', label='Filtered signal')
        plt.legend()
        plt.show()
        References
        ----------
        .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
        Data by Simplified Least Squares Procedures. Analytical
        Chemistry, 1964, 36 (8), pp 1627-1639.
        .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
        W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
        Cambridge University Press ISBN-13: 9780521880688
        """
        try:
            window_size = np.abs(np.int(window_size))
            order = np.abs(np.int(order))
        except ValueError, msg:
            raise ValueError("window_size and order have to be of type int")
        if window_size % 2 != 1 or window_size < 1:
            raise TypeError("window_size size must be a positive odd number")
        if window_size < order + 2:
            raise TypeError("window_size is too small for the polynomials order")
        order_range = range(order+1)
        half_window = (window_size -1) // 2
        # precompute coefficients
        b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
        # pad the signal at the extremes with
        # values taken from the signal itself
        firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
        lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
        y = np.concatenate((firstvals, y, lastvals))
        return np.convolve( m[::-1], y, mode='valid')


    def checkAmpliconMeltCurve(self, amplicon_seq, melt_curve_failures_cache):
        melt_curve_shows_one_product = False
        cached_result_used = False
        
        # If a sequence very similar to 'amplicon_seq' has been previously evaluated, use results from that previous evaluation
        containing_seqs = filter(lambda x: amplicon_seq in x, melt_curve_failures_cache)
        if (len(containing_seqs)>0):
            # How much longer any one previous failed sequence is on one end of the amplicon_seq match
            the_min_of_the_maxes_for_containing_seq = min(map(lambda w: max(map(len,w)), map(lambda y:y.split(amplicon_seq), containing_seqs)))
            if (the_min_of_the_maxes_for_containing_seq < 10):
                cached_result_used = True
                melt_curve_failures_cache.add(amplicon_seq)
                
        contained_seqs = filter(lambda x: x in amplicon_seq, melt_curve_failures_cache)    
        if (not cached_result_used and len(contained_seqs)>0):
            # How much longer the amplicon_seq is on one end for each match to one of the contained, previously-failed sequences
            the_min_of_the_maxes_for_contained_seq = min(map(lambda w: max(map(len,w)), map(lambda y: amplicon_seq.split(y), contained_seqs)))
            if (the_min_of_the_maxes_for_contained_seq < 10):
                cached_result_used = True
                melt_curve_failures_cache.add(amplicon_seq)

        if (not cached_result_used):
            base_url = "https://www.dna.utah.edu/db/services/cgi-bin/udesign_ucsd.cgi"
            request_data = {'seq':amplicon_seq, 'rs':0, 'cation':20, 'mg':2, 'dmso':0}
            retry_count = 0
            while (retry_count < 10):
                try:
                    response = requests.get(base_url, params=request_data, timeout=30)
                    retry_count = 100
                except requests.exceptions.Timeout:
                    print >> sys.stderr, "ERROR: timeout (10 sec) exception URL https://www.dna.utah.edu. Retry count %d" % retry_count
                    time.sleep(5)
                    retry_count += 1
                except requests.exceptions.RequestException as e:
                    print >> sys.stderr, "ERROR: general exception URL https://www.dna.utah.edu. Retry count %d" % retry_count
                    print >> sys.stderr, e
                    time.sleep(5)
                    retry_count += 1

            if (retry_count == 10):
                print >> sys.stderr, "ERROR: giving up for sequence %s" % amplicon_seq
            else:
                doc = lxml.etree.fromstring( str(response.text) )
                D = dict(((elt.tag,elt.text) for elt in doc.getchildren()[0]))
                T = map(float, D['temperature'].split())
                H = np.array(map(float, D['helicity'].split()), 'f')

                neg_dH = -1.0 * self.savitzky_golay(y=H, window_size=7, order=5, deriv=1)

                # Determine approximate temperature difference between adjacent melting datapoints
                T_diffs = np.array(T,'f')[1:] - np.array(T,'f')[:-1]
                max_T_diff = max(map(lambda y:y[1], filter(lambda x:x[0]>67.0 and x[0]<80.0, zip(list(T[1:]), list(T_diffs)))))
                adj_T_cutoff = 1.1 * max_T_diff

                # For the noise cutoff, find the 10 consecutive values with the lowest median value
                Tm = T[np.argmax(neg_dH)]
                T_and_deriv = filter(lambda x:x[0]>=67.0 and x[0]<=min(82.0,Tm-3*max_T_diff), zip(T,neg_dH))
                if (len(T_and_deriv)<=10):
                    T_and_deriv    = filter(lambda x:x[0]>=67.0, zip(T,neg_Dh))[0:11]
                lowest_median, stddev_for_lowest_median = 1e10, None
                for i in xrange(len(T_and_deriv)-10):
                    region_median = np.median(map(lambda y:y[1], T_and_deriv[i:i+10]))
                    if (region_median < lowest_median):
                        lowest_median = region_median
                        stddev_for_lowest_median = np.std(map(lambda y:y[1], T_and_deriv[i:i+10]))

                baseline = lowest_median + 10.0*stddev_for_lowest_median
                neg_dH_noise_region_cutoff = baseline + 0.05 * (max(neg_dH) - baseline)
                lenient_neg_dH_noise_region_cutoff = 0.15 * max(neg_dH) # Allow small bumps. This cutoff really overules the other, making the other superfluous.
                
                neg_dH_local_maxima = []
                for i in xrange(1, len(neg_dH)-1):
                    if (neg_dH[i] > neg_dH_noise_region_cutoff and neg_dH[i] > lenient_neg_dH_noise_region_cutoff and neg_dH[i]>neg_dH[i-1] and neg_dH[i]>neg_dH[i+1]):
                        neg_dH_local_maxima.append( (T[i],neg_dH[i]) )

                melt_curve_shows_one_product = len(neg_dH_local_maxima)==1
                if (not melt_curve_shows_one_product):
                    melt_curve_failures_cache.add(amplicon_seq)

                if (False and not melt_curve_shows_one_product):
                    print >> sys.stderr, "WARNING: melt curve indicates multiple products"
                    print >> sys.stderr, amplicon_seq
                    plt.plot(T, neg_dH, 'k')
                    plt.plot([65,100], [neg_dH_noise_region_cutoff, neg_dH_noise_region_cutoff], 'r--', lw=2)
                    plt.title("%s" % self.parent_transcript)
                    plt.savefig("%s_two_melt_peaks.jpg" % self.parent_transcript)
                    plt.show()
                    import pdb
                    pdb.set_trace()
            
        return melt_curve_shows_one_product


    # TODO: modify so that only returns a primer pair with better penalty and shorter amplicon, OR if yields fewer unwanted amplicons
    def findPrimerPair(self, transcriptomic_MFE, genomic_MFE, mfe_options, target_isoform_ID, oligo_thermo,
                       fwd_tag='', rev_tag='', allow_unintended_targets=False, allow_ue_part_type=False, verbose=False, skip_melt_curve=True):
            
        assert (self.mRNA_seq != None)
        part_type = self.mRNA_coords[-1]
        message, genome_message, transcriptome_message = "", "", ""
        found_primer_pair = 0
        # Extract target/exclude/include region based on part type
        region_specifications = self.composePrimer3RegionSpec(allow_ue_part_type)

        melt_curve_failures_cache = set()
        while (not self.hasPrimerPair()):
            all_candidate_primer_pairs = []
            for region_specification in region_specifications:
                candidate_primer_pairs = self.getCandidatePrimers(oligo_thermo, region_specification, fwd_tag, rev_tag)
                all_candidate_primer_pairs.extend(candidate_primer_pairs)

            if (len(all_candidate_primer_pairs) > 0):
                for pp in ranked_primer_pairs:
                    if (not self.hasPrimerPair()):
                        # Check for unintended amplifications from the transcriptome and from the genome reference sequence.
                        ampl_success = transcriptomic_MFE.processPrimer(pp.fwd_seq, pp.rev_seq, mfe_options)
                        assert (ampl_success), "Primers unexpectedly did not amplify"
                        transcriptome_specificity_okay, transcriptome_message, transcriptome_primer_targets = transcriptomic_MFE.verifyPrimerSpecificity(set([target_isoform_ID]))

                        if (not transcriptome_specificity_okay and verbose):
                            print >> sys.stderr, "Transcriptome specificity failed: %s" % transcriptome_message
                        elif (transcriptome_specificity_okay or allow_unintended_targets):
                            if (part_type != "ue" or allow_ue_part_type):
                                genomic_MFE.processPrimer(pp.fwd_seq, pp.rev_seq, mfe_options)
                                #print >> sys.stderr, "Genome Hit:", genomic_MFE.getTargetNames()                                
                                genome_specificity_okay, genome_message, genome_primer_targets = genomic_MFE.verifyPrimerSpecificity()
                                if (not genome_specificity_okay and verbose):
                                    print >> sys.stderr, "Genome specificity failed: %s" % genome_message
                                elif (genome_specificity_okay):
                                    if (not skip_melt_curve):
                                        product_melt_curve_okay = self.checkAmpliconMeltCurve(pp.product, melt_curve_failures_cache)
                                    if (skip_melt_curve or product_melt_curve_okay):
                                        found_primer_pair = pp.product_size
                                        fp_name = "%d_fp" % pp.fwd_start
                                        rp_name = "%d_rp" % pp.rev_start
                                        self.setPrimerPair(param_set_rank, fp_name, pp.fwd_seq, pp.fwd_Tm, pp.fwd_gc_perc,
                                                           rp_name, pp.rev_seq, pp.rev_Tm, pp.rev_gc_perc, pp.primer_pair_penalty, target_isoform_ID,
                                                           transcriptome_primer_targets)

                                    elif (verbose):
                                        print >> sys.stderr, "INFO: failed melt curve"
                            elif (allow_ue_part_type==True):
                                if (not skip_melt_curve):
                                    product_melt_curve_okay = self.checkAmpliconMeltCurve(pp.product, melt_curve_failures_cache)
                                if (skip_melt_curve or product_melt_curve_okay):
                                    found_primer_pair = pp.product_size
                                    fp_name = "%d_fp" % pp.fwd_start
                                    rp_name = "%d_rp" % pp.rev_start
                                    self.setPrimerPair(param_set_rank, fp_name, pp.fwd_seq, pp.fwd_Tm, pp.fwd_gc_perc,
                                                       rp_name, pp.rev_seq, pp.rev_Tm, pp.rev_gc_perc, pp.primer_pair_penalty, target_isoform_ID, transcriptome_primer_targets)


        if (found_primer_pair > 0):
            message = "INFO: Found appropriate primer pair with %dbp amplicon and %d unintended targets. Num with same target seq = %d" % \
              (found_primer_pair, len(self.unintended_targets), self.num_unintended_w_target_sequence)
        elif (transcriptome_message != ""):
            message = transcriptome_message
        elif (genome_message != ""):
            message = genome_message

        return message


    def isBetterThan(self, other_primered_signature):
        '''Is better if:
              1) Lower number of unintended target sequences that are the same as the intended target's sequence, while maintaining similar Primer3 penalty
              2) Fewer unintended targets, while maintaining similar Primer3 penalty and same number of same-as-target sequences
        '''
        this_better_than_other, reason = False, ""

        penalties_similar = other_primered_signature.penalty <= 2.0 #1.2 * self.penalty
        this_fewer_unintended = len(self.unintended_targets) < len(other_primered_signature.unintended_targets)
        target_sequence = filter(lambda x:x[0]==self.targeted_isoforms, self.primer_targets)[0][1]

        num_this_equal_to_intended = len(filter(lambda y:y[1]==target_sequence, self.unintended_targets))
        num_other_equal_to_intended = len(filter(lambda y:y[1]==target_sequence, other_primered_signature.unintended_targets))
        this_fewer_to_intended = num_this_equal_to_intended < num_other_equal_to_intended
        this_same_to_intended = num_this_equal_to_intended == num_other_equal_to_intended

        if (this_fewer_to_intended and penalties_similar):
            this_better_than_other = True
            reason = "New best - fewer (%d from %d) unintended isoforms with identical target sequence and with similar Primer3 penalty." % \
              (num_this_equal_to_intended, num_other_equal_to_intended)
        elif (this_fewer_unintended and penalties_similar and this_same_to_intended):
            assert (num_this_equal_to_intended == num_other_equal_to_intended)
            this_better_than_other = True
            reason = "New best - fewer unintended isoforms (%d from %d), but same number same-as-target sequences (%d) and with similar Primer3 penalty." % \
              (len(self.unintended_targets), len(other_primered_signature.unintended_targets), num_this_equal_to_intended)

        return (this_better_than_other, reason)


