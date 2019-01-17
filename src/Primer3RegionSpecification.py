class Primer3RegionSpecification:
    """Input coordinate positions must be 1-based and inclusive. Region specs produced are in 1-based inclusive coordinates"""
    def __init__(self, part_type, part_left_index, part_right_index_internal, part_left_index_internal, part_right_index, 
                 mRNA_seq_len, fp_overlap, rp_overlap, min_primer_len=18, max_primer_len=23, min_amplicon_len=60, max_amplicon_len=450):
        self.part_type = part_type
        self.pli = part_left_index
        self.part_right_index_internal = part_right_index_internal
        self.part_left_index_internal = part_left_index_internal
        self.pri = part_right_index
        self.seq_len = mRNA_seq_len
        self.fp_overlap = fp_overlap
        self.rp_overlap = rp_overlap
        self.m5po = 7 # PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION parameter for primer3_core
        self.m3po = 4 # PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION parameter for primer3_core
        self.min_primer_len = min_primer_len
        self.max_primer_len = max_primer_len
        self.min_amplicon_len = min_amplicon_len
        self.max_amplicon_len = max_amplicon_len

        if (self.part_type in ["sj,ue", "ue,sj"]):
            assert (self.part_right_index_internal >=1 and self.part_left_index_internal >= 1)

        self.region_spec = None
        self.FP_check_fxn = None
        self.RP_check_fxn = None
        
        self.composeSpecification()

        
    def isFeasible(self):
        return self.region_spec != None


    def getRegionSpecificationPrimer3Format(self):
        if (self.region_spec == None):
            print >> sys.stderr, "ERROR attempting to get a non-existent region specification. Exiting."
            sys.exit(1)
        return self.region_spec
    

    def primersOkay(self, fp, rp):
        """Each of fp and rp are a tuple of (primer start position, primer end position) in 1-based coordinates.
        The primer pair start and end are in the coordinates of the relevant strand, so while
        forward primer start < forward primer end, reverse primer start > reverse primer end."""
        return self.FP_check_fxn(fp[0],fp[1]) and self.RP_check_fxn(rp[0],rp[1])
    
        
    def composeSpecification(self):
        seq_len, part_type, li, ri, m5po, m3po = self.seq_len, self.part_type, self.pli, self.pri, self.m5po, self.m3po
        P_min, P_max, A_min, A_max  = self.min_primer_len, self.max_primer_len, self.min_amplicon_len, self.max_amplicon_len
        fpr_li, fpr_ri, rpr_li, rpr_ri = None, None, None, None
        sojl = None
        
        self.part_right_index_internal
        self.part_left_index_internal
        uniqueness_constraint = ""
        unique_primer_region_len = None
        
        if (part_type == "ue" and not self.fp_overlap and not self.rp_overlap):
            self.region_spec = "PRIMER_MIN_THREE_PRIME_DISTANCE=2\nPRIMER_NUM_RETURN=500\n"
            self.FP_check_fxn = lambda ps,pe: True
            self.RP_check_fxn = lambda ps,pe: True            
        elif (self.fp_overlap and self.rp_overlap):
            if (part_type == "sjue"):
                fpr_li = max(1, li - P_max + m3po + 1)
                fpr_ri = li + P_max - m5po
                rpr_li = li + m3po + 1
                rpr_ri = min(seq_len, ri + P_max - m3po)
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    li - m5po + 1 >= fpr_li and li + m3po + P_min <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = " ".join( map(str, range(rpr_li + m3po - 1, rpr_ri - m5po + 1)) )
                    uniqueness_constraint = "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=0\n"
                    unique_primer_region_len = fpr_ri - fpr_li + 1
                    self.FP_check_fxn = lambda ps,pe: ps <= li - m5po + 1 and pe >= li + m3po
                    self.RP_check_fxn = lambda ps,pe: True
            elif (part_type == "sj,ue"):
                fpr_li = max(1, li - P_max + m3po + 1)
                fpr_ri = li + P_max - m5po
                rpr_li = max(self.part_left_index_internal - P_max + m5po, li + m3po + 1)
                rpr_ri = min(seq_len, ri + P_max - m3po)
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and 
                    li - m5po + 1 >= fpr_li and self.part_left_index_internal + m5po - 1 <= rpr_ri and
                    fpr_ri - fpr_li + 1 >= P_min):
                    sojl = " ".join(map(str, range(rpr_li + m3po - 1, rpr_ri - m5po + 1)))
                    uniqueness_constraint = "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=0\n"
                    unique_primer_region_len = fpr_ri - fpr_li + 1
                    self.FP_check_fxn = lambda ps,pe: ps <= li - m5po + 1 and pe >= li + m3po
                    self.RP_check_fxn = lambda ps,pe: True
            elif (part_type == "uesj"):
                fpr_li = max(1, li - P_max + m3po)
                fpr_ri = ri - m3po - 1
                rpr_li = ri - P_max + m5po
                rpr_ri = min(seq_len, ri + P_max - m3po - 1)
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and 
                    ri - m3po - P_min >= fpr_li and ri + m5po - 1 <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = " ".join(map(str, range(fpr_li + m5po - 1, fpr_ri - m3po + 1)))
                    uniqueness_constraint = "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=0\n"
                    unique_primer_region_len = rpr_ri - rpr_li + 1
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: ps >= ri + m5po - 1 and pe <= ri - m3po
            elif (part_type == "ue,sj"):
                fpr_li = max(1, li - P_max + m3po)
                fpr_ri = min(self.part_right_index_internal - P_max + m3po + 1, ri - P_max + m5po - 1)
                rpr_li = ri - P_max + m5po
                rpr_ri = min(seq_len, ri + P_max - m3po - 1)
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    self.part_right_index_internal - m5po + 1 >= fpr_li and ri + m5po - 1 <= rpr_ri and
                    fpr_ri - fpr_li + 1 >= P_min):
                    sojl = " ".join(map(str, range(fpr_li + m5po - 1, fpr_ri - m3po + 1)))
                    uniqueness_constraint = "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=0\n"
                    unique_primer_region_len = rpr_ri - rpr_li + 1
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: ps >= ri + m5po - 1 and pe <= ri - m3po
            elif (part_type == "sj,sj"):
                fpr_li = max(1, li - P_max + m3po + 1)
                fpr_ri = li + P_max - m5po
                rpr_li = ri - P_max + m5po
                rpr_ri = min(seq_len, ri + P_max - m3po - 1)
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    li - m5po + 1 >= fpr_li and ri + m5po - 1 <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = str(ri-1)
                    uniqueness_constraint = "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=0\n"
                    unique_primer_region_len = fpr_ri - fpr_li + 1
                    self.FP_check_fxn = lambda ps,pe: ps <= li - m5po + 1 and pe >= li + m3po
                    self.RP_check_fxn = lambda ps,pe: True                
            elif (part_type != "sj"):
                print >> sys.stderr, "ERROR: unrecognized part type %s. Exiting." % part_type
                sys.exit(1)
                
        elif (self.fp_overlap and not self.rp_overlap):
            if (part_type == "sj"):
                fpr_li = max(1, li - P_max + 1 + m3po)
                fpr_ri = li + P_max - m5po
                rpr_li = ri + m3po
                rpr_ri = min(seq_len, ri - m5po + A_max - 1)    
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    li - m5po + 1 >= fpr_li and    li + m3po + P_min <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = str(li)
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: True
            elif (part_type == "sjue"):
                fpr_li = max(1, li - P_max + m3po + 1)
                fpr_ri = li + P_max - m5po
                rpr_li = ri + 1
                rpr_ri = min(seq_len, li - m5po + A_max)
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    li - m5po + 1 >= fpr_li and ri + P_min <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = str(li)
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: True
            elif (part_type == "sj,ue"):
                fpr_li = max(1, li - P_max + m3po + 1)
                fpr_ri = li + P_max - m5po
                rpr_li = ri + 1
                rpr_ri = min(seq_len, li - m5po + A_max)        
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and 
                    li - m5po + 1 >= fpr_li and ri + P_min <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = str(li)
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: True                
            elif (part_type == "uesj"):
                fpr_li = max(1, li - P_max + m3po)
                fpr_ri = ri + m3po - 1
                rpr_li = ri + 1
                rpr_ri = min(seq_len, ri - m5po + A_max - 1)
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and 
                    ri - m5po >= fpr_li and ri + P_min <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = " ".join(map(str, range(fpr_li + m5po - 1, ri)))
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: True
            elif (part_type == "ue,sj"):
                fpr_li = max(1, li - P_max + m3po)
                fpr_ri = self.part_right_index_internal - P_max + m3po + 1
                rpr_li = ri + 1
                rpr_ri = min(seq_len, self.part_right_index_internal - m5po + A_max)        
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    self.part_right_index_internal - m5po + 1 >= fpr_li and ri + P_min <= rpr_ri and
                    fpr_ri - fpr_li + 1 >= P_min):
                    sojl = " ".join(map(str, range(fpr_li + m5po - 1, self.part_right_index_internal + 1)))
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: True
            elif (part_type == "sj,sj"):
                fpr_li = max(1, li - P_max + m3po + 1)
                fpr_ri = li + P_max - m5po
                rpr_li = ri + 1
                rpr_ri = min(seq_len, li - m5po + A_max)        
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    li - m5po + 1 >= fpr_li and ri + P_min <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = str(li)
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: True
            else:
                print >> sys.stderr, "ERROR: unrecognized part type %s. Exiting." % part_type
                sys.exit(1)

        elif (not self.fp_overlap and self.rp_overlap):
            if (part_type == "sj"):
                fpr_li = max(1, li + m5po - A_max + 1)
                fpr_ri = ri - m3po - 1
                rpr_li = ri - P_max + m5po
                rpr_ri = min(seq_len, li + P_max - m3po)
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    li - m3po + 1 - P_min >= fpr_li and li + m5po <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = str(li)
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: True
            elif (part_type == "sjue"):
                fpr_li = max(1, li + m5po - A_max + 1)
                fpr_ri = li - 1
                rpr_li = li - m3po + 1
                rpr_ri = min(seq_len, ri + m5po)
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    li - P_min >= fpr_li and li + m5po <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = " ".join(map(str, range(li, rpr_ri - m5po + 1)))
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: True
            elif (part_type == "sj,ue"):
                fpr_li = max(1, self.part_left_index_internal + m5po - A_max)
                fpr_ri = li - 1
                rpr_li = self.part_left_index_internal - P_max + m5po
                rpr_ri = min(seq_len, ri + P_max - m3po)
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    li - P_min >= fpr_li and self.part_left_index_internal + m5po - 1 <= rpr_ri and
                    fpr_ri - fpr_li + 1 >= P_min):
                    sojl = " ".join(map(str, range(rpr_li + m3po - 1, rpr_ri - m5po + 1)))
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: True                
            elif (part_type == "uesj"):
                fpr_li = max(1, li - P_max + m5po - A_max + 1)
                fpr_ri = li - 1
                rpr_li = li - P_max + m5po
                rpr_ri = min(seq_len, ri + m5po - 1)
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    li - P_min >= fpr_li and ri + m5po - 1 <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = " ".join(map(str, range(rpr_li + m3po - 1, ri)))
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: True
            elif (part_type == "ue,sj"):
                fpr_li = max(1, li - P_max + m5po - A_max + 1)
                fpr_ri = li - 1
                rpr_li = ri - P_max + m5po
                rpr_ri = min(seq_len, ri + P_max - 1 - m3po)
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    li - P_min >= fpr_li and ri + m5po - 1 <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = str(ri-1)
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: True
            elif (part_type == "sj,sj"):
                fpr_li = max(1, li + m5po - A_max + 1)
                fpr_ri = li - 1
                rpr_li = ri - P_max + m5po
                rpr_ri = min(seq_len, ri + P_max - m3po - 1)
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    li - P_min >= fpr_li and ri + m5po - 1 <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = str(ri-1)
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: True
            else:
                print >> sys.stderr, "ERROR: unrecognized part type %s. Exiting." % part_type
                sys.exit(1)

        elif (not self.fp_overlap and not self.rp_overlap):
            if (part_type in ["sj", "sjue", "sj,ue", "uesj", "ue,sj", "sj,sj"]):
                fpr_li = max(1, ri + P_max - A_max + 1)
                fpr_ri = li - 1
                rpr_li = ri + 1
                rpr_ri = min(seq_len, li - P_max + A_max - 1)        
                if (rpr_ri - fpr_li + 1 - 2*P_min >= A_min and rpr_li - fpr_ri + 2*P_min <= A_max and
                    li - P_min >= fpr_li and ri + P_min <= rpr_ri and fpr_ri - fpr_li + 1 >= P_min):
                    sojl = ""
                    self.FP_check_fxn = lambda ps,pe: True
                    self.RP_check_fxn = lambda ps,pe: True
            else:
                print >> sys.stderr, "ERROR: unrecognized part type %s. Exiting." % part_type
                sys.exit(1)

        else:
            print >> sys.stderr, "ERROR: forward/reverse primer overlap specifications improperly set. Exiting."
            sys.exit(1)


        if (sojl != None): # fpr_li != None and fpr_ri != None and rpr_li != None, rpr_ri != None):
            self.region_spec = "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=%d,%d,%d,%d\n" % (fpr_li, fpr_ri-fpr_li+1, rpr_li, rpr_ri-rpr_li+1)
            if (uniqueness_constraint != ""):
                self.region_spec += uniqueness_constraint
                num_return = np.sum(map(lambda x: unique_primer_region_len - x + 1, range(self.min_primer_len, self.max_primer_len+1)))
                self.region_spec += "PRIMER_NUM_RETURN=%d\n" % num_return
            else:
                self.region_spec += "PRIMER_NUM_RETURN=10\n"
            if (sojl != ""):
                self.region_spec += "SEQUENCE_OVERLAP_JUNCTION_LIST=%s\n" % sojl

