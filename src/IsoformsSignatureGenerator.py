class IsoformSignatureGenerator:
    def __init__(self, transcript_id, mRNA_sequence, otp, strand):
        self.ID = transcript_id
        self.mRNA_seq = mRNA_sequence
        self.my_parts = set(otp)
        self.coord_correspondence = None

        self.max_fusionparts_separation = 500
        self.mRNA_conv = None
        self.potential_signatures = []
        self.potential_signature_index = 0
        self.potential_signatures_lookup = {}
        
        fusion_parts = []
        allowed_part_types = None
        if (any(map(lambda x: "sj" in x[3], otp))):
            allowed_part_types = set(["sj","sj,sj","uesj","sjue","ue,sj","sj,ue"])
            fusion_parts = self.mergeNearbyParts(otp)
            for genomic_coords_part, otp_indices in fusion_parts:
                self.my_parts.add( genomic_coords_part )
        else:
            allowed_part_types = set(["ue"])

        self.mapGenomicCoordsTomRNA(otp, fusion_parts, strand)
        
        # Create a 0-based mRNA index to 1-based genomic index tuple list
        if (mRNA_sequence != None):
            self.coord_correspondence = commonCGDB.correspondmRNACoordsToGenomicCoords(mRNA_sequence, otp, strand)

        for tup in filter(lambda x: x[-1] in allowed_part_types, self.my_parts):
            new_signature = Signature(self.ID, self.mRNA_seq, otp[0][0], strand, tup, self.mRNA_conv[tup], self.coord_correspondence)
            self.potential_signatures.append( new_signature )
            self.potential_signatures_lookup[tup] = new_signature


    def myID(self):
        return self.ID


    def lookupSignature(self, target_tup):
        ret_val = None
        if (self.potential_signatures_lookup.has_key(target_tup)):
            ret_val = self.potential_signatures_lookup[target_tup]
        else:
            # For backwards compatibility reasons
            for tup, S in self.potential_signatures_lookup.items():
                if (target_tup == (tup[0], tup[1], tup[-2], tup[-1])):
                    ret_val = S
                    break

        if (ret_val == None):
            # This is necessary because the signature stuff is now distinguishing between
            # adjoining merged parts and non-adjoint merged parts with a comma.
            part_type_w_comma = target_tup[-1][0:2] + "," + target_tup[-1][2:4]
            target_tup_w_comma = (target_tup[0], target_tup[1], target_tup[2], part_type_w_comma)
            if (self.potential_signatures_lookup.has_key(target_tup_w_comma)):
                ret_val = self.potential_signatures_lookup[target_tup_w_comma]
            else:
                # For backwards compatibility reasons
                for tup, S in self.potential_signatures_lookup.items():
                    if (target_tup_w_comma == (tup[0], tup[1], tup[-2], tup[-1])):
                        ret_val = S
                        break

        return ret_val
    
    
    def signatureIndicatesMe(self, aSignature):
        return any(map(lambda x: aSignature.isSameAs(x), self.potential_signatures))


    def getAllSignatures(self):
        return self.potential_signatures

    
    def nextPotentialSignature(self):
        ret_val = None
        if (self.potential_signature_index < len(self.potential_signatures)):
            ret_val = self.potential_signatures[self.potential_signature_index]
            self.potential_signature_index += 1
        return ret_val


    def mergeNearbyParts(self, otp):
        fusion_parts = []
        # Try to merge closeby splice junctions and unique exonic parts to form "new" parts
        otp_sj_ue_indices = map(lambda y:y[0], filter(lambda x:x[1][3] in ["sj","ue"], zip(range(len(otp)),otp)))
        if (len(otp_sj_ue_indices) > 1):
            for i in xrange(len(otp_sj_ue_indices)-1):
                i1 = otp_sj_ue_indices[i]
                for j in xrange(i+1,len(otp_sj_ue_indices)):
                    i2 = otp_sj_ue_indices[j]
                    intervening_seq_len = 0
                    for index in xrange(i1+1,i2):
                        if (otp[index][3] != "sj"):
                            intervening_seq_len += (otp[index][2] - otp[index][1])
                    if (otp[i1][3] == "sj"):
                        intervening_seq_len += 1
                    if (otp[i2][3] == "sj"):
                        intervening_seq_len += 1
                    if (intervening_seq_len <= self.max_fusionparts_separation):
                        if (i2-i1 == 1):
                            merged_types = otp[i1][3] + otp[i2][3]
                        else:
                            merged_types = otp[i1][3] + "," + otp[i2][3]

                        #if (merged_types in ["sj,ue", "ue,sj"]):
                        fusion_parts.append( ((otp[i1][0], otp[i1][1], otp[i1][2], otp[i2][1], otp[i2][2], merged_types), (i1,i2)) )
                        #else:
                        #    fusion_parts.append( ((otp[i1][0], otp[i1][1], otp[i2][2], merged_types), (i1,i2)) )

        return fusion_parts


    def mapGenomicCoordsTomRNA(self, otp, fusion_parts, strand):
        """Returns a dictionary for converting between genomic coordinate based part definitions to mRNA coordinate based part definitions.
        "otp" is the ordered transcript parts definitions. "fusion_parts" described (potentially) non-contiguous parts combinations, and is a list of
        ((chromosome, part 1 start, part 2 end, merged_types), (part 1 index, part 2 index)). Maintains coordinates in BED format"""
        mRNA_conv = {}
        assert (strand in ["-","+"])
        if (strand == "-"):
            otp.reverse()

        start, end = (0, otp[0][2] - otp[0][1])
        mRNA_conv[otp[0]] = (otp[0][0], start, end, otp[0][3])
        prev_end = end
        for tup in otp[1:]:
            chrom, part_start, part_end, part_type = tup
            if (part_type == "sj"):
                new_end = prev_end
                mRNA_conv[tup] = (chrom, prev_end-1, prev_end+1, part_type)
            else:
                new_end = prev_end + part_end - part_start
                mRNA_conv[tup] = (chrom, prev_end, new_end, part_type)
            prev_end = new_end

        if (strand == "-"):
            otp.reverse()

        for genomic_coords_part, otp_indices in fusion_parts:
            if (strand == "-"):
                i2, i1 = otp_indices
            else:
                i1, i2 = otp_indices    

            if (abs(i2-i1)==1): # If the parts are adjacent...
                merged_types = otp[i1][3] + otp[i2][3]
            else:
                merged_types = otp[i1][3] + "," + otp[i2][3]

            mRNA_conv[genomic_coords_part] = (mRNA_conv[otp[i1]][0], mRNA_conv[otp[i1]][1], mRNA_conv[otp[i1]][2], mRNA_conv[otp[i2]][1], mRNA_conv[otp[i2]][2], merged_types)

        self.mRNA_conv = mRNA_conv
