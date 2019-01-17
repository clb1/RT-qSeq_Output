class SeqRead(object):
    '''Sequencing reads are kept or placed in their 5'->3' orientation for the strand on which they occur'''
    DNA_complement_table = maketrans("ACGTNacgtn","TGCANtgcan")

    def __init__(self, ID, nucleotide_sequence, quality_string, quality_scores, from_reverse_strand):
        self.ID = ID
        self.original_sequencing_readID = ID
        self.seq = nucleotide_sequence
        self.qual_string = quality_string
        self.qual_scores = quality_scores

        # Start and stop positions are in the context of the 5' and 3' position of the primer/sequencing adapter molecules.
        # So if a primer is sense with the sequencing read, start < stop. But if antisense, start > stop.
        self.orientation_5p = None
        self.start_5p = None
        self.stop_5p = None
        self.umi_5p = None
        self.umi_5p_qual = None

        self.orientation_3p = None
        self.start_3p = None
        self.stop_3p = None
        self.umi_3p = None
        self.umi_3p_qual = None

        self.orientation_udtd5 = None
        self.stop_udtd5 = None
        self.start_udtd5 = None

        self.orientation_udtd7 = None
        self.stop_udtd7 = None
        self.start_udtd7 = None        

        if (from_reverse_strand):
            self.reverse_complement()

            
    def resetID(self, new_ID):
        self.ID = new_ID
        

    def getID(self):
        return self.ID


    def getSeqLen(self):
        return len(self.seq)


    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result


    def reverse_complement(self):
        self.seq = translate(self.seq, self.DNA_complement_table)[::-1]
        self.qual_string = self.qual_string[::-1]
        self.qual_scores = self.qual_scores[::-1]


    def getSequence(self, as_reverse_complement=False):
        ret_seq = self.seq
        if (as_reverse_complement):
            ret_seq = translate(self.seq, self.DNA_complement_table)[::-1]
        return ret_seq


    def getQualString(self, return_reversed=False):
        ret_string = self.qual_string
        if (return_reversed):
            ret_string = self.qual_string[::-1]
        return ret_string


    def getAdapterTrimmedSequence(self):
        reject_reason = None
        
        if (self.orientation_udtd5 != None and self.orientation_udtd7 != None): # Should only contain one adapter, at most
            reject_reason = "Read has both adapters"
            ret_seq = ""
        else:
            ret_seq = self.seq

            try:
                if (self.orientation_udtd5 != None):
                    assert (self.orientation_udtd5 == "antisense")
                    ret_seq = ret_seq[0:self.stop_udtd5]
                elif (self.orientation_udtd7 != None):
                    assert (self.orientation_udtd7 == "antisense")
                    ret_seq = ret_seq[0:self.stop_udtd7]
            except AssertionError:
                reject_reason = "Read has 'sense' adapter"
                
        return ret_seq, reject_reason

    
    def asFasta(self, as_reverse_complement=False):
        assert (len(self.seq) > 0 and len(self.qual_string)>0)
        ret_string = None
        if (as_reverse_complement):
            revcomp_seq = translate(self.seq, self.DNA_complement_table)[::-1]
            ret_string = ">%s\n%s" % (self.ID, revcomp_seq)
        else:
            ret_string = ">%s\n%s" % (self.ID, self.seq[::-1])
        return ret_string

    
    def asFastq(self, as_reverse_complement=False):
        assert (len(self.seq) > 0 and len(self.qual_string)>0)
        ret_string = None
        if (as_reverse_complement):
            revcomp_seq = translate(self.seq, self.DNA_complement_table)[::-1]
            ret_string = "@%s\n%s\n+\n%s" % (self.ID, revcomp_seq, self.qual_string[::-1])
        else:
            ret_string = "@%s\n%s\n+\n%s" % (self.ID, self.seq, self.qual_string)
        return "@%s\n%s\n+\n%s" % (self.ID, self.seq, self.qual_string)

    
    def addMatch(self, match_spec, match_start_pos, match_stop_pos, start_pos_in_reversed_coords, umi_len):
        """Record the match as it occurs in the sequencing read's relative context (ie native 5'->3' orientation)."""
        assert (match_spec in ["5p","3p"] and match_spec not in ["UDTD5", "UDTD7"])
        reject_reason = None

        try:
            assert (match_start_pos >= 0 and match_stop_pos <= len(self.seq))
        except AssertionError:
            pdb.set_trace()
            
        if (start_pos_in_reversed_coords):
            orientation = "antisense"
            match_start_pos, match_stop_pos = (len(self.seq)-match_start_pos, len(self.seq)-match_stop_pos)      # TODO: didn't include the "-1" here...
            if (match_spec in ["5p","3p"]):
                umi = translate(self.seq[match_start_pos:match_start_pos+umi_len][::-1], self.DNA_complement_table)  # ...so that "+1" doesn't need to be included here
                umi_qual = self.qual_scores[match_start_pos:match_start_pos+umi_len][::-1]
        else:
            orientation = "sense"
            if (match_spec in ["5p","3p"]):
                umi = self.seq[0:match_start_pos][-umi_len:]
                umi_qual = self.qual_scores[0:match_start_pos][-umi_len:]

        if (match_spec == "5p"):
            assert (self.start_5p == None)
            self.orientation_5p = orientation
            self.start_5p = match_start_pos
            self.stop_5p = match_stop_pos
            self.umi_5p = umi
            self.umi_5p_qual = umi_qual
        elif (match_spec == "3p"):
            assert (self.start_3p == None)
            self.orientation_3p = orientation
            self.start_3p = match_start_pos
            self.stop_3p = match_stop_pos
            self.umi_3p = umi
            self.umi_3p_qual = umi_qual
#        elif (match_spec == "UDTD5" and self.orientation_udtd5 == None):
#            self.orientation_udtd5 = orientation
#            self.start_udtd5 = match_start_pos
#            self.stop_udtd5 = match_stop_pos
#        elif (match_spec == "UDTD7" and self.orientation_udtd7 == None):
#            self.orientation_udtd7 = orientation
#            self.start_udtd7 = match_start_pos
#            self.stop_udtd7 = match_stop_pos

        return reject_reason

    
    def hasPrimer(self, primer_spec, orientation):
        assert (orientation in ["sense", "antisense"])
        assert (primer_spec in ["5p", "3p"])
        return ((primer_spec == "5p" and self.orientation_5p == orientation) or
                (primer_spec == "3p" and self.orientation_3p == orientation))


    def hasAPrimer(self):
        return (self.orientation_5p != None or self.orientation_3p != None)
    

#    def hasAdapter(self, adapter):
#        assert (adapter == "UDTD5" or adapter == "UDTD7")
#        return ((adapter == "UDTD5" and self.stop_udtd5 != None) or (adapter == "UDTD7" and self.stop_udtd7 != None))


#    def primerOrder(self):
#        assert (self.start_5p != None and self.start_3p != None)
#        primer_order = "5p,3p"
#        if (self.start_5p > self.start_3p):
#            primer_order = "3p,5p"
#        return primer_order


    def trimToPrimers(self):
        assert (self.start_5p != None or self.start_3p != None)
        reject_reason = None

        #print >> sys.stderr, "Before trim: %s" % self.seq
        if (self.orientation_5p == "sense"):
            start = self.start_5p
            if (self.orientation_3p != None):
                if (self.orientation_3p == "antisense"):
                    stop = self.start_3p
                else:
                    reject_reason = "Both primers in same read have same sense"
            else:
                stop = len(self.seq)
        elif (self.orientation_3p == "sense"):
            start = self.start_3p
            if (self.orientation_5p != None):
                if (self.orientation_5p == "antisense"):
                    stop = self.start_5p
                else:
                    reject_reason = "Both primers in same read have same sense"
            else:
                stop = len(self.seq)
        else:
            print >> sys.stderr, "ERROR: tried to trim a read that does not have a sense primer"
            sys.exit(1)
            
        if (reject_reason == None):
            if (stop <= start):
                reject_reason = "Antisense primer occurs before sense primer in read"
            else:
                self.seq = self.seq[start:stop]
                self.qual_string = self.qual_string[start:stop]
                self.qual_scores = self.qual_scores[start:stop]
                #print >> sys.stderr, "After  trim: %s\n" % self.seq

        return reject_reason


    def numPrimersFound(self):
        return int(self.start_5p != None) + int(self.start_3p != None)


    def getUMI(self, primer_spec):
        ret_val = []
        assert (primer_spec == "5p" or primer_spec == "3p")
        if (primer_spec == "5p" and self.umi_5p != None and len(self.umi_5p) > 0):
            ret_val.append(self.umi_5p)
        elif (primer_spec == "3p" and self.umi_3p != None and len(self.umi_3p) > 0):
            ret_val.append(self.umi_3p)
        return ret_val


    def getUMIQual(self, primer_spec):
        ret_val = []
        assert (primer_spec == "5p" or primer_spec == "3p")
        if (primer_spec == "5p" and self.umi_5p_qual != None and len(self.umi_5p_qual) > 0):
            ret_val.append(self.umi_5p_qual)
        elif (primer_spec == "3p" and self.umi_3p_qual != None and len(self.umi_3p_qual) > 0):
            ret_val.append(self.umi_3p_qual)
        return ret_val


    def numNs(self):
        return self.seq.count('N')
    
