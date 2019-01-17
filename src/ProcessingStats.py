class ProcessingStats(object):
    def __init__(self):
        # Total number of read pairs read from reads-to-targets SAM/BAM file
        self.num_pairs_read = 0
        self.num_merged_pair_written = 0

        self.num_pairs_fail_merge = 0
        self.counts_failed_merge_per_target = defaultdict(int)

        # For read pairs where both mates unaligned to any target sequences
        self.num_pairs_none_aligned = 0
        self.num_pairs_none_aligned_merged = 0

        # For read pairs in which the mates align to different target sequences
        self.num_pairs_discordant_target = 0
        self.num_pairs_discordant_target_merged = 0
        
        # For reads pairs in which both pairs align to the same target
        self.num_pairs_both_aligned = 0
        self.num_pairs_both_aligned_merged = 0

        # For reads pairs in which one mate unaligned
        self.num_pairs_one_aligned = 0
        self.num_pairs_one_aligned_merged = 0

        # For read pairs that had internal region soft clipping in their alignments to the target sequence
        # and that did not have any merged-pair or paired-end alignments to the genome
        self.num_mates_to_target_align_clipped_and_unsupported = 0

        # For read pairs that align completely and contiguously (ie CIGAR is just match, '#M')
        self.num_CandC = 0
        self.num_CandC_merged = 0

        # For read pairs that do not align completely and contiguously
        # (ie CIGAR has indel and/or softclipped, or one unmapped to a target)
        self.num_not_CandC = 0
        self.num_not_CandC_merged = 0

        # For read pairs that were found to have a better alignment to a genomic location other
        # that the source of the expected target sequence
        self.better_genomic_for_incomplete = 0
        self.better_genomic_for_complete = 0


    def updateForMergeAttempt(self, target_seq_ID, merges, qc_report):
        merge_successful = len(merges) > 0

        # TODO: handle alt alignments that are as good (ie merges has both target and alt genomic that are equally good...should this every happen?!)
        #best_is_genomic

        if (qc_report.both_mates_unaligned_to_target):
            self.num_pairs_none_aligned += 1
        elif (qc_report.mates_aligned_different_targets):
            self.num_pairs_discordant_target += 1
        elif (qc_report.one_mate_unaligned):
            self.num_pairs_one_aligned += 1
        else:
            self.num_pairs_both_aligned += 1                

        if (qc_report.both_aligned_CandC):
            self.num_CandC += 1
        else:
            self.num_not_CandC += 1

        if (qc_report.mates_to_target_align_clipped_and_unsupported):
            self.num_mates_to_target_align_clipped_and_unsupported += 1
        
        self.num_pairs_read += 1
        if (merge_successful):
            self.num_merged_pair_written += 1
            if (qc_report.both_aligned_CandC):
                self.num_CandC_merged += 1
                self.num_pairs_both_aligned_merged += 1
                #if (best_is_genomic):
                #    self.better_genomic_for_complete += 1
            else:
                self.num_not_CandC_merged += 1
                if (qc_report.both_mates_unaligned_to_target):
                    self.num_pairs_none_aligned_merged += 1 # TODO: currently not reported in output?
                elif (qc_report.mates_aligned_different_targets):
                    self.num_pairs_discordant_target_merged += 1 # TODO: currently not reported in output
                    pdb.set_trace()
                elif (qc_report.one_mate_unaligned):
                    self.num_pairs_one_aligned_merged += 1  # TODO: currently not reported in output?
                else:
                    self.num_pairs_both_aligned_merged += 1                        

                #if (best_is_genomic):
                #    self.better_genomic_for_incomplete += 1
        else:
            self.num_pairs_fail_merge += 1
            #self.counts_failed_merge_per_target[target_seq_ID] += 1 # TODO: currently not handling duplicates. Use set of UMIs instead of integers.

        if (self.num_pairs_read % 100000 == 0):
            print >> sys.stderr, self.num_pairs_read


    def checkAndReport(self, to_filename):
        report_str = ""

        # Confirm that numbers balance
        # self.num_pairs_fail_merge + 
        num_not_written = self.num_pairs_discordant_target + (self.num_pairs_none_aligned - self.num_pairs_none_aligned_merged) + \
          (self.num_pairs_one_aligned - self.num_pairs_one_aligned_merged) + (self.num_pairs_both_aligned - self.num_pairs_both_aligned_merged)

        assert (self.num_pairs_read - self.num_merged_pair_written == num_not_written)
        #assert (self.num_pairs_both_aligned_merged == self.num_not_CandC_merged + self.num_CandC_merged)
        #assert (self.num_merged_pair_written == self.num_pairs_both_aligned_merged + self.num_pairs_one_aligned_merged + self.num_pairs_none_aligned_merged)

        # Report statistics
        perc = 100.0 * float(self.num_merged_pair_written) / float(self.num_pairs_read)
        report_str += "Pairs total, written/read: %d / %d (%5.3f%%)\n" % (self.num_merged_pair_written, self.num_pairs_read, perc)

        perc = 100.0 * float(self.num_pairs_fail_merge) / float(self.num_pairs_read)
        report_str += "Pairs failed merge: %d (%5.3f%% of %d)\n" % (self.num_pairs_fail_merge, perc, self.num_pairs_read)
        
        perc = 100.0 * float(self.num_CandC) / float(self.num_pairs_read)
        report_str += "Pairs CandC: %d (%5.3f%% of %d)\n" % (self.num_CandC, perc, self.num_pairs_read)

        perc = 100.0 * float(self.num_not_CandC) / float(self.num_pairs_read)
        report_str += "Pairs not CandC: %d (%5.3f%% of %d)\n" % (self.num_not_CandC, perc, self.num_pairs_read)

        report_str += "\n"
        
        if (self.num_CandC > 0):
            perc = 100.0 * float(self.num_CandC_merged) / float(self.num_CandC)
            report_str += "Pairs CandC, written/read: %d / %d (%5.3f%%)\n" % (self.num_CandC_merged, self.num_CandC, perc)
        else: 
            report_str += "Pairs CandC, written/read: 0 / 0\n"

        if (self.num_not_CandC > 0):
            perc = 100.0 * float(self.num_not_CandC_merged)/float(self.num_not_CandC)
            report_str += "Pairs not CandC, written/read: %d / %d (%5.3f%%)\n" % (self.num_not_CandC_merged, self.num_not_CandC, perc)
        else:
            report_str += "Pairs not CandC, written/read: 0 / 0\n"            

        # TODO: Expand to include equally-good alt genomic alignments
        #report_str += "Number both aligned - CandC with a better genomic alignment: %d\n" % self.better_genomic_for_complete
        #report_str += "Number both aligned - not CandC with a better genomic alignment: %d\n" % self.better_genomic_for_incomplete

        report_str += "\n"

        if (self.num_pairs_one_aligned > 0):
            perc = 100.0 * float(self.num_pairs_one_aligned_merged) / float(self.num_pairs_one_aligned)
            report_str += "Pairs one unaligned, written/read: %d / %d (%5.3f%%)\n" % (self.num_pairs_one_aligned_merged, self.num_pairs_one_aligned, perc)
        else:
            report_str += "Pairs one unaligned, written/read: 0 / 0\n"

        if (self.num_pairs_none_aligned > 0):
            perc = 100.0 * float(self.num_pairs_none_aligned_merged) / float(self.num_pairs_none_aligned)
            report_str += "Pairs both unaligned, written/read: %d / %d (%5.3f%%)\n" % (self.num_pairs_none_aligned_merged, self.num_pairs_none_aligned, perc)
        else:
            report_str += "Pairs both unaligned, written/read: 0 / 0\n"

        if (self.num_pairs_discordant_target > 0):
            #perc = 100.0 * 0 / float(self.num_pairs_discordant_target)
            report_str += "Pairs with discordant mate targets, written/read: 0 / %d (0%%)\n" % (self.num_pairs_discordant_target) #, perc
        else:
            report_str += "Pairs with discordant mate targets, written/read: 0 / 0\n"

        if (self.num_mates_to_target_align_clipped_and_unsupported > 0):
            #perc = 100.0 * self.num_mates_to_target_align_clipped_and_unsupported
            report_str += "Pairs with clipped target alignments and no genome support, written/read: 0 / %d (0%%)\n" % \
             self.num_mates_to_target_align_clipped_and_unsupported
        else:
            report_str += "Pairs with clipped target alignments and no genome support, written/read: 0 / 0\n"

            
        with open(to_filename, 'w') as op:
            op.write(report_str)

                
    def writeMergeFailureCounts(self, failed_merge_counts):
        target_seq_ID_and_count = sorted(self.counts_failed_merge_per_target.items(), key=itemgetter(1), reverse=True)
        op = open(failed_merge_counts, 'w')
        for tup in target_seq_ID_and_count:
            op.write("%s\t%d\n" % tup)
        op.close()
        
