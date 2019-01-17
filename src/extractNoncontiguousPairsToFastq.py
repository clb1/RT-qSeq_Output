#!/usr/bin/env python

from collections import defaultdict
from string import maketrans, translate
import gzip
import os
import sys
import pysam

import pdb

if (__name__ == "__main__"):
    input_bam, incomplete_R1_fastq, incomplete_R2_fastq, complete_R1_fastq, complete_R2_fastq = sys.argv[1:] # , log_file

    return_code = 0

    DNA_complement_table = maketrans("ACGTNacgtn","TGCANtgcan")

    ip_bam = pysam.AlignmentFile(input_bam, 'rb')
    op_incompl_r1 = gzip.open(incomplete_R1_fastq, 'wb')
    op_incompl_r2 = gzip.open(incomplete_R2_fastq, 'wb')
    op_compl_r1 = gzip.open(complete_R1_fastq, 'wb')
    op_compl_r2 = gzip.open(complete_R2_fastq, 'wb')

    maxN = 7
    #Ncount_hist = defaultdict(int)
    num_maxN_filtered = 0
    
    try:
        prev_readID = None
        while True:
            # For the two while loops below, can also use the bwa mem parameters
            # -Y            use soft clipping for supplementary alignments
            # -M            mark shorter split hits as secondary
            # and then do
            # while (sam_line1.is_secondary or sam_line1.is_supplementary)
            sam_line1 = ip_bam.next()
            #while (not sam_line1.is_unmapped and 'H' in sam_line1.cigarstring):
            while (sam_line1.is_secondary or sam_line1.is_supplementary):
                sam_line1 = ip_bam.next()

            sam_line2 = ip_bam.next()
            #while (not sam_line2.is_unmapped and 'H' in sam_line2.cigarstring):
            while (sam_line2.is_secondary or sam_line2.is_supplementary):
                sam_line2 = ip_bam.next()

            assert (sam_line1.qname == sam_line2.qname)
            curr_readID = sam_line1.qname

            # Expecting just one target match per read pair
            assert (curr_readID != prev_readID)
             alignment_problem = sam_line1.is_unmapped or sam_line2.is_unmapped or 'S' in sam_line1.cigarstring or 'S' in sam_line2.cigarstring or \
              'I' in sam_line1.cigarstring or 'I' in sam_line2.cigarstring  or 'D' in sam_line1.cigarstring or 'D' in sam_line2.cigarstring  

            Ncount = max(sam_line1.query_sequence.count('N'), sam_line2.query_sequence.count('N'))
            #Ncount_hist[Ncount] += 1

            if (Ncount > maxN):
                num_maxN_filtered += 1
            elif (alignment_problem):
                if (sam_line1.is_read1):
                    op_incompl_r1.write("@%s\n%s\n+\n%s\n" % (curr_readID, sam_line1.query_sequence, sam_line1.qual))
                    seq2_rc = translate(sam_line2.query_sequence, DNA_complement_table)[::-1]
                    qual2_rev = sam_line2.qual[::-1]
                    op_incompl_r2.write("@%s\n%s\n+\n%s\n" % (curr_readID, seq2_rc, qual2_rev))
                else:
                    op_incompl_r1.write("@%s\n%s\n+\n%s\n" % (curr_readID, sam_line2.query_sequence, sam_line2.qual))
                    seq2_rc = translate(sam_line1.query_sequence, DNA_complement_table)[::-1]
                    qual2_rev = sam_line1.qual[::-1]
                    op_incompl_r2.write("@%s\n%s\n+\n%s\n" % (curr_readID, seq2_rc, qual2_rev))
            else:
                if (sam_line1.is_read1):
                    op_compl_r1.write("@%s\n%s\n+\n%s\n" % (curr_readID, sam_line1.query_sequence, sam_line1.qual))
                    seq2_rc = translate(sam_line2.query_sequence, DNA_complement_table)[::-1]
                    qual2_rev = sam_line2.qual[::-1]
                    op_compl_r2.write("@%s\n%s\n+\n%s\n" % (curr_readID, seq2_rc, qual2_rev))
                else:
                    op_compl_r1.write("@%s\n%s\n+\n%s\n" % (curr_readID, sam_line2.query_sequence, sam_line2.qual))
                    seq2_rc = translate(sam_line1.query_sequence, DNA_complement_table)[::-1]
                    qual2_rev = sam_line1.qual[::-1]
                    op_compl_r2.write("@%s\n%s\n+\n%s\n" % (curr_readID, seq2_rc, qual2_rev))

            prev_readID = curr_readID
            sam_line1, sam_line2 = None, None

    except StopIteration:
        ip_bam.close()
        op_incompl_r1.close()
        op_incompl_r2.close()
        op_compl_r1.close()
        op_compl_r2.close()

    if (sam_line1 != None and sam_line2 != None):
        print >> sys.stderr, "ERROR: ended in unexpected state. Removing split fastq files."
        os.unlink(incomplete_R1_fastq)
        os.unlink(incomplete_R2_fastq)
        os.unlink(complete_R1_fastq) 
        os.unlink(complete_R2_fastq)
        return_code = 1
    #else:
    #    total_num = sum(Ncount_hist.values())
    #    perc = 100.0 * float(num_maxN_filtered)/float(total_num)

        #op_log = open(log_file, 'w')
        #op_log.write("INFO: filtered out %d (%5.4f%%) read pairs that exceeded %d Ns per read\n" % (num_maxN_filtered, perc, maxN))
    
        #op_log.write("\nHistogram of read pair N counts:\n")
        #op_log.write("\t N  Num_Pairs\n")
        #Nmax = max(Ncount_hist.keys())
        #for i in xrange(Nmax+1):
        #    if (Ncount_hist[i] > 0):
        #        op_log.write("\t%d: %d\n" % (i, Ncount_hist[i]))
        #    
        #op_log.close()

    sys.exit(return_code)
