#!/usr/bin/env python

# Standardizing read order
# ------------------------
# 1) Place both reads in their forward orientation
# 2) Mark locations of 5p, 3p, UDTD5, and UDTD7 in each read in its forward orientation
# 3) Determine which read is the upstream read and which is the downstream read
# 4) Trim the upstream read from its 5p. Return 5p UMI
# 5) Trim the upstream read from its 3p, if it contains a 3p. Return 3p UMI.
# 6) Trim the downstream read from its 3p. Return 3p UMI.
# 7) Trim the downstream read from its 5p, if it contains a 5p. Return 5p UMI.
# 8) Reverse complement the downstream read
# 9) Error-correct the UMIs
# 10) Output: 5p UMI, upstream read, revcomp downstream read, 3p UMI
#
#


from Bio import SeqIO
from collections import defaultdict, Counter
from copy import copy
from operator import itemgetter, xor
import gzip
import os
from string import maketrans, translate
import sys
import pysam
import pdb

import SeqRead


def readPrimerSeqs(primers_fasta):
    all_primer_seqs = {}
    #DNA_complement_table = maketrans("ACGTNacgtn","TGCANtgcan")

    primers_5p, primers_3p = {}, {}
    for record in SeqIO.parse(primers_fasta,'fasta'):
        if (record.id.endswith("_5p") or record.id.endswith("_3p")): # For skipping Illumina adapters
            primer_label, fwd_or_rev = record.id.rsplit('_',1)
            if (fwd_or_rev == "5p"):
                primers_5p[primer_label] = record.seq
            elif (fwd_or_rev == "3p"):
                primers_3p[primer_label] = record.seq

    #revcomp_5p = translate(elems[8], DNA_complement_table)[::-1]
    #revcomp_3p = translate(elems[9], DNA_complement_table)[::-1]

    for primer_label, fwd in primers_5p.items():
        fwd_rc = fwd.reverse_complement()
        rev = primers_3p[primer_label]
        rev_rc = rev.reverse_complement()
        all_primer_seqs[primer_label] = (str(fwd), str(rev), str(fwd_rc), str(rev_rc))

    return all_primer_seqs

#
# upstream read:
#  o Has 5p primer on fwd strand and no (revcomp) 3p or UDTD7
#       - Assert that there is a 3p if there is a UDTD7
#  o Has 5p and then revcomp 3p
#
# downstream read:
#  o Has 3p primer on fwd strand and no (revcomp) 5p or UDTD5
#       - Assert that there is a 5p if there is a UDTD5
#  o Has 3p and then revcomp 5p
#
def standardizeMateOrder(mate1, mate2):
    reject_reason = None
    # Standardize here to make mate1 be the mate with the 5' primer.
    # If both mates have both primers, order does not matter.
    mate1_has_sense_5p = mate1.hasPrimer("5p", "sense")
    mate1_has_sense_3p = mate1.hasPrimer("3p", "sense")
    mate1_has_antisense_5p = mate1.hasPrimer("5p", "antisense")
    mate1_has_antisense_3p = mate1.hasPrimer("3p", "antisense")

    mate2_has_sense_5p = mate2.hasPrimer("5p", "sense")
    mate2_has_sense_3p = mate2.hasPrimer("3p", "sense")
    mate2_has_antisense_5p = mate2.hasPrimer("5p", "antisense")
    mate2_has_antisense_3p = mate2.hasPrimer("3p", "antisense")

    # TODO
    # -If has revcomp 5p/revcomp 3p, then must have 3p/5p. Use for for recovery
    #if (reject_reason == None):
    #    reject_reason = assessAttemptRecovery(mate1, mate2)

    if (mate1_has_sense_5p and mate2_has_sense_3p):
        if (mate1_has_antisense_3p or mate2_has_antisense_5p): # Both should have the other primer if one does
            if (not (mate1_has_antisense_3p and mate2_has_antisense_5p)):
                reject_reason = "Second primer not detected in read(s)"
    elif (mate2_has_sense_5p and mate1_has_sense_3p):
        if (mate2_has_antisense_3p or mate1_has_antisense_5p): # Both should have the other primer if one does
            if (not (mate2_has_antisense_3p and mate1_has_antisense_5p)):
                reject_reason = "Second primer not detected in read(s)"
        mate1, mate2 = mate2, mate1
    #elif (mate1.numPrimersFound() != mate2.numPrimersFound()):
    #    reject_reason = "Mates have unequal number of primer matches"
    else:
        reject_reason = "Uncategorized standardizeMateOrder() failure"
        
    return mate1, mate2, reject_reason


def findMatchStartStop(match_spec, bam_line, reference_length):
    match_start_pos = bam_line.qstart
    match_stop_pos = bam_line.qend
        
    if (bam_line.alen != reference_length):
        # Maximally extend primer
        aligned_pairs = filter(lambda x: x[1] != None, bam_line.get_aligned_pairs())
        if (aligned_pairs[0][1] <= aligned_pairs[0][0] > 0): # 5' end
            match_start_pos -= aligned_pairs[0][1]
        else:
            match_start_pos = 0

        if (aligned_pairs[-1][1] != reference_length-1): # 3' end
            match_stop_pos += reference_length - aligned_pairs[-1][1]
            match_stop_pos = min(match_stop_pos, bam_line.qlen-1)
                
    #elif (match_spec in ["UDTD5", "UDTD7"] and bam_line.alen != reference_length):
    #    # Fix cases where there is up to a 2bp mismatch on the sequencing adapter's 3' end
    #    aligned_pairs = filter(lambda x: x[1] != None, bam_line.get_aligned_pairs())
    #    if (aligned_pairs[-1][1] in [reference_length-3, reference_length-2]):
    #        match_stop_pos += aligned_pairs[-1][1]

    return (match_start_pos, match_stop_pos)


# Primer dimers will have Illumina sequencing adapter in the sequencing reads.
# The non-adapter (e.g. UDTD5 & UDTD7) portion of a read must be long enough to have
# both umis+primers and an intervening sequence longer than 2*umi (~length than possible hairpin).
# If not, it's considered to be a primer dimer.
def checkForPrimerDimer(mate1, mate2, primer_seqs, umi_len):
    reject_reason = None
    #aux_info = None
    
    fwd_primer, rev_primer, revcomp_fwd_primer, revcomp_rev_primer = primer_seqs

    len_fwd_primer = len(fwd_primer)
    len_rev_primer = len(rev_primer)
    len_primers_and_umis = len_fwd_primer + len_rev_primer + 2*umi_len
    
    trimmed_read1_seq, reject_reason = mate1.getAdapterTrimmedSequence()
    if (reject_reason == None):
        trimmed_read2_seq, reject_reason = mate2.getAdapterTrimmedSequence()
    
    if (reject_reason == None):
        mate1_has_sense_5p = mate1.hasPrimer("5p", "sense")
        mate1_has_sense_3p = mate1.hasPrimer("3p", "sense")
        mate1_has_antisense_5p = mate1.hasPrimer("5p", "antisense")
        mate1_has_antisense_3p = mate1.hasPrimer("3p", "antisense")

        mate2_has_sense_5p = mate2.hasPrimer("5p", "sense")
        mate2_has_sense_3p = mate2.hasPrimer("3p", "sense")
        mate2_has_antisense_5p = mate2.hasPrimer("5p", "antisense")
        mate2_has_antisense_3p = mate2.hasPrimer("3p", "antisense")

        mate1_only_5p = (mate1_has_sense_5p or mate1_has_antisense_5p) and (not (mate1_has_sense_3p or mate1_has_antisense_3p))
        mate2_only_5p = (mate2_has_sense_5p or mate2_has_antisense_5p) and (not (mate2_has_sense_3p or mate2_has_antisense_3p))

        mate1_only_3p = (not (mate1_has_sense_5p or mate1_has_antisense_5p)) and (mate1_has_sense_3p or mate1_has_antisense_3p)
        mate2_only_3p = (not (mate2_has_sense_5p or mate2_has_antisense_5p)) and (mate2_has_sense_3p or mate2_has_antisense_3p)

        mate1_has_both_improper = ((mate1_has_sense_5p and mate1_has_sense_3p) or (mate1_has_antisense_3p and mate1_has_antisense_5p))
        mate2_has_both_improper = ((mate2_has_sense_5p and mate2_has_sense_3p) or (mate2_has_antisense_3p and mate2_has_antisense_5p))

        mate1_has_both_proper = ((mate1_has_sense_5p and mate1_has_antisense_3p) or (mate1_has_sense_3p and mate1_has_antisense_5p))
        mate2_has_both_proper = ((mate2_has_sense_5p and mate2_has_antisense_3p) or (mate2_has_sense_3p and mate2_has_antisense_5p))

        # Determine the "primer-combination class"
        if (mate1_only_5p and mate2_only_5p):
            reject_reason = "Primer dimer, only 5' primers in reads"
            #umis_5p = mate1.getUMI("5p") + mate2.getUMI("5p")
            #umis_5p_qual = mate1.getUMIQual("5p") + mate2.getUMIQual("5p")
            #aux_info = mergeUMIs(umis_5p, umis_5p_qual)
        elif (mate1_only_3p and mate2_only_3p):
            reject_reason = "Primer dimer, only 3' primers in reads"
            #umis_3p = mate1.getUMI("3p") + mate2.getUMI("3p")
            #umis_3p_qual = mate1.getUMIQual("3p") + mate2.getUMIQual("3p")
            #aux_info = mergeUMIs(umis_3p, umis_3p_qual)
        elif (mate1_has_both_improper or mate2_has_both_improper):
            reject_reason = "Primer dimer, same sense 5' & 3' primers in read(s)"
        elif (mate1_has_both_proper or mate2_has_both_proper):
            if (min(len(trimmed_read1_seq), len(trimmed_read2_seq)) - len_primers_and_umis <= 2*umi_len):
                reject_reason = "Short primer dimer, proper 5' & 3' primers in read(s)"
        elif (not((mate1_has_sense_5p and mate2_has_sense_3p) or (mate1_has_sense_3p and mate2_has_sense_5p))):
            print >> sys.stderr, "WARNING: unhandled case in checkForPrimerDimer() for %s" % mate1.original_sequencing_readID
            reject_reason = "unhandled case in checkForPrimerDimer()"
            # TODO: RESOLVE: If it's true that UDTD5/UDTD7 is always 5p/3p and UDTD5 is mate1 and UDTD7 is mate2, how can it be that
            # both conditions of this elif() are encountered????

    return reject_reason #, aux_info


def processBamLinesForTarget(target_seq_ID, bam_lines, references, reference_lengths, umi_len):
    reject_reason = None
    mate1, mate2 = None, None

    # Order the bam lines so that primary and better matches are first
    bam_lines = map(lambda x: (x, int(x.seq!=None), x.get_tag("AS"), int(not(x.is_supplementary or x.is_secondary))), bam_lines)
    bam_lines = sorted(bam_lines, key=itemgetter(1,2,3), reverse=True)
    
    for bam_line in map(lambda y: y[0], bam_lines):
        match_spec = references[bam_line.tid]
        reference_length = reference_lengths[bam_line.tid]

        assert (not bam_line.is_unmapped)
        assert (match_spec not in ["UDTD5", "UDTD7"])

        #if (bam_line.qname.startswith("K00180:182:H7VLMBBXX:8:1101:2311:1068") and match_spec[-2:] == "3p"):
            #import pdb
            #pdb.set_trace()

        match_start_pos, match_stop_pos = findMatchStartStop(match_spec, bam_line, reference_length)
        
        assert (target_seq_ID == references[bam_line.tid][0:-3])
        match_spec = references[bam_line.tid][-2:]  

        if (bam_line.qname.endswith("_1")):
            if (mate1 == None):
                mate1 = SeqRead(bam_line.qname[0:-2], bam_line.seq, bam_line.qual, bam_line.query_qualities, bam_line.is_reverse)
            reject_reason = mate1.addMatch(match_spec, match_start_pos, match_stop_pos, bam_line.is_reverse, umi_len)
        elif (bam_line.qname.endswith("_2")):
            if (mate2 == None):
                mate2 = SeqRead(bam_line.qname[0:-2], bam_line.seq, bam_line.qual, bam_line.query_qualities, bam_line.is_reverse)
            reject_reason = mate2.addMatch(match_spec, match_start_pos, match_stop_pos, bam_line.is_reverse, umi_len)
                
        if (reject_reason != None):
            break

    #assert (mate1.hasAPrimer() and mate2.hasAPrimer())

    if (reject_reason == None):
        if (mate1 == None or mate2 == None):
            reject_reason = "mate is None"
        elif (not (mate1.hasPrimer("5p", "sense") and mate2.hasPrimer("3p", "sense")) and    # TODO: This should be reduceable to just the first case
            not (mate1.hasPrimer("3p", "sense") and mate2.hasPrimer("5p", "sense"))):      # TODO: Some cases should be rescueable if one mate has both primer
            reject_reason = "Sense primer(s) missing in mates"
            mate1, mate2 = None, None
        elif (reject_reason == None and (mate1.numNs() > 7 or mate2.numNs() > 7)):
            reject_reason = "Mate w/ too many Ns"
    else:
        assert (mate1 != None and mate2 != None)

    return (reject_reason, mate1, mate2)


def mergeUMIs(umis, umis_qual):
    merged_umi = None
    
    if (len(umis[0]) < len(umis[1])): # one is truncated
        merged_umi = umis[1]
    elif (len(umis[0]) > len(umis[1])): # the other is truncated, 
        merged_umi = umis[0]
    elif (any(map(lambda x: 'N' in x, umis))):
        N_filtered_umis = filter(lambda x: not 'N' in x, umis)
        if (len(N_filtered_umis) == 1):
            merged_umi = N_filtered_umis[0]
        else:
            # Both of the umi's have N's in them. Try to find the correct UMI that doesn't have any N's in it.
            assert (len(N_filtered_umis) == 0)
            corrected_umi = list(umis[0])
            for diff_pos in [i for i in xrange(len(umis[0])) if umis[0][i] != umis[1][i]]:
                if (umis_qual[1][diff_pos] == 'N' or (int(umis_qual[0][diff_pos]) >= int(umis_qual[1][diff_pos]))):
                    corrected_umi[diff_pos] = umis[0][diff_pos]
                else:
                    corrected_umi[diff_pos] = umis[1][diff_pos]
            merged_umi = ''.join(corrected_umi)
    else:
        # Determine the sequence positions that differ and form a corrected UMI based on quality values
        corrected_umi = list(umis[0])
        for diff_pos in [i for i in xrange(len(umis[0])) if umis[0][i] != umis[1][i]]:
            corrected_umi[diff_pos] = umis[0][diff_pos] if (int(umis_qual[0][diff_pos]) >= int(umis_qual[1][diff_pos])) else umis[1][diff_pos]
        merged_umi = ''.join(corrected_umi)

    return merged_umi


def formUMIbasedPairID(target_seq_ID, umis_5p, umis_5p_qual, umis_3p, umis_3p_qual, umi_len):
    final_umis = []
    reject_reason = None

    for curr_umis, curr_umis_qual in [(umis_5p, umis_5p_qual), (umis_3p, umis_3p_qual)]:
        num_distinct_umis = len(set(curr_umis))
        if (num_distinct_umis == 0):
            final_umis.append( 'N' * umi_len )
        elif (num_distinct_umis == 1):
            final_umis.append(curr_umis[0])
        else:
            assert (num_distinct_umis == 2)
            merged_umi = mergeUMIs(curr_umis, curr_umis_qual)
            final_umis.append(merged_umi)

    assert (len(final_umis) == 2)
    return (reject_reason, "%s_%s" % (target_seq_ID, ",".join(final_umis)))


def extractUMIsTrimAndStandardize(target_seq_ID, bam_lines, references, reference_lengths, umi_len, primer_seqs):
    """1) Trim reads upstream of primer sequence(s), relative to primer sequences(s)
       2) Relabel reads <umi_5p>_<umi_3p>
    """
    aux_info = None
    reject_reason, mate1, mate2 = processBamLinesForTarget(target_seq_ID, bam_lines, references, reference_lengths, umi_len)

    if (reject_reason == None):
        mate1, mate2, reject_reason = standardizeMateOrder(mate1, mate2)

    if (reject_reason == None):
        umis_5p = mate1.getUMI("5p") + mate2.getUMI("5p")
        umis_5p_qual = mate1.getUMIQual("5p") + mate2.getUMIQual("5p")
        umis_3p = mate1.getUMI("3p") + mate2.getUMI("3p")
        umis_3p_qual = mate1.getUMIQual("3p") + mate2.getUMIQual("3p")

        reject_reason, umi_based_ID = formUMIbasedPairID(target_seq_ID, umis_5p, umis_5p_qual, umis_3p, umis_3p_qual, umi_len)

        if (reject_reason == None):
            mate1.resetID(umi_based_ID)
            mate2.resetID(umi_based_ID)

            reject_reason = mate1.trimToPrimers()
            if (reject_reason == None):
                mate2.trimToPrimers()
                if (reject_reason == None):
                    reject_reason = performMatesQC(mate1, mate2, primer_seqs, umi_len)
        
    return (reject_reason, mate1, mate2) # aux_info, 

        
def performMatesQC(mate1, mate2, primer_seqs, umi_len):
    """To be called after order of the mates has been standardized and the mates' sequences
    have been trimmed to their primer starts"""
    reject_reason = None
    
    fwd_primer, rev_primer, revcomp_fwd_primer, revcomp_rev_primer = primer_seqs

    read1_seq = mate1.getSequence()
    read2_seq = mate2.getSequence()

    len_fwd_primer = len(fwd_primer)
    len_rev_primer = len(revcomp_rev_primer)
    sum_primer_lens = len_fwd_primer + len_rev_primer

    if (len(read1_seq) < sum_primer_lens or len(read2_seq) < sum_primer_lens):
        reject_reason = "QC fail: trimmed mate too short"

    # TODO: not sure about this
    #else:
    #    read1_starts_w_fwd_primer = len([i for i in xrange(len_fwd_primer) if read1_seq[i] != fwd_primer[i]]) <= 4
    #    read2_starts_w_rev_primer = len([i for i in xrange(len_rev_primer) if read2_seq[i] != rev_primer[i]]) <= 4
    #    
    #    if (not read1_starts_w_fwd_primer or not read2_starts_w_rev_primer):
    #        reject_reason = "QC fail: primer doesn't initiate read"
            
    return reject_reason


# Helper function to selectTwoBestNonOverlappingMatchesPerRead(), below
def selectBestPrimerMatchInOlapGroup(match_group, most_common_target, is_read1, references):
    # Criteria to consider:
    #    a) proper 3p/5p, which depends on which side of mid-read
    #    b) higher alignment score
    #    c) most frequent primer spec
    #
    # Sort on tuple (sum of (not most common target + AS rank + not expected 3p/5p), not most common target, AS rank, not expected 3p/5p)
    
    query_len = -1
    max_AS = -1
    for mg in match_group:
        this_query_len = mg[2].infer_query_length() if (mg[2].query_length == 0) else mg[2].query_length # Infer because secondary alignments have "*" for the sequence
        query_len = this_query_len if (this_query_len > query_len) else query_len
        
        this_AS = mg[2].get_tag('AS')
        max_AS = this_AS if (this_AS > max_AS) else max_AS
        
    assert (query_len > 0)
    mid_read_pos = query_len/2

    matches_with_criteria = []
    for mg in match_group:
        expected_primer_type = "5p" if ((is_read1 and mg[0] < mid_read_pos) or (not is_read1 and mg[1] > mid_read_pos)) else "3p"
        has_expected_primer_type = references[mg[2].tid][-2:] == expected_primer_type
        has_max_AS = mg[2].get_tag('AS') == max_AS
        is_most_common = references[mg[2].tid][0:-3] in most_common_target
        factor_sum = int(is_most_common) + int(has_max_AS) + int(has_expected_primer_type)
        matches_with_criteria.append( (factor_sum, is_most_common, has_max_AS, has_expected_primer_type, mg) )

    matches_with_criteria = sorted(matches_with_criteria, key=itemgetter(0,1,2,3), reverse=True)

    return map(itemgetter(4), matches_with_criteria)[0]


# Helper function to filterSpuriousMatches(), below.
def selectTwoBestNonOverlappingMatchesPerRead(aug_read1_non_udtd_matches, aug_read2_non_udtd_matches, most_common_target, references):
    '''The aug match lists are composed of tuples (start pos, end pos, bam line), where start/end pos are
    relative to the read in its 5'->3' orientation. Lists are expected to be sorted by read alignment end position.
    At least one of the lists is expected to have more than two matches.'''
    read1_two_best, read2_two_best = [], []

    # Group the bam lines/matches that overlap the leftmost and rightmost match of each read and then select from those
    # the bam line(s) that are for the most common target
    if (len(aug_read1_non_udtd_matches) > 1):
        read1_left_group = set(filter(lambda x: x[0] <= aug_read1_non_udtd_matches[0][1], aug_read1_non_udtd_matches))
        read1_right_group = set(filter(lambda x: x[1] >= aug_read1_non_udtd_matches[-1][0], aug_read1_non_udtd_matches))
        if (read1_left_group.isdisjoint(read1_right_group)):
            # Select one from each group
            if (len(read1_left_group) > 1):
                best_match = selectBestPrimerMatchInOlapGroup(read1_left_group, most_common_target, True, references)
                read1_two_best.append( best_match )
            else:
                read1_two_best.extend( read1_left_group )

            if (len(read1_right_group) > 1):
                best_match = selectBestPrimerMatchInOlapGroup(read1_right_group, most_common_target, True, references)
                read1_two_best.append( best_match )
            else:
                read1_two_best.extend( read1_right_group )

        else:
            # Select one from merged group
            best_match = selectBestPrimerMatchInOlapGroup(read1_left_group | read1_right_group, most_common_target, True, references)
            read1_two_best.append( best_match )                

        assert (len(read1_two_best) <= 2)

    else:
        read1_two_best = aug_read1_non_udtd_matches


    if (len(aug_read2_non_udtd_matches) > 1):
        read2_left_group = set(filter(lambda x: x[0] <= aug_read2_non_udtd_matches[0][1], aug_read2_non_udtd_matches))
        read2_right_group = set(filter(lambda x: x[1] >= aug_read2_non_udtd_matches[-1][0], aug_read2_non_udtd_matches))
        if (read2_left_group.isdisjoint(read2_right_group)):
            # Select one from each group
            if (len(read2_left_group) > 1):
                best_match = selectBestPrimerMatchInOlapGroup(read2_left_group, most_common_target, False, references)
                read2_two_best.append( best_match )
            else:
                read2_two_best.extend( read2_left_group )

            if (len(read2_right_group) > 1):
                best_match = selectBestPrimerMatchInOlapGroup(read2_right_group, most_common_target, False, references)
                read2_two_best.append( best_match )
            else:
                read2_two_best.extend( read2_right_group )

        else:
            # Select one from merged group
            best_match = selectBestPrimerMatchInOlapGroup(read2_left_group | read2_right_group, most_common_target, False, references)
            read2_two_best.append( best_match )                

        assert (len(read2_two_best) <= 2)

    else:
        read2_two_best = aug_read2_non_udtd_matches

    read1_two_best = map(itemgetter(2), read1_two_best)
    read2_two_best = map(itemgetter(2), read2_two_best)

    assert (len(read1_two_best) > 0 or len(read2_two_best) > 0)
        
    return (read1_two_best, read2_two_best)


# Helper function to filterSpuriousMatches(), below.
def removeInternalPrimersOfSameType(read1_non_udtd_matches, read2_non_udtd_matches, references):
    '''Lists are expected to be sorted by read alignment end position.
    At least one of the lists is expected to have more than two matches.'''

    # For read1, remove all 5p primers after the first instance and all 3p primers before the last
    first_5p_bam_line, last_3p_bam_line = (None, None), (None, None)
    for i, bam_line in enumerate(read1_non_udtd_matches):
        if (references[bam_line.tid][-2:] == "5p" and first_5p_bam_line[0] == None):
            first_5p_bam_line = (i, bam_line)
        elif (references[bam_line.tid][-2:] == "3p"):
            last_3p_bam_line = (i, bam_line)

    filt_read1_non_udtd_matches = map(itemgetter(1), sorted(filter(lambda x: x[0] != None, [first_5p_bam_line, last_3p_bam_line]), key=itemgetter(0)))

    # For read2, remove all 3p primers after the first instance and all 5p primers before the last
    first_3p_bam_line, last_5p_bam_line = (None, None), (None, None)
    for i, bam_line in enumerate(read2_non_udtd_matches):
        if (references[bam_line.tid][-2:] == "3p" and first_3p_bam_line[0] == None):
            first_3p_bam_line = (i, bam_line)
        elif (references[bam_line.tid][-2:] == "5p"):
            last_5p_bam_line = (i, bam_line)

    filt_read2_non_udtd_matches = map(itemgetter(1), sorted(filter(lambda x: x[0] != None, [first_3p_bam_line, last_5p_bam_line]), key=itemgetter(0)))

    return filt_read1_non_udtd_matches, filt_read2_non_udtd_matches


# All sequencing-read coordinates are relative to the non-reversed read. That is, how it is sequenced 5'->3'.
def filterSpuriousMatches(bam_lines, references):
    reject_reason = None
    passfilter_bam_lines = []
    non_udtd_bam_lines = []
    unmapped_bam_lines = []
    filtered_bam_lines = []    
    max_5p_start_pos = 11 # 0-based
    max_allowed_primer_NM = 3
    max_allowed_UDTD_NM = 5
    read_data_cache = {}      # For bwa alignments that don't include the sequence and quality data
    
    most_common_target = set()
    primer_target_counts = Counter( map(lambda x: references[x.tid][0:-3], bam_lines) )
    most_common_targets = filter(lambda x: not x[0].startswith("UDTD"), primer_target_counts.most_common())
    if (len(most_common_targets)>=1):
        most_common_target.add( most_common_targets[0][0] )
    if (len(most_common_targets)>1):
        i = 1
        while (i < len(most_common_targets) and most_common_targets[i][1]==most_common_targets[0][1]):
            most_common_target.add( most_common_targets[i][0] )    
            i += 1
    try:
        assert (len(most_common_targets) == 2 and most_common_targets[0][1] > most_common_targets[1][1])
    except AssertionError:
        pdb.set_trace()

    # Establish UDTD match for each read, one at most.
    # o If more than one UDTD in read, pick the one with highest AS and remove remainder.
    # o Flag cases where only one read has a UDTD.
    read1_udtd_bam_line, read2_udtd_bam_line = None, None
    read1_udtd_read_start_pos, read1_udtd_read_end_pos = 1e10, 1e10
    read2_udtd_read_start_pos, read2_udtd_read_end_pos = 1e10, 1e10
    for bam_line in bam_lines:
        read_num = bam_line.qname[-2:]
        read_rev = bam_line.is_reverse
        if (bam_line.seq != None and not read_data_cache.has_key( (read_num, read_rev) )):
            read_data_cache[(read_num, read_rev)] = (bam_line.seq, bam_line.qual, bam_line.query_qualities)
        
        if (bam_line.is_unmapped):
            unmapped_bam_lines.append(bam_line)

        elif (references[bam_line.tid].startswith("UDTD")):
            
            if (bam_line.qname.endswith("_1")):
                if (not isinstance(read1_udtd_bam_line, pysam.calignedsegment.AlignedSegment) or bam_line.get_tag('AS') > read1_udtd_bam_line.get_tag('AS')):
                    read1_udtd_bam_line = bam_line
                    read1_udtd_read_start_pos = bam_line.query_length - bam_line.query_alignment_end if (bam_line.is_reverse) else bam_line.query_alignment_start
                    read1_udtd_read_end_pos = bam_line.query_length - bam_line.query_alignment_start if (bam_line.is_reverse) else bam_line.query_alignment_end
            else:
                if (not isinstance(read2_udtd_bam_line, pysam.calignedsegment.AlignedSegment) or bam_line.get_tag('AS') > read2_udtd_bam_line.get_tag('AS')):
                    read2_udtd_bam_line = bam_line
                    read2_udtd_read_start_pos = bam_line.query_length - bam_line.query_alignment_end if (bam_line.is_reverse) else bam_line.query_alignment_start
                    read2_udtd_read_end_pos = bam_line.query_length - bam_line.query_alignment_start if (bam_line.is_reverse) else bam_line.query_alignment_end

        # references[bam_line.tid][0:-3] not in most_common_target and 
        elif (bam_line.get_tag('NM') > max_allowed_primer_NM or (not bam_line.is_reverse and bam_line.qstart - bam_line.reference_start > max_5p_start_pos)):
            filtered_bam_lines.append(bam_line)

        else:
            non_udtd_bam_lines.append(bam_line)
            

    # If only one read has an adapter sequence and the match has many mismatches, consider it to be spurious and discard.
    if (isinstance(read1_udtd_bam_line, pysam.calignedsegment.AlignedSegment) and
        not isinstance(read2_udtd_bam_line, pysam.calignedsegment.AlignedSegment) and read1_udtd_bam_line.get_tag('NM') > max_allowed_UDTD_NM):
        read1_udtd_bam_line = None
        read1_udtd_read_start_pos = 1e10
        read1_udtd_read_end_pos = 1e10
    elif (not isinstance(read1_udtd_bam_line, pysam.calignedsegment.AlignedSegment) and
          isinstance(read2_udtd_bam_line, pysam.calignedsegment.AlignedSegment) and read2_udtd_bam_line.get_tag('NM') > max_allowed_UDTD_NM):
        read2_udtd_bam_line = None        
        read2_udtd_read_start_pos = 1e10
        read2_udtd_read_end_pos = 1e10
    try:
        assert (isinstance(read1_udtd_bam_line, pysam.calignedsegment.AlignedSegment) == isinstance(read2_udtd_bam_line, pysam.calignedsegment.AlignedSegment))
    except AssertionError:
        reject_reason = "Unmatched UDTD in read pair"
        #print >> sys.stderr, "Filtering for %s: %d -> %d" % (bam_lines[0].qname[0:-2], len(bam_lines), len(passfilter_bam_lines))
        return passfilter_bam_lines, reject_reason


    # o If a UDTD in read, keep primer matches that are between read start and beginning of UDTD. Eliminate any primer matches that occur after UDTD start.
    read1_non_udtd_matches = []
    read2_non_udtd_matches = []
    for bam_line in non_udtd_bam_lines:
        #print bam_line.qname[-2:], references[bam_line.tid]
        query_len = bam_line.infer_query_length() if (bam_line.query_length == 0) else bam_line.query_length # Because secondary alignments have "*" for the sequence
        align_start_pos = query_len - bam_line.query_alignment_end if (bam_line.is_reverse) else bam_line.query_alignment_start
        align_end_pos = query_len - bam_line.query_alignment_start if (bam_line.is_reverse) else bam_line.query_alignment_end
        if (bam_line.qname.endswith("_1")  and  align_start_pos < read1_udtd_read_start_pos):
            read1_non_udtd_matches.append( (align_start_pos, align_end_pos, bam_line) )
        elif (bam_line.qname.endswith("_2")  and align_start_pos < read2_udtd_read_start_pos):
            read2_non_udtd_matches.append( (align_start_pos, align_end_pos, bam_line) )


    # With UDTD matches, could have primer dimer, chimera, or short (but real) amplicon. Determine whether it's primer dimer,
    # and let the other two cases be handled below.
    if (isinstance(read1_udtd_bam_line, pysam.calignedsegment.AlignedSegment)): # Both reads have a UDTD match if one does
        # Almost certainly a primer dimer, but check to make sure that there isn't a proper primer combination that would come from a very short amplicon
        reject_reason = "Primer dimer"    
        if (len(read1_non_udtd_matches) > 0 or len(read2_non_udtd_matches) > 0):
            target_primer_instances = defaultdict(set)
            for targetID, primer_5p3p in map(lambda x: (references[x[2].tid][0:-3],references[x[2].tid][-2:]), read1_non_udtd_matches + read2_non_udtd_matches):
                target_primer_instances[targetID].add(primer_5p3p)
            if (any(map(lambda x: len(x)>1, target_primer_instances.values()))):
                reject_reason = None
                #print >> sys.stderr, "UNHANDLED CASE 1"
                #pdb.set_trace()            
                # reject_reason = None

    if (reject_reason == None):
        read1_non_udtd_matches = sorted(read1_non_udtd_matches, key=itemgetter(1))
        read2_non_udtd_matches = sorted(read2_non_udtd_matches, key=itemgetter(1))

        # Make certain there are at most two (nonoverlapping) primer matches per read 
        if (len(read1_non_udtd_matches) > 2 or len(read2_non_udtd_matches) > 2 or
            (len(read1_non_udtd_matches) == 2 and read1_non_udtd_matches[0][1] >= read1_non_udtd_matches[1][0]) or # two, but overlapping
            (len(read2_non_udtd_matches) == 2 and read2_non_udtd_matches[0][1] >= read2_non_udtd_matches[1][0])):  # two, but overlapping
            read1_non_udtd_matches, read2_non_udtd_matches = selectTwoBestNonOverlappingMatchesPerRead(read1_non_udtd_matches, read2_non_udtd_matches,
                                                                                                       most_common_target, references)
        else:
            read1_non_udtd_matches = map(itemgetter(2), read1_non_udtd_matches)
            read2_non_udtd_matches = map(itemgetter(2), read2_non_udtd_matches)
            
        if (len(read1_non_udtd_matches) == 2 or len(read2_non_udtd_matches) == 2):
            read1_non_udtd_matches, read2_non_udtd_matches = removeInternalPrimersOfSameType(read1_non_udtd_matches, read2_non_udtd_matches, references)

        # Handle all the case that can be encountered with a maximum of two non-UDTD matches per read
        if (len(read1_non_udtd_matches) == 0 and len(read2_non_udtd_matches) == 0):
            reject_reason = "Both reads unmapped"

        elif (len(read1_non_udtd_matches) + len(read2_non_udtd_matches) == 1):
            reject_reason = "Read 2 unmapped" if (len(read1_non_udtd_matches)==1) else "Read 1 unmapped"

        elif (len(read1_non_udtd_matches) == 1 and len(read2_non_udtd_matches) == 1):
            read1_match = (references[read1_non_udtd_matches[0].tid][0:-3], references[read1_non_udtd_matches[0].tid][-2:])
            read2_match = (references[read2_non_udtd_matches[0].tid][0:-3], references[read2_non_udtd_matches[0].tid][-2:])

            # If same target and with 3p/5p, then good. If same target with both 3p or both 5p, problem
            if (read1_match[0] == read2_match[0]):
                if (read1_match[1] == "5p" and read2_match[1] == "3p"):
                    passfilter_bam_lines = read1_non_udtd_matches + read2_non_udtd_matches
                else:
                    reject_reason = "Missing complementary primer"
                    
            else: # If different targets and mate1 has 5p and mate2 has 3p, then chimera. Otherwise fallback to mate(s) unmapped?
                if (read1_match[1] == "5p" and read2_match[1] == "3p"):
                    reject_reason = "Chimera"   # TODO? Add info about which primers made Chimera 
                else:
                    reject_reason = "CASE 1,1  B"

        elif (len(read1_non_udtd_matches) == 2 and len(read2_non_udtd_matches) == 2):
            read1_0_match = (references[read1_non_udtd_matches[0].tid][0:-3], references[read1_non_udtd_matches[0].tid][-2:])
            read1_1_match = (references[read1_non_udtd_matches[1].tid][0:-3], references[read1_non_udtd_matches[1].tid][-2:])

            read2_0_match = (references[read2_non_udtd_matches[0].tid][0:-3], references[read2_non_udtd_matches[0].tid][-2:])
            read2_1_match = (references[read2_non_udtd_matches[1].tid][0:-3], references[read2_non_udtd_matches[1].tid][-2:])

            # If both have same target and with 3p/5p, then good.
            if (read1_0_match[0] == read1_1_match[0] and read2_0_match[0] == read2_1_match[0]):
                if (read1_0_match[0] != read2_0_match[0]):
                    reject_reason = "Mates have different targets"
                elif (read1_0_match[1] == "5p" and read1_1_match[1] == "3p" and read2_0_match[1] == "3p" and read2_1_match[1] == "5p"):
                    passfilter_bam_lines.extend( read1_non_udtd_matches )
                    passfilter_bam_lines.extend( read2_non_udtd_matches )
                elif (read1_0_match[1] == "3p" and read1_1_match[1] == "5p" and read2_0_match[1] == "5p" and read2_1_match[1] == "3p"):
                    # TODO: Rescue these cases
                    #passfilter_bam_lines.extend( read1_non_udtd_matches )
                    #passfilter_bam_lines.extend( read2_non_udtd_matches )
                    reject_reason = "Flipped adapters (can rescue)"
                else:
                    reject_reason = "UNHANDLED CASE X"
                    #pdb.set_trace()
                    
            elif (read1_0_match[0] == read2_1_match[0] and read1_1_match[0] == read2_0_match[0]): # TODO? Check 5p/3p
                # If both have the different targets and mate1 has 5p and mate2 has 3p, then chimera. 
                reject_reason = "Chimera"
                #pdb.set_trace()
                
            # TODO: Instead of spiralling into all the cases, maybe just check if len(most_common_target)==1 and
            # any in read1/2_non_udtd_matchs not in most_common and remove those.

            #elif (read1_0_match[0] != read1_1_match[0] and read2_0_match[0] == read2_1_match[0] and read1_0_match[0] == read2_0_match[0]):
            #    # Second match in read1 is spurious. 
            #    assert (read1_0_match[1] == "5p" and read2_0_match[1] == "3p" and read2_1_match[1] == "5p")
            #    passfilter_bam_lines.append( read1_non_udtd_matches[0] )
            #    passfilter_bam_lines.extend( read2_non_udtd_matches )
            #
            #elif (read2_0_match[0] != read2_1_match[0] and read1_0_match[0] == read1_1_match[0] and read2_0_match[0] == read1_0_match[0]):
            #    # Second match in read2 is spurious.
            #    assert (read1_0_match[1] == "5p" and read1_1_match[1] == "3p" and read2_0_match[1] == "3p")
            #    passfilter_bam_lines.extend( read1_non_udtd_matches )
            #    passfilter_bam_lines.append( read2_non_udtd_matches[0] )
            #
            #elif (read1_0_match[0] == read2_0_match[0] and read1_0_match[1] == "5p" and read2_0_match[1] == "3p"):
            #    # Maybe second matches in both reads are both spurious. This is the case if first match on each read are a
            #    # complementary primer pair. Remove bam lines for the spurious.
            #    passfilter_bam_lines.append( [read1_non_udtd_matches[0], read2_non_udtd_matches[0]] )
            #    print >> sys.stderr, "CASE 2,2  A"
            #    pdb.set_trace()

            else:
                reject_reason = "UNHANDLED CASE 3"
                #print >> sys.stderr, "UNHANDLED CASE 3"
                #pdb.set_trace()
                
        elif (len(read1_non_udtd_matches) == 2 and len(read2_non_udtd_matches) == 0):
            read1_0_match = (references[read1_non_udtd_matches[0].tid][0:-3], references[read1_non_udtd_matches[0].tid][-2:])
            read1_1_match = (references[read1_non_udtd_matches[1].tid][0:-3], references[read1_non_udtd_matches[1].tid][-2:])

            # If both have same target and with 3p/5p, then good.
            if (read1_0_match[0] == read1_1_match[0] and read1_0_match[1] == "5p" and read1_1_match[1] == "3p"):
                passfilter_bam_lines.extend( read1_non_udtd_matches )
            else: # Otherwise fallback to mate unmapped?
                reject_reason = "One or both reads unmapped"
                #print >> sys.stderr, "UNHANDLED CASE 5"
                #pdb.set_trace()
                
        elif (len(read1_non_udtd_matches) == 0 and len(read2_non_udtd_matches) == 2): # Conversely
            read2_0_match = (references[read2_non_udtd_matches[0].tid][0:-3], references[read2_non_udtd_matches[0].tid][-2:])
            read2_1_match = (references[read2_non_udtd_matches[1].tid][0:-3], references[read2_non_udtd_matches[1].tid][-2:])

            # If both have same target and with 3p/5p, then good.
            if (read2_0_match[0] == read2_1_match[0] and read2_0_match[1] == "3p" and read2_1_match[1] == "5p"):
                passfilter_bam_lines.extend( read2_non_udtd_matches )
            else: # Otherwise fallback to mate unmapped?
                reject_reason = "One or both reads unmapped"
                #print >> sys.stderr, "UNHANDLED CASE 4"
                #pdb.set_trace()

        elif (len(read1_non_udtd_matches) == 2 and len(read2_non_udtd_matches) == 1):
            read1_0_match = (references[read1_non_udtd_matches[0].tid][0:-3], references[read1_non_udtd_matches[0].tid][-2:])
            read1_1_match = (references[read1_non_udtd_matches[1].tid][0:-3], references[read1_non_udtd_matches[1].tid][-2:])

            read2_0_match = (references[read2_non_udtd_matches[0].tid][0:-3], references[read2_non_udtd_matches[0].tid][-2:])

            # If read1 matches have same target and with 3p/5p, then good. If read2 matches is same target, keep. Otherwise disgard.
            if (read1_0_match[0] == read1_1_match[0] and read1_0_match[1] == "5p" and read1_1_match[1] == "3p"):
                passfilter_bam_lines.extend( read1_non_udtd_matches )
                if (read2_0_match[0] == read1_0_match[0]):
                    passfilter_bam_lines.extend( read2_non_udtd_matches )
            elif (read1_0_match[0] != read1_1_match[0] and read1_0_match[0] == read2_0_match[0]):
                # Maybe second match in read1 is spurious. This is the case if the first matchs in both read are same target and appropriate 5p/3p.
                #     Disgard bam line for second match of read 1
                if (read1_0_match[1] == "3p" and read2_0_match[1] == "5p"):
                    reject_reason = "Flipped adapters (can rescue)"
                else:
                    assert (read1_0_match[1] == "5p" and read2_0_match[1] == "3p")
                    passfilter_bam_lines.extend( [read1_non_udtd_matches[0], read2_non_udtd_matches[0]] )

                # TODO?: Check 5p/3p
                #print >> sys.stderr, "CASE 2,1  A"
                #pdb.set_trace()

            else:
                reject_reason = "Unresolvable"
                #print >> sys.stderr, "UNHANDLED CASE 6"
                #pdb.set_trace()

        elif (len(read1_non_udtd_matches) == 1 and len(read2_non_udtd_matches) == 2): # Conversely
            read1_0_match = (references[read1_non_udtd_matches[0].tid][0:-3], references[read1_non_udtd_matches[0].tid][-2:])

            read2_0_match = (references[read2_non_udtd_matches[0].tid][0:-3], references[read2_non_udtd_matches[0].tid][-2:])
            read2_1_match = (references[read2_non_udtd_matches[1].tid][0:-3], references[read2_non_udtd_matches[1].tid][-2:])

            # If read1 matches have same target and with 3p/5p, then good. If read2 matches is same target, keep. Otherwise disgard.
            if (read2_0_match[0] == read2_1_match[0] and read2_0_match[1] == "3p" and read2_1_match[1] == "5p"):
                passfilter_bam_lines.extend( read2_non_udtd_matches )
                if (read1_0_match[0] == read2_0_match[0]):
                    passfilter_bam_lines.extend( read1_non_udtd_matches )
            elif (read2_0_match[0] != read2_1_match[0] and read2_0_match[0] == read1_0_match[0]):
                # Maybe second match in read2 is spurious. This is the case if the first matchs in both read are same target and appropriate 5p/3p.
                #     Disgard bam line for second match of read 1
                if (read1_0_match[1] == "3p" and read2_0_match[1] == "5p"):
                    reject_reason = "Flipped adapters (can rescue)"
                else:
                    assert (read1_0_match[1] == "5p" and read2_0_match[1] == "3p")
                    passfilter_bam_lines.extend( [read2_non_udtd_matches[0], read1_non_udtd_matches[0]] )
                #print >> sys.stderr, "CASE 1,2  A"
                #pdb.set_trace()
            else:
                reject_reason = "Unresolvable"
                #print >> sys.stderr, "UNHANDLED CASE 7"
                #pdb.set_trace()
            
        else:
            print >> sys.stderr, "SHOULD NOT GET TO THIS CASE"
            pdb.set_trace()

    # Place missing data into bam line objects
    for bam_line in passfilter_bam_lines:
        read_num = bam_line.qname[-2:]
        read_rev = bam_line.is_reverse
        if (bam_line.seq == None):
            if (read_data_cache.has_key( (read_num, read_rev) ) ):
                bam_line.seq, bam_line.qual, bam_line.query_qualities = read_data_cache[(read_num, read_rev)]
            else:
                assert (read_data_cache.has_key( (read_num, not read_rev) ) )
                seq, qual, query_qualities = read_data_cache[(read_num, not read_rev)]
                DNA_complement_table = maketrans("ACGTNacgtn","TGCANtgcan")
                bam_line.seq = translate(seq[::-1], DNA_complement_table)
                bam_line.qual = qual[::-1]
                bam_line.query_qualities = query_qualities[::-1]

    #print >> sys.stderr, "Filtering for %s: %d -> %d" % (bam_lines[0].qname[0:-2], len(bam_lines), len(passfilter_bam_lines))
    return passfilter_bam_lines, reject_reason


def processReadPairAlignments(curr_readID_bam_lines, references, reference_lengths, cooccuring_target_seq_IDs, reject_reason_target_counts, umi_len, all_primer_seqs):
    reject_reason = None # , aux_info, None
    mate1, mate2 = None, None
    target_seq_ID = None
    
    try:
        bam_lines, reject_reason = filterSpuriousMatches(curr_readID_bam_lines, references)
    except AssertionError:
        reject_reason = "Unhandled AssertionError"
        
    if (reject_reason == None):
        primer_targets = set(map(lambda y: references[y.tid][0:-3], filter(lambda x: references[x.tid].endswith("_5p") or references[x.tid].endswith("_3p"), bam_lines)))
        assert (len(primer_targets) == 1)
        target_seq_ID = list(primer_targets)[0]
        primer_seqs = all_primer_seqs[target_seq_ID]
        reject_reason, mate1, mate2 = extractUMIsTrimAndStandardize(target_seq_ID, bam_lines, references, reference_lengths, umi_len, primer_seqs) # aux_info        

    return (reject_reason, target_seq_ID, mate1, mate2) #  aux_info,


def modifyReadPairsUsingAlignments(sam_stream, output_tsv, umi_len, all_primer_seqs): # output_aux, 
    reject_reason_statistics = defaultdict(int)
    reject_reason_target_counts = defaultdict(lambda: defaultdict(int))

    ip_bam = pysam.AlignmentFile(sam_stream, 'r')
    op = gzip.open(output_tsv, 'wb')
    #op_aux = gzip.open(output_aux, 'wb')
    
    references = ip_bam.references
    reference_lengths = ip_bam.lengths
    cooccuring_target_seq_IDs = Counter()
    curr_readID_bam_lines = []
    num_read_pairs_read, num_read_pairs_written = 0,0
    num_read_pairs_w_no_primer_matches, num_read_pairs_w_incomplete_primer_matches = 0,0

    bam_line = ip_bam.next()

    curr_readID = bam_line.qname[0:-2]
    curr_readID_bam_lines.append(bam_line)

    counter = 0
    try:
        while True:
            next_bam_line = ip_bam.next()
            next_readID = next_bam_line.qname[0:-2]

            #if (next_readID != "K00180:182:H7VLMBBXX:8:2216:8410:21729"):
            #    continue
            
            if (next_readID != curr_readID):
                counter += 1
                if (counter%100000 == 0):
                    print >> sys.stderr, "INFO: processed %d read pairs" % counter

                #if (curr_readID.startswith("K00180:182:H7VLMBBXX:8:1101:11221:1068")):
                #    pdb.set_trace()

                # aux_info, 
                reject_reason, target_seq_ID, mate1, mate2 = processReadPairAlignments(curr_readID_bam_lines, references, reference_lengths,
                                                                                       cooccuring_target_seq_IDs, reject_reason_target_counts,
                                                                                       umi_len, all_primer_seqs)
                num_read_pairs_read += 1
                reject_reason_statistics[reject_reason] += 1
                
                if (reject_reason == None):
                    UMI_ID = mate1.getID()
                    mate1_seq = mate1.getSequence()
                    mate1_qual = mate1.getQualString()
                    mate2_seq = mate2.getSequence(True)
                    mate2_qual = mate2.getQualString(True)
                    op.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (curr_readID, UMI_ID, mate1_seq, mate1_qual, mate2_seq, mate2_qual))
                    num_read_pairs_written += 1
                else:
                    pass
                    #print >> sys.stderr, reject_reason
                #elif (aux_info != None):
                #    op_aux.write("%s\t%s\t%s\n" % (reject_reason, target_seq_ID, aux_info))
                    
                curr_readID_bam_lines = [next_bam_line] #.clear()
                curr_readID = next_readID

            else:
                curr_readID_bam_lines.append(next_bam_line)

    except StopIteration:
        # aux_info, 
        reject_reason, target_seq_ID, mate1, mate2 = processReadPairAlignments(curr_readID_bam_lines, references, reference_lengths,
                                                                               cooccuring_target_seq_IDs, reject_reason_target_counts,
                                                                               umi_len, all_primer_seqs)
        num_read_pairs_read += 1
        reject_reason_statistics[reject_reason] += 1
        
        if (reject_reason == None):
            UMI_ID = mate1.getID()
            mate1_seq = mate1.getSequence()
            mate1_qual = mate1.getQualString()
            mate2_seq = mate2.getSequence(True)
            mate2_qual = mate2.getQualString(True)
            op.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (curr_readID, UMI_ID, mate1_seq, mate1_qual, mate2_seq, mate2_qual))
            num_read_pairs_written += 1
        else:
            pass
            #print >> sys.stderr, reject_reason

        #elif (aux_info != None):
        #    op_aux.write("%s\t%s\t%s\n" % (reject_reason, target_seq_ID, aux_info))

    ip_bam.close()
    op.close()
    #op_aux.close()

    return (reject_reason_statistics, num_read_pairs_read, num_read_pairs_written, cooccuring_target_seq_IDs, reject_reason_target_counts)


def writeLogfile(num_read_pairs_read, num_read_pairs_written, reject_reason_statistics, cooccuring_target_seq_IDs, reject_reason_target_counts):
    consistency_count = 0
    
    op = open(logfile, 'w')
    op.write("Total number of read pairs as input: %d\n" % num_read_pairs_read)

    perc_written = 100.0 * float(num_read_pairs_written)/float(num_read_pairs_read)
    op.write("Number of good read pairs written: %d (%4.2f%%)\n" % (num_read_pairs_written, perc_written))
    consistency_count += num_read_pairs_written
    
    rejections_to_list_out = []
    reason_count_tups = reject_reason_statistics.items()
    reason_count_tups = sorted(reason_count_tups, key=itemgetter(1), reverse=True)
    op.write("\nRejection Statistics:\n")
    for reason, count in reason_count_tups:
        if (reason != None):
            perc = 100.0 * float(count)/float(num_read_pairs_read)
            op.write("%s: %d (%4.2f%%)\n" % (reason, count, perc))
            consistency_count += count
            if (perc > 1.0):
                rejections_to_list_out.append( (perc, reason) )

    if (consistency_count != num_read_pairs_read):
        op.write("\nWARNING: Cannot account for all of the read pairs read. Consistency count = %d.\n" % consistency_count)

    rejections_to_list_out = sorted(rejections_to_list_out, key=itemgetter(0), reverse=True)

    for reason in map(lambda x:x[1], rejections_to_list_out):
        op.write("\n\nREJECTION REASON: %s\n" % reason)
        if (len(reject_reason_target_counts[reason]) > 0):
            for target_seq_ID, num_times in sorted(reject_reason_target_counts[reason].items(), key=itemgetter(1), reverse=True):
                op.write("%s\t%d\n" % (target_seq_ID, num_times))

    op.write("\n\nTargets Cooccurring in Same Read Pair\tNumber of Times Seen\n")
    for target_seq_IDs, num_times_seen in sorted(cooccuring_target_seq_IDs.items(), key=itemgetter(1), reverse=True):
        op.write("%s\t%d\n" % (target_seq_IDs, num_times_seen))
    #
    op.close()


if (__name__ == "__main__"):
    umi_len, primers_fasta, sam_stream, output_tsv, logfile = sys.argv[1:] # output_aux, 

    all_primer_seqs = readPrimerSeqs(primers_fasta)
    
    (reject_reason_statistics, num_read_pairs_read, num_read_pairs_written, cooccuring_target_seq_IDs, reject_reason_target_counts) = \
      modifyReadPairsUsingAlignments(sam_stream, output_tsv, int(umi_len), all_primer_seqs) #output_aux, 

    writeLogfile(num_read_pairs_read, num_read_pairs_written, reject_reason_statistics, cooccuring_target_seq_IDs, reject_reason_target_counts)

    sys.exit(0)
    
