#!/usr/bin/env python

from Bio import AlignIO, SeqIO, pairwise2
from StringIO import StringIO
from collections import defaultdict, namedtuple
from itertools import groupby
from operator import attrgetter, itemgetter
from os import unlink
import gzip
import numpy as np
from pyfaidx import Fasta
import pysam
import subprocess
import sys
from tempfile import NamedTemporaryFile

import cPickle

import ProcessingStats

from string import maketrans, translate
DNA_complement_table = maketrans("ACGTNacgtn","TGCANtgcan")

MergeMatesResult = namedtuple('MergeMatesResult', 'targetID, matesID, seq, qual, cigar, edit_dist')
MergeTarget = namedtuple('MergeTarget', 'ID, seq')
ReadPairQC = namedtuple('QC', 'both_mates_unaligned_to_target, mates_aligned_different_targets, one_mate_unaligned, \
                        mates_not_aligned_to_primers, both_aligned_CandC, mates_to_target_align_clipped_and_unsupported, \
                        standardized_reads_aligned_ok')


# Used in a couple of places below:
# CIGAR types:
#     0 = Match
#     1 = Insertion
#     2 = Deletion
#     3 = Skip
#     4 = Soft clipping
#     5 = Hard clipping
#     6 = Padding
# Format of 'cigar' (cigar_type, cigar_len)
conv_cigar = {'0':'M', '1':'I', '2':'D', '3':'N', '4':'S', '5':'H'}    



def parseTargetSeqData(primers_fasta, products_fasta):
    target_seq_data = {}
    primers_5p, primers_3p = {}, {}
    for record in SeqIO.parse(primers_fasta,'fasta'):
        if (record.id.endswith("_5p") or record.id.endswith("_3p")): # For skipping Illumina adapters
            primer_label, fwd_or_rev = record.id.rsplit('_',1)
            if (fwd_or_rev == "5p"):
                primers_5p[primer_label] = record.seq
            elif (fwd_or_rev == "3p"):
                primers_3p[primer_label] = record.seq

    products = {}
    for record in SeqIO.parse(products_fasta,'fasta'):
        product_label = record.id.split(" ")[0]
        #products[record.id] = record.seq
        products[product_label] = record.seq

    assert (set(primers_5p.keys()) == set(primers_3p.keys()))
    try:
        assert (set(primers_5p.keys()) == set(products.keys()))
    except AssertionError, ae:
        print >> sys.stderr, "WARNING: primers and products don't match, probably due to ERCC primers but no ERCC products. FIX."

    for key, product in products.items():
        fwd = primers_5p[key]
        fwd_rc = fwd.reverse_complement()
        rev = primers_3p[key]
        rev_rc = rev.reverse_complement()
        target_seq_data[key] = (str(product).upper(), str(fwd), str(rev), str(fwd_rc), str(rev_rc), len(fwd), len(rev))

    return target_seq_data


def readGenomicReadPairAlignments(bam_files_list, target_seq_data):
    '''Return read pair genomic alignments that are mapped and that offer even a
    faint hope of being real (eg may indicate unknown sequence in the genome).'''
    alt_alignments = defaultdict(list)

    for bam_file in bam_files_list:
        ip_bam = pysam.AlignmentFile(bam_file, 'rb')

        try:
            while True:
                genome_align1 = ip_bam.next()
                while (genome_align1.is_secondary or genome_align1.is_supplementary):
                    genome_align1 = ip_bam.next()

                genome_align2 = ip_bam.next()
                while (genome_align2.is_secondary or genome_align2.is_supplementary):
                    genome_align2 = ip_bam.next()

                assert (genome_align1.qname == genome_align2.qname), "Expected paired reads have different names: %s and %s" % (genome_align1.qname, genome_align2.qname)
                curr_readID = genome_align1.qname
                
                assert (genome_align1.is_read1 != genome_align2.is_read1), "Paired reads not aligned FR, for %s" % curr_readID
                assert (not (genome_align1.is_unmapped or genome_align2.is_unmapped)), "Unmapped pairs should not be in input, for %s" % curr_readID
                if (not genome_align1.is_read1):
                    genome_align1, genome_align2 = genome_align2, genome_align1

                is_softclipped = 'S' in genome_align1.cigarstring or 'S' in genome_align2.cigarstring
                from_same_chrom = genome_align1.tid == genome_align2.tid

                if (not is_softclipped and from_same_chrom):
                    # An alt alignment is only potentially real if the genomic region of the alignment starts
                    # match the expected primer sequences. So beginning of read _1 should have match and end 
                    # of read _2 should have match.
                    # TODO: The primer-matching logic needs to be updated to account for the possibility of a splice in the middle of a primer.
                    #       So, accumulate the contiguous matches across splice junctions from each end and make certain that they are >= fp_len/rp_len.
                    #fp_len = target_seq_data[target_seqID][5]
                    #rp_len = target_seq_data[target_seqID][6]
                    # TODO: Confirm that orientation logic is correct.
                    # and genome_align1.cigar[0][1]>=fp_len,  and genome_align1.cigar[-1][1]>=fp_len , Doesn't work when splice junction in primer
                    #genome_align1_ok = ((not genome_align1.is_reverse and genome_align1.cigar[0][0]==0) or (genome_align1.is_reverse and genome_align1.cigar[-1][0]==0))
                    # and genome_align2.cigar[0][1]>=rp_len,  and genome_align2.cigar[-1][1]>=rp_len
                    #genome_align2_ok = ((not genome_align2.is_reverse and genome_align2.cigar[0][0]==0) or (genome_align2.is_reverse and genome_align2.cigar[-1][0]==0))

                    # Check whether alignment is already recorded (from a previous bam file)
                    aln_already_recorded = False
                    if (alt_alignments.has_key(curr_readID)):
                        genome_align1_positions = genome_align1.positions
                        genome_align2_positions = genome_align2.positions
                        mate1_already_recorded = any(map(lambda x: x[0].positions == genome_align1_positions, alt_alignments[curr_readID]))
                        mate2_already_recorded = any(map(lambda x: x[1].positions == genome_align2_positions, alt_alignments[curr_readID]))
                        aln_already_recorded = mate1_already_recorded or mate2_already_recorded

                    if (not aln_already_recorded): # genome_align1_ok and genome_align2_ok and 
                        # TODO: need to revcomp one of the bam lines?
                        genome_align1.set_tag("XR", genome_align1.reference_name, 'Z')
                        genome_align2.set_tag("XR", genome_align2.reference_name, 'Z')
                        alt_alignments[curr_readID].append( (genome_align1, genome_align2) )

        except StopIteration:
            ip_bam.close()

    return alt_alignments


def readSingletonAlignmentsBam(bam_files_list):
    all_alignments = defaultdict(list)

    for bam_file in bam_files_list:
        ip_bam = pysam.AlignmentFile(bam_file, 'rb')

        try:
            while True:
                alignment = ip_bam.next()
                while ('S' in alignment.cigarstring or alignment.is_secondary or alignment.is_supplementary):
                    alignment = ip_bam.next()

                curr_readID = alignment.qname
                #target_seqID = curr_readID.split(',',1)[1].split('_')[0]

                assert (not alignment.is_unmapped)

                # An alt alignment is only potentially real if the genomic region of the alignment starts
                # match the expected primer sequences. So beginning of read _1 should have match and end 
                # of read _2 should have match.
                # TODO: The primer-matching logic needs to be updated to account for the possibility of a splice in the middle of a primer.
                #       So, accumulate the contiguous matches across splice junctions from each end and make certain that they are >= fp_len/rp_len.
                #fp_len = target_seq_data[target_seqID][5]
                #rp_len = target_seq_data[target_seqID][6]

                # and alignment.cigar[0][1]>=fp_len,  and alignment.cigar[-1][1]>=fp_len , Doesn't work when splice junction in primer
                #alignment_ok = ((not alignment.is_reverse and alignment.cigar[0][0]==0) or (alignment.is_reverse and alignment.cigar[-1][0]==0))

                # Check whether alignment is already recorded (from a previous bam file)
                aln_already_recorded = False
                if (all_alignments.has_key(curr_readID)):
                    alignment_positions = alignment.positions
                    aln_already_recorded = any(map(lambda x: x.positions == alignment_positions, all_alignments[curr_readID]))

                if (not aln_already_recorded): # alignment_ok and 
                    alignment.set_tag("XR", alignment.reference_name, 'Z')
                    all_alignments[curr_readID].append( alignment )

        except StopIteration:
            ip_bam.close()

    return all_alignments


def readSingletonAlignmentsFasta(fasta_files_list):
    '''Assumes USEARCH fasta alignment output format'''
    all_alignments = defaultdict(list)

    for fasta_file in fasta_files_list:
        with gzip.open(fasta_file, 'rb') as ip:
            L = ip.read().split("\n\n")
            for aln in L[0:-1]:
                readID, read_aln, targetID, target_aln = aln.strip().split("\n")
                all_alignments[ (readID[1:], targetID[1:]) ].append( (read_aln, target_aln) )

    return all_alignments


def mergeOverlappingAlignmentPair(curr_readID, target_align1, target_align2, target):
    assert (target_align1.reference_end > target_align2.reference_start) # Do overlap
    align_merge_edit_dist = 0

    # Error-correct the overlapping portions
    olap_len = target_align1.reference_end - target_align2.reference_start
    seqqual1 = zip(target_align1.query_sequence[-olap_len:], target_align1.qual[-olap_len:])
    seqqual2 = zip(target_align2.query_sequence[0:olap_len], target_align2.qual[0:olap_len])
    assert (len(seqqual1) == len(seqqual2)) # Expect this to always be true?

    target_seq_olap_seq = target.seq[target_align2.reference_start:target_align1.reference_end]
    target_seq_olap_qual = []
    assert(len(target_seq_olap_seq)==len(seqqual1))

    # If both reads agree that the target sequence should be different, then record it as changed.
    nuc_changes = []
    for i in xrange(len(target_seq_olap_seq)):
        max_qual = seqqual1[i][1] if (seqqual1[i][1] > seqqual2[i][1]) else seqqual2[i][1]
        target_seq_olap_qual.append( max_qual )
        if (seqqual1[i][0] != target_seq_olap_seq[i] and seqqual2[i][0] != target_seq_olap_seq[i]):
            if (seqqual1[i][0] == seqqual2[i][0]):
                nuc_changes.append( (i, seqqual1[i][0]) )

    target_seq_olap_qual = ''.join(target_seq_olap_qual)

    if (len(nuc_changes) > 0):
        L = list(target_seq_olap_seq)
        for pos, nuc in nuc_changes:
            L[pos] = nuc
        target_seq_olap_seq = ''.join(L)
                
    seq_merged = target_align1.query_sequence[0:-olap_len] + target_seq_olap_seq + target_align2.query_sequence[olap_len:]
    qual_merged = target_align1.qual[0:-olap_len] + target_seq_olap_qual + target_align2.qual[olap_len:]

    if (len(seq_merged) == len(target.seq)):
        align_merge_cigar = "%dM" % len(seq_merged)
        for i in xrange(len(seq_merged)):
            if (seq_merged[i] != target.seq[i]):
                align_merge_edit_dist += 1
    else:
        align_merge_cigar, align_merge_edit_dist = alignToTarget(curr_readID, target.seq, seq_merged)

    return MergeMatesResult(target.ID, curr_readID, seq_merged, qual_merged, align_merge_cigar, align_merge_edit_dist)


def mergeSeparatedAlignmentPair(curr_readID, target_align1, target_align2, target, fake_qual_letter='I'):
    align_merge_edit_dist = 0

    intervening_target_seq = target.seq[target_align1.reference_end:target_align2.reference_start]
    seq_merged = target_align1.query_sequence + intervening_target_seq + target_align2.query_sequence
    qual_merged = target_align1.qual + fake_qual_letter * len(intervening_target_seq) + target_align2.qual

    if (len(seq_merged) == len(target.seq)):
        align_merge_cigar = "%dM" % len(seq_merged)
        for i in xrange(len(seq_merged)):
            if (seq_merged[i] != target.seq[i]):
                align_merge_edit_dist += 1
    else:
        align_merge_cigar, align_merge_edit_dist = alignToTarget(curr_readID, target.seq, seq_merged)

    return MergeMatesResult(target.ID, curr_readID, seq_merged, qual_merged, align_merge_cigar, align_merge_edit_dist)


def mergeAbuttingAlignmentPair(curr_readID, target_align1, target_align2, target):
    align_merge_edit_dist = 0

    # If fail this assertion, have I encountered a situation where reads exactly abutt at a splice junction?"
    assert (target_align1.reference_end == target_align2.reference_start), "ERROR: mates do not exactly abutt on target for %s" % curr_readID
        
    seq_merged = target_align1.query_sequence + target_align2.query_sequence
    qual_merged = target_align1.qual + target_align2.qual

    if (len(seq_merged) == len(target.seq)):
        align_merge_cigar = "%dM" % len(seq_merged)
        for i in xrange(len(seq_merged)):
            if (seq_merged[i] != target.seq[i]):
                align_merge_edit_dist += 1
    else:
        align_merge_cigar, align_merge_edit_dist = alignToTarget(curr_readID, target.seq, seq_merged)
        
    return MergeMatesResult(target.ID, curr_readID, seq_merged, qual_merged, align_merge_cigar, align_merge_edit_dist)


def mergeCompleteTargetAlignments(curr_readID, target_align1, target_align2, target, genome_seq):
    # The two alignments do one of the following, relative to the target sequence: 
    # overlap, leave an intervening gap, or exactly abutt
    merge_result = None

    try:
        assert (target_align1.reference_start==0 and target_align2.reference_end==len(target.seq))

        if (target_align1.reference_end > target_align2.reference_start): # Do overlap
            merge_result = mergeOverlappingAlignmentPair(curr_readID, target_align1, target_align2, target)
        elif (target_align1.reference_end < target_align2.reference_start): # Intervening gap >= 1bp
            merge_result = mergeSeparatedAlignmentPair(curr_readID, target_align1, target_align2, target)
        else: # Exactly abutt
            merge_result = mergeAbuttingAlignmentPair(curr_readID, target_align1, target_align2, target)
    except AssertionError:
        print >> sys.stderr, "WARNING: complete alignments to target didn't align to ends of target sequence, %s" % curr_readID
        #pdb.set_trace()

    return merge_result


def alignToTarget(curr_readID, target_seq, other_seq):
    '''The returned CIGAR is for other_seq relative to target_seq.'''

    edit_distance = 0
    cigar = ''
    
    alignments = pairwise2.align.localms(other_seq, target_seq, 1,-2,-6,-1)
    other_alignment = alignments[0][0]
    target_alignment = alignments[0][1]

    for i in xrange(len(other_alignment)):
        if (other_alignment[i] == '-'):
            edit_distance += 1
            cigar += 'D'
        elif (target_alignment[i] == '-'):
            edit_distance += 1
            cigar += 'I'
        else:
            cigar += 'M'
            if (other_alignment[i] != target_alignment[i]):
                edit_distance += 1
                
    cigar_tuples = [(len(list(group)),cigar_type) for cigar_type, group in groupby(cigar)]
    assert (cigar_tuples[0][1]=='M' and cigar_tuples[-1][1]=='M'), "ERROR: candidate merged seq for %s expected to match ends of target sequence" % curr_readID
    cigar_string = ''.join(map(lambda x: "%d%s" % x, cigar_tuples))

    return (cigar_string, edit_distance)


def reconstructSequenceFromMSA(targetID, curr_readID, msa, mate1_qual, mate2_qual, patch_intermate_gap=False):
    '''Expects the alignments to be ordered mate1, mate2, target.'''
    # Merging heuristics:
    #   o Require both mates to agree on the nucleotide for an insertion when in overlap
    #   o Require both mates to agree on the position of a deletion when in overlap
    #   o NO LONGER APPLIED: Don't accept indels within primer region
    #   o NO LONGER APPLIED: In final alignment, must be the case if # resolved discrepancies > # unresolved discrepancies. Otherwise alignment deemed over-forced.
    global conv_cigar
    merge_result = None
    S, Q, C = '', '', []
    m1_ind, m2_ind, t_ind = -1, -1, -1
    edit_distance = 0

    # Find columns delimiting overlap region
    num_cols = msa.get_alignment_length()
    olap_cols = map(lambda w: w[0], filter(lambda y: y[1][0]!='-' and y[1][1]!='-', map(lambda x: (x,msa[:,x]), range(num_cols))))
    olap_start, olap_stop = (olap_cols[0], olap_cols[-1]) if (len(olap_cols)>0) else (-1,-1)
    mate_alignments_overlap = len(olap_cols) > 0

    #primer_cols = set(range(18) + range(num_cols-18, num_cols))
    for c in xrange(num_cols):
        in_olap = c >= olap_start and c <= olap_stop
        #in_primer = c in primer_cols
        assert (msa[:,c] != '---'), "MSA in reconstructSequenceFromMSA() has column that is all '-', for readID %s" % curr_readID
        m1,m2,t = list(msa[:,c])
        if (m1 != '-'): m1_ind += 1  # These three are indices into the nucleotide sequences, not into MSA positions. That would be 'c'.
        if (m2 != '-'): m2_ind += 1
        if (t != '-'):   t_ind += 1
        if (t != '-'):
            if (m1 == m2 == t):
                S += t
                assert (mate1_qual[m1_ind]!=None and mate2_qual[m2_ind]!=None)
                Q += max(mate1_qual[m1_ind], mate2_qual[m2_ind])
                C.append('M')
            elif (m1 == t):
                S += t
                assert (mate1_qual[m1_ind]!=None)
                Q += mate1_qual[m1_ind]
                C.append('M')
            elif (m2 == t):
                S += t
                assert (mate2_qual[m2_ind]!=None)
                Q += mate2_qual[m2_ind]
                C.append('M')
            elif (m1 == m2 == '-'):
                if (not in_olap and not mate_alignments_overlap and patch_intermate_gap):    # Fill in the gap between mates whose alignments do not overlap
                    S += t
                    Q += 'I'
                    C.append('M')
                else:
                    C.append('D')
                    edit_distance += 1
            elif (m1 == m2):
                S += m1
                assert (mate1_qual[m1_ind]!=None and mate2_qual[m2_ind]!=None)
                Q += max(mate1_qual[m1_ind], mate2_qual[m2_ind])
                C.append('M')
                edit_distance += 1
            elif (m1 == '-'):
                assert (m2 != t), "reconstructSequenceFromMSA(), (m1 == '-'), c = %d" % c
                S += m2
                assert (mate2_qual[m2_ind]!=None)
                Q += mate2_qual[m2_ind]
                C.append('M')
                edit_distance += 1
            elif (m2 == '-'):
                assert (m1 != t), "reconstructSequenceFromMSA(), (m2 == '-'), c = %d" % c
                S += m1
                assert (mate1_qual[m1_ind]!=None)
                Q += mate1_qual[m1_ind]
                C.append('M')
                edit_distance += 1
            elif (m1 != t and m2 != t and m1 != m2):
                S += t
                Q += 'I'
                C.append('M')
            else:
                print >> sys.stderr, "ERROR: shouldn't have reached this case in reconstructSequenceFromMSA() for %s. Exiting." % curr_readID
                sys.exit(1)
        else: # t == '-'
            if (m1 == m2):
                S += m1
                assert (mate1_qual[m1_ind]!=None and mate2_qual[m2_ind]!=None)
                Q += max(mate1_qual[m1_ind], mate2_qual[m2_ind])
                C.append('I')
                edit_distance += 1
            elif (m1 == '-'): 
                if (not in_olap): # Skip if m1 and target agree on '-', but m2 has an insert while in overlap region
                    S += m2
                    assert (mate2_qual[m2_ind]!=None)
                    Q += mate2_qual[m2_ind]
                    C.append('I')
                    edit_distance += 1
            elif (m2 == '-'):
                if (not in_olap): # Skip if m2 and target agree on '-', but m1 has an insert while in overlap region
                    S += m1
                    assert (mate1_qual[m1_ind]!=None)
                    Q += mate1_qual[m1_ind]
                    C.append('I')
                    edit_distance += 1
            elif (m1 != m2):
                assert (in_olap), "ERROR: should be  in overlap if mate alignment symbols are both not '-'. curr_readID = %s" % curr_readID
                assert (mate1_qual[m1_ind]!=None and mate2_qual[m2_ind]!=None)
                S += m1 if (mate1_qual[m1_ind] > mate2_qual[m2_ind]) else m2
                Q += max(mate1_qual[m1_ind], mate2_qual[m2_ind])
                C.append('I')
                edit_distance += 1
            else:
                print >> sys.stderr, "ERROR: shouldn't have reached this case in reconstructSequenceFromMSA() for %s. Exiting." % qname
                sys.exit(1)

    cigartuples = [(len(list(group)),cigar_type) for cigar_type, group in groupby(C)]
    cigarstring = ''.join(map(lambda x: "%d%s" % x, cigartuples))
    merge_result = MergeMatesResult(targetID, curr_readID, S, Q, cigarstring, edit_distance)

    return merge_result


def createBed12Line(chrom, strand, all_genome_positions):
    all_genome_positions.sort()
    first_pos = all_genome_positions[0]

    block_sizes = []
    block_starts = []
    for k, g in groupby(enumerate(all_genome_positions), lambda (i,x):i-x):
        group = map(itemgetter(1), g)
        block_starts.append( str(group[0] - first_pos) )
        block_sizes.append( str(group[-1] - group[0] + 1) )
        
    bed12_line = "%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\%d\t%s\t%s" % \
      (chrom, , , , strand, , , "0,0,0", len(block_starts), ','.join(block_starts), ','.join(block_sizes))

    
    return bed12_line


def mergeOverlappingPairwiseAlignments(target_align1, target_align2):
    S, Q, C, edit_dist = '','','',0
    expanded_cigar1 = dict(zip(target_align1.get_aligned_pairs(), list(''.join(map(lambda x: conv_cigar[str(x[0])] * x[1], target_align1.cigartuples)))))
    expanded_cigar2 = dict(zip(target_align2.get_aligned_pairs(), list(''.join(map(lambda x: conv_cigar[str(x[0])] * x[1], target_align2.cigartuples)))))

    # TODO: Make sure this is done in the upstream functions()
    #assert (mate1_chrom == mate2_chrom), "ERROR: unhandled case of better genomic alignment being interchromosomal in altGenomicAlignments()"

    mate1_seq, mate1_qual = target_align1.query_sequence, target_align1.qual
    mate2_seq, mate2_qual = target_align2.query_sequence, target_align2.qual

    ref_nucs = dict(zip(target_align1.get_reference_positions(), target_align1.get_reference_sequence().upper()))
    ref_nucs.update( dict(zip(target_align2.get_reference_positions(), target_align2.get_reference_sequence().upper())) )

    D1,D2 = {},{}
    for D, aligned_pairs, expanded_cigar, mate_seq, mate_qual in [(D1, target_align1.get_aligned_pairs(), expanded_cigar1, mate1_seq, mate1_qual),
                                                                  (D2, target_align2.get_aligned_pairs(), expanded_cigar2, mate2_seq, mate2_qual)]:
        prev_ref_pos = None
        for query_pos, ref_pos in aligned_pairs:
            cigar_letter = expanded_cigar[ (query_pos,ref_pos) ]
            if (cigar_letter == 'M'): # query and reference match
                D[ref_pos] = [(mate_seq[query_pos], mate_qual[query_pos], ref_nucs[ref_pos], cigar_letter)]
                prev_ref_pos = ref_pos
            elif (cigar_letter in ['D','N']): # Deletion or intron skip in query relative to reference
                assert (query_pos == None)
                D[ref_pos] = [ (None, None, 'X', cigar_letter) ] # 'X' is placeholder. There is a nucleotide, but don't need to know it, below.
                prev_ref_pos = ref_pos
            elif (cigar_letter == 'I'): # Insert in query relative to reference
                assert(ref_pos == None)
                D[prev_ref_pos].append( (mate_seq[query_pos], mate_qual[query_pos], None, cigar_letter) )
            elif (cigar_letter != 'S'):
                print >> sys.stderr, "ERROR: unexpected case in preprocessing step of mergeOverlappingPairwiseAlignments(). Exiting."
                sys.exit(1)
                
    all_ref_pos = sorted(set(D1.keys() + D2.keys())) # [None] + 
    for ref_pos in all_ref_pos:
        if (D1.has_key(ref_pos) and D2.has_key(ref_pos)):
            # TODO: This case will need to be modified to handle 'S'/'M' conflicts and probably different number of different number of indels.
            # Solution is probably to upfront resolve differences into one list of [(nuc, qual, ref_nuc, cigar)], and then proceed as in elif() cases below.
            try:
                assert (len(D1[ref_pos]) == len(D2[ref_pos]))
            except AssertionError:
                pdb.set_trace()
            for i in xrange(len(D1[ref_pos])):
                nuc1, qual1, ref_nuc, cigar1 = D1[ref_pos][i]
                nuc2, qual2, ref_nuc, cigar2 = D2[ref_pos][i]
                try:
                    assert (cigar1 == cigar2), "CIGARs not equal"
                except AssertionError:
                    pdb.set_trace()

                C += cigar1
                if (ref_nuc == None): # Insert in queries relative to reference
                    assert (nuc1 != None and nuc2 != None), "Neither query alignments should be None"
                    edit_dist += 1
                    if (qual1 > qual2):
                        S += nuc1
                        Q += qual1
                    else:
                        S += nuc2
                        Q += qual2
                elif (nuc1 == None or nuc2 == None): # Deletion or intron skip in queries relative to reference
                    assert (nuc1 == nuc2 == None), "Both query alignments should be None"
                    if (cigar1 == 'D'):
                        edit_dist += 1
                else:
                    if (nuc1 == nuc2):
                        S += nuc1
                        Q += max(qual1,qual2)
                        if (nuc1 != ref_nuc):
                            edit_dist += 1
                    elif (nuc1 == ref_nuc):
                        S += nuc1
                        Q += qual1
                    elif (nuc2 == ref_nuc):
                        S += nuc2
                        Q += qual2
                    else:
                        S += ref_nuc
                        Q += min(qual1,qual2)
                        edit_dist += 1

        elif (D1.has_key(ref_pos)):
            for nuc1, qual1, ref_nuc, cigar1 in D1[ref_pos]:
                C += cigar1
                if (cigar1 == 'S'):
                    assert (nuc1 != None and ref_nuc == None)
                    edit_dist += 1
                elif (cigar1 == 'I'):
                    assert (ref_nuc == None)
                    S += nuc1
                    Q += qual1
                    edit_dist += 1
                elif (cigar1 == 'D'):
                    assert (nuc1 == None and ref_nuc != None)
                    edit_dist += 1
                elif (cigar1 == 'M'):
                    assert (nuc1 != None and ref_nuc != None)
                    S += nuc1
                    Q += qual1
                    if (nuc1 != ref_nuc):
                        edit_dist += 1
                elif (cigar1 != 'N'):
                    print >> sys.stderr, "ERROR: unexpected cigar1 in mergeOverlappingPairwiseAlignments(). Exiting."
                    pdb.set_trace()
                    sys.exit(1)

        elif (D2.has_key(ref_pos)):
            for nuc2, qual2, ref_nuc, cigar2 in D2[ref_pos]:
                C += cigar2
                if (cigar2 == 'S'):
                    assert (nuc2 != None and ref_nuc == None)
                    edit_dist += 1
                elif (cigar2 == 'I'):
                    assert (ref_nuc == None)
                    S += nuc2
                    Q += qual2
                    edit_dist += 1
                elif (cigar2 == 'D'):
                    assert (nuc2 == None and ref_nuc != None)
                    edit_dist += 1
                elif (cigar2 == 'M'):
                    assert (nuc2 != None and ref_nuc != None)
                    S += nuc2
                    Q += qual2
                    if (nuc2 != ref_nuc):
                        edit_dist += 1
                elif (cigar2 != 'N'):
                    print >> sys.stderr, "ERROR: unexpected cigar2 in mergeOverlappingPairwiseAlignments(). Exiting."
                    sys.exit(1)
        elif (ref_pos != None):
            print >> sys.stderr, "ERROR: unexpected case in mergeOverlappingPairwiseAlignments(). Exiting."
            pdb.set_trace()
            sys.exit(1)
            
    # Add soft-clipping to cigar and to edit distance.
    # For left:
    #   - Get number of soft-clipped positions at beginning of left mate
    #   - Subtract the number of (match) positions in right mate that are less than the smallest position in left mate.
    #   - Add that number of 'S' to beginning of CIGAR and edit distance.
    if (target_align1.cigartuples[0][0] == 4):
        numS = target_align1.cigartuples[0][1]
        numS_covered = max(0, min(target_align1.positions) - min(target_align2.positions))
        numS_to_add = 0 if (numS_covered >= numS) else numS - numS_covered
        C = 'S' * numS_to_add + C
        edit_dist += numS_to_add

    # For right:
    #   - Get number of soft-clipped positions at end of right mate
    #   - Subtract the number of (match) positions in left mate that are greater then the largest position in right mate.
    #    - Add that number of 'S' to end of CIGAR and edit distance.
    if (target_align2.cigartuples[-1][0] == 4):
        numS = target_align2.cigartuples[-1][1]
        numS_covered = max(0, max(target_align1.positions) - max(target_align2.positions))
        numS_to_add = 0 if (numS_covered >= numS) else numS - numS_covered
        C = C + 'S' * numS_to_add
        edit_dist += numS_to_add
    
    compact_cigar = ''
    for k, g in groupby(C):
        compact_cigar += "%d%s" % (len(list(g)), k)

    return (S, Q, compact_cigar, edit_dist)



def altGenomicAlignments2(targetID, curr_readID, target_align1, target_align2, curr_edit_dist, alt_alignments, alt_merged_to_genome, genome_seq):
    ''' Nonstandard CIGAR format is <chr>:<strand>:<reference left pos>-<reference right pos>:<CIGAR>
        Sequence, quality string, and CIGAR are *NOT* reported in the 5'->3' orientation for those on the '-' strand.

        <reference start pos> and <reference stop pos> are in the context of the forward strand regardless of the actual strand of the genomic_target_strand

        Assumption: The alternative alignments (alt_alignments) provided to this function are to genomic regions that are not the source
                    of the expected target.
    '''
    alt_genomic_alignments = []

    have_alt_merged = alt_merged_to_genome.has_key(curr_readID)
    have_alt_paired = alt_alignments.has_key(curr_readID)
    assert (have_alt_merged or have_alt_paired)

    alt_paired = alt_alignments[curr_readID]
    alt_merged = alt_merged_to_genome[curr_readID]

    # TODO: Confirm logic used here.
    # If merged and paired genome alignments are disjoint, record both
    # Else
    #    If both alignments are same genome positions, only use merged
    #    Else If the paired and merged alignments are to different places in the genome, use both
    #    Else If one is a subset of the other, somehow reconcile with realignment. Or, just use the one with the smaller edit distance.
    #    Else the two alignments must only partially overlap. Somehow reconcile, or use the one with the smaller edit distance.
    if (have_alt_merged and have_alt_paired):
        # Try to remove paired-end alignments that are already accounted for by merged-pair alignments
        merged_positions = map(lambda x: set(x.positions), alt_merged)
        paired_positions = map(lambda x: set(x[0].positions + x[1].positions), alt_paired)

        to_remove_indices = set()
        for i, pp in enumerate(paired_positions):
            if (any(map(lambda mp: pp==mp, merged_positions))): # paired alignment is identical to a merged alignment
                to_remove_indices.add(i)
            elif (all(map(lambda mp: pp.isdisjoint(mp), merged_positions))): # paired alignment has no positions in common with a merged alignment
                pass
                #print >> sys.stderr, "Separate genome alignments found for paired vs merged, %s" % curr_readID                
            elif (('S' in alt_paired[i][0].cigarstring or 'S' in alt_paired[i][1].cigarstring) and any(map(lambda mp: pp.is_subset(mp), merged_positions))):
                # Remove soft-clipped paired-end alignments that are subsets of a merged-pair alignment.
                to_remove_indices.add(i)

        alt_paired = [y for index, y in enumerate(alt_paired) if index not in to_remove_indices]
        alt_paired_edit_dist = [y for index, y in enumerate(alt_paired_edit_dist) if index not in to_remove_indices]
        paired_positions = [y for index, y in enumerate(paired_positions) if index not in to_remove_indices]
        have_alt_paired = len(alt_paired) > 0

        for j, mp in enumerate(merged_positions):
            assert (not any(map(lambda x: mp == x, paired_positions)))
            if (any(map(lambda x: mp < x, paired_positions))):
                print >> sys.stderr, "WARNING: Merged alignment is proper subset of paired-end alignment for %s" % curr_readID

    # Rank the alignments by edit distance
    # Get the best paired-end alignment(s)
    if (have_alt_paired):
        L = []
        for pe_align1, pe_align2 in alt_paired:
            edit_dist = sum(map(lambda y: y[1], filter(lambda x: x[0]==4, pe_align1.cigartuples))) + pe_align1.get_tag("NM") + \
              sum(map(lambda y: y[1], filter(lambda x: x[0]==4, pe_align2.cigartuples))) + pe_align2.get_tag("NM")
            L.append( ((pe_align1, pe_align2), edit_dist) )
        min_edit_dist = min(map(itemgetter(1), L))
        L = filter(lambda x: x[1]==min_edit_dist, L)
        #assert (len(L) == 1), "Multiple equally good alternate paired-end alignments in different genome locations for %s in altGenomicAlignments2()" % curr_readID
        alt_paired = map(itemgetter(0), L)
        alt_paired_edit_dist = map(itemgetter(1), L)

    # Get the best merged alignment(s)
    if (have_alt_merged):
        L = []
        for align in alt_merged:
            edit_dist = sum(map(lambda y: y[1], filter(lambda x: x[0]==4, align.cigartuples))) + align.get_tag("NM") 
            L.append( (align, edit_dist) )
        min_edit_dist = min(map(itemgetter(1), L))
        L = filter(lambda x: x[1]==min_edit_dist, L)
        #assert (len(L) == 1), "Multiple equally good alternate merged alignments in different genome locations for %s in altGenomicAlignments2()" % curr_readID
        alt_merged = map(itemgetter(0), L)
        alt_merged_edit_dist = map(itemgetter(1), L)

    if (have_alt_merged):
        for i, genome_align in enumerate(alt_merged):
            if (genome_align.cigartuples[0][0]==0 and genome_align.cigartuples[-1][0]==0):
                S, Q, C = genome_align.query_sequence, genome_align.qual, genome_align.cigarstring

                chrom = genome_align.get_tag("XR")
                genome_strand = '-' if (genome_align.is_reverse) else '+'                
                cigarstring = "%s:%s:%d-%d:%s" % (chrom, genome_strand, min(genome_align.positions)+1, max(genome_align.positions)+1, C)

                merge_result = MergeMatesResult(targetID, curr_readID, S, Q, cigarstring, alt_merged_edit_dist[i])
                alt_genomic_alignments.append(merge_result)

    if (have_alt_paired):
        for i, (genome_align1, genome_align2) in enumerate(alt_paired):
            pairwise_alignments_overlap = len(set(genome_align1.positions) & set(genome_align2.positions)) > 0
            if (pairwise_alignments_overlap):
                #print >> sys.stderr, "Paired read genomic positions overlap, but no merged pair genome alignment, %s" % curr_readID
                # This happens when legitimate alternate genome alignments to pseudogenes are found for the pairs but not for the merged. It's at least partly
                # due to HISAT parameterizations and the interplay of mismatch splicing score penalties for both the paired-end and merged-pair alignments.
                #
                # TODO: employ same logic as in resolveIncompleteAlignmentsToTarget() to create merged sequence, and then add to alt_genome_alignments?
                # In all cases so far this seems to be a HISAT2 bug. (Not sure this comment now relevant.)
                #pdb.set_trace()
                pass
            if (genome_align1.cigartuples[0][0]==0 and genome_align2.cigartuples[-1][0]==0):
                S = "%s,%s" % (genome_align1.query_sequence, genome_align2.query_sequence)
                Q = "%s,%s" % (genome_align1.qual, genome_align2.qual)

                mate1_chrom = genome_align1.get_tag("XR")
                mate2_chrom = genome_align2.get_tag("XR")
                try:
                    assert (mate1_chrom == mate2_chrom), "Unhandled case of better genomic alignment being interchromosomal"
                    assert (genome_align1.is_reverse != genome_align2.is_reverse), "Paired-end reads not mapped to opposite strands as expected"

                    genome_strand = '-' if (genome_align1.is_reverse) else '+'

                    paired_positions = set(genome_align1.positions + genome_align2.positions)
                    cigarstring = "%s:%s:%d-%d:%s,%s" % (mate1_chrom, genome_strand, min(paired_positions)+1, max(paired_positions)+1,
                                                        genome_align1.cigarstring, genome_align2.cigarstring)

                    merge_result = MergeMatesResult(targetID, curr_readID, S, Q, cigarstring, alt_paired_edit_dist[i])
                    alt_genomic_alignments.append(merge_result)
                    # TODO: For paired-end alignments that don't overlap but do end within the same exon will be stitched together. Here is where that needs to happen.
                except AssertionError, ae:
                    print >> sys.stderr, "AssertionError: %s. For %s. Continuing." % (ae.message, curr_readID)
                    
    # Reconcile merged and paired alignments
    if (len(alt_genomic_alignments) > 1):
        min_edit_dist = min(map(attrgetter("edit_dist"), alt_genomic_alignments))
        alt_genomic_alignments = filter(lambda x: x.edit_dist <= min_edit_dist, alt_genomic_alignments)

    return alt_genomic_alignments


def altGenomicAlignmentsFromPairedReads(targetID, curr_readID, target_align1, target_align2, curr_edit_dist, alt_alignments, genome_seq):
    ''' Nonstandard CIGAR format is <chr>:<strand>:<reference left pos>-<reference right pos>:<CIGAR>
        Sequence, quality string, and CIGAR are *NOT* reported in the 5'->3' orientation for those on the '-' strand.

        <reference start pos> and <reference stop pos> are in the context of the forward strand regardless of the actual strand of the genomic_target_strand

        Assumption: The alternative alignments (alt_alignments) provided to this function are to genomic regions that are not the source
                    of the expected target.
    '''
    alt_genomic_alignments = []

    for genome_align1, genome_align2 in alt_alignments[curr_readID]:
        pairwise_alignments_overlap = len(set(genome_align1.positions) & set(genome_align2.positions)) > 0
        all_genome_positions = set(genome_align1.positions + genome_align2.positions)

        mate1_seq, mate1_qual = genome_align1.query_sequence, genome_align1.qual
        mate2_seq, mate2_qual = genome_align2.query_sequence, genome_align2.qual

        if (pairwise_alignments_overlap):
            pdb.set_trace()
            # 1) Extract genomic sequence based on min/max alignment positions
            #mate1_chrom = alt_alignments['references'][genome_align1.tid]
            #mate2_chrom = alt_alignments['references'][genome_align2.tid]
            assert (mate1_chrom == mate2_chrom), "Unhandled case of better genomic alignment being interchromosomal"
            assert (genome_align1.is_reverse != genome_align2.is_reverse)
            genome_strand = '-' if (genome_align1.is_reverse) else '+'

            genome_start, genome_stop = min(all_genome_positions), max(all_genome_positions)
            nuc_seq_fasta = genome_seq[mate1_chrom][genome_start-1:genome_stop]
            nuc_seq = nuc_seq_fasta.seq.upper()

            if (genome_strand == '-'):
                print >> sys.stderr, "On - genome strand"
                pdb.set_trace()
                nuc_seq = translate(nuc_seq, DNA_complement_table)[::-1]
                mate1_seq = translate(mate1_seq, DNA_complement_table)[::-1]
                mate1_qual = mate1_qual[::-1]
                mate2_seq = translate(mate2_seq, DNA_complement_table)[::-1]
                mate2_qual = mate2_qual[::-1]

            a_genome = NamedTemporaryFile(prefix='tmp', dir=None, delete=False)
            a_genome.write(">genome\n%s\n" % nuc_seq)
            a_genome.close()
            
            # 2) create pairwise alignment between mate1 and genome and profile-align to genome
            S = []
            for (query_pos, genome_pos) in genome_align1.get_aligned_pairs():
                q, t = '-', '-'
                if (query_pos != None):
                    q = mate1_seq[query_pos]
                if (genome_pos != None):
                    t = nuc_seq[genome_pos-genome_start]
                S.append( (q,t) )
            a2_in = NamedTemporaryFile(prefix='tmp', dir=None, delete=False)
            a2_out = NamedTemporaryFile(prefix='tmp', dir=None, delete=False)
            a2_in.write(">mate1\n%s\n>target1\n%s\n" % (''.join(map(itemgetter(0), S)), ''.join(map(itemgetter(1), S))))
            a2_in.close()
            
            cmd_line = "muscle -quiet -profile -in1 %s -in2 %s -out %s" % (a2_in.name, a_genome.name, a2_out.name)
            err_code = subprocess.check_call(cmd_line)
            #child = subprocess.Popen(cmd_line, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            #stdout_data, stderr_data = child.communicate()
            #msa2 = AlignIO.read(StringIO(stdout_data), "fasta")
            #msa2.sort()
            
            # 3) create pairwise alignment between mate2 and genome and profile-align to genome
            S = []
            for (query_pos, genome_pos) in genome_align1.get_aligned_pairs():
                q, t = '-', '-'
                if (query_pos != None):
                    q = mate1_seq[query_pos]
                if (genome_pos != None):
                    t = nuc_seq[genome_pos-genome_start]
                S.append( (q,t) )
            a3_in = NamedTemporaryFile(prefix='tmp', dir=None, delete=False)
            a3_out = NamedTemporaryFile(prefix='tmp', dir=None, delete=False)
            a3_in.write(">mate2\n%s\n>target2\n%s\n" % (''.join(map(itemgetter(0), S)), ''.join(map(itemgetter(1), S))))
            a3_in.close()

            cmd_line = "muscle -quiet -profile -in1 %s -in2 %s -out %s" % (a3_in.name, a_genome.name, a3_out.name)
            child = subprocess.Popen(cmd_line, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            stdout_data, stderr_data = child.communicate()
            msa3 = AlignIO.read(StringIO(stdout_data), "fasta")
            #msa3.sort()

            # 4) muscle profile_align alignments from 2) and 3)
            cmd_line = "muscle -quiet -profile -in1 %s -in2 %s" % (a2_out.name, a3_out.name)
            child = subprocess.Popen(cmd_line, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            stdout_data, stderr_data = child.communicate()
            msa = AlignIO.read(StringIO(stdout_data), "fasta")

            # 5) extract mate1, mate2, and one of the target
            msa.sort() # Places the alignments in the order mate1, mate2, genome, genome, target1, target2
            assert (msa[2].seq == msa[3].seq), "Same genome sequences don't have same alignment in MSA"
            msa = msa[0:3]

            unlink(a_genome)
            unlink(a2_in.name)
            unlink(a2_out.name)
            unlink(a3_in.name)
            unlink(a3_out.name)
            
            # 6) remove '---' columns ?

            try:
                merge_result = reconstructSequenceFromMSA(target.ID, curr_readID, msa, target_align1.qual, target_align2.qual, False)
            except AssertionError, ae:
                print >> sys.stderr, ae.message
                pdb.set_trace()
                
            if (merge_result.edit_dist < curr_edit_dist):
                new_cigar = "%s:%s:%d-%d:%s" % (mate1_chrom, genome_strand, genome_start, genome_stop, merge_result.cigar)
                merge_result = merge_result._replace(cigar=new_cigar)
                alt_genomic_alignments.append(merge_result)

        else:
            print >> sys.stderr, "WARNING: non-overlapping genomic alignment found for %s." % curr_readID
            edit_distance = -1    # TODO: this is nM tag value + #S,I,D. Be careful not to double-count the S,I,D
            if (candidate_merge != None and edit_dist < curr_edit_dist):
                S = "%s,%s" % (mate1_seq, mate2_seq)
                Q = "%s,%s" % (mate1_qual, mate2_qual)
                cigarstring = "%s:%s:%d-%d:%s,%s" % (mate1_chrom, genomic_target_strand, min(all_genome_positions)+1, max(all_genome_positions)+1,
                                                     genome_align1.cigarstring, genome_align2.cigarstring)
                candidate_merge = MergeMatesResult(targetID, curr_readID, S, Q, cigarstring, edit_distance)
                alt_genomic_alignments.append(candidate_merge)

    return alt_genomic_alignments


def altGenomicAlignments(targetID, curr_readID, target_align1, target_align2, curr_edit_dist, alt_alignments, genome_seq):
    ''' Nonstandard CIGAR format is <chr>:<strand>:<reference left pos>-<reference right pos>:<CIGAR>
        Sequence, quality string, and CIGAR are *NOT* reported in the 5'->3' orientation for those on the '-' strand.

        <reference start pos> and <reference stop pos> are in the context of the forward strand regardless of the actual strand of the genomic_target_strand

        Assumption: The alternative alignments (alt_alignments) provided to this function are to genomic regions that are not the source
                    of the expected target.
    '''

    global conv_cigar
    global DNA_complement_table
    alt_genomic_merges = []

    assert (len(alt_alignments[curr_readID]) <= 1), "ERROR: unhandled case of having multiple alt genomic alignments, for " % curr_readID
        
    for genome_align1, genome_align2 in alt_alignments[curr_readID]:
        cigar_strings, cigar_tuples = [], []
        for genome_align in [genome_align1, genome_align2]:
            if ('N' in genome_align.cigarstring):
                expanded_cigarstring = ''.join(map(lambda x: str(x[0]) * x[1], filter(lambda x: x[0]!=3, genome_align.cigartuples))) # Strip 'N' for intron skips
                cigartuples = ["%d%s" % (len(list(group)),conv_cigar[cigar_type]) for cigar_type, group in groupby(expanded_cigarstring)]
                cigar_strings.append( ''.join(cigartuples) )
                cigar_tuples.append(cigartuples)
            else:
                cigar_strings.append(genome_align.cigarstring)
                cigar_tuples.append(genome_align.cigartuples)

        if (target_align1.cigarstring != cigar_strings[0] or target_align2.cigarstring != cigar_strings[1]):
            # Evaluate based on number of exact match positions (ie number of 'M' in cigar minus basepair mismatches)
            target_align1_numM = sum(map(lambda y: y[1], filter(lambda x: x[0]==0, target_align1.cigartuples)))
            target_align2_numM = sum(map(lambda y: y[1], filter(lambda x: x[0]==0, target_align2.cigartuples)))
            genome_align1_numM = sum(map(lambda y: y[1], filter(lambda x: x[0]==0, genome_align1.cigartuples)))
            genome_align2_numM = sum(map(lambda y: y[1], filter(lambda x: x[0]==0, genome_align2.cigartuples)))

            target_align1_exact_match = target_align1_numM - target_align1.get_tag("NM")
            target_align2_exact_match = target_align2_numM - target_align2.get_tag("NM")
            genome_align1_exact_match = genome_align1_numM - genome_align1.get_tag("NM")
            genome_align2_exact_match = genome_align2_numM - genome_align2.get_tag("NM")

            target_align_exact_matches = target_align1_exact_match + target_align2_exact_match
            genome_align_exact_matches = genome_align1_exact_match + genome_align2_exact_match

            if (genome_align_exact_matches > target_align_exact_matches):
                # Need to return a pseudo-CIGAR that includes target (because using same primers) but also indicates a genomic
                # alignment (because the target is not the appropriate reference)
                #mate1_chrom = alt_alignments['references'][genome_align1.tid]
                #mate2_chrom = alt_alignments['references'][genome_align2.tid]
                mate1_seq, mate1_qual = genome_align1.query_sequence, genome_align1.qual
                mate2_seq, mate2_qual = genome_align2.query_sequence, genome_align2.qual
                assert (mate1_chrom == mate2_chrom), "ERROR: unhandled case of better genomic alignment being interchromosomal in altGenomicAlignments()"

                # Merge non-overlapping reads with intervening sequence if gap has no known intron (or just a splice junction?) and
                #      1) not too long (limit?), or
                #      2) within (potentially long) exon
                # Otherwise, join sequences with a '+'.

                genome_align1_positions = set(genome_align1.positions)
                genome_align2_positions = set(genome_align2.positions)
                pairwise_alignments_overlap = len(genome_align1_positions & genome_align2_positions) > 0
                all_genome_positions = genome_align1_positions.union( genome_align2_positions )

                genomic_target_strand = '-' if (genome_align1.is_reverse) else '+'

                # TODO: change to using the new Fasta() method
                #region_spec = "%s:%d-%d" % (mate1_chrom, min(all_genome_positions)+1, max(all_genome_positions)+1)
                #nuc_seq_fasta = pysam.faidx(genome_seq, region_spec)
                #genomic_target_seq = ''.join(map(lambda x: x.strip().upper(), nuc_seq_fasta[1:]))

                #if (genomic_target_strand == '-'):
                #    genomic_target_seq = translate(genomic_target_seq, DNA_complement_table)[::-1]
                #    mate1_seq = translate(mate1_seq, DNA_complement_table)[::-1]
                #    mate1_qual = translate(mate1_qual, DNA_complement_table)[::-1]
                #    mate2_seq = translate(mate2_seq, DNA_complement_table)[::-1]
                #    mate2_qual = translate(mate2_qual, DNA_complement_table)[::-1]
                    
                if (pairwise_alignments_overlap):
                    # 1) Extract genomic sequence based on min/max alignment positions
                    # 2) create pairwise alignment between mate1 and its local genomic
                    # 3) muscle profile_align 1)&2)
                    # 4) Do same as 2) & 3) for mate2
                    # 5) muscle profile_align alignments from 3) and 4)
                    # 6) extract mate1, mate2, and target (the sequence from 1).
                    # 7) remove '---' columns
                    # 8) call reconstruct...() to get S, Q, C, edit_dist

                    S, Q, cigarstring, edit_dist = mergeOverlappingPairwiseAlignments(genome_align1, genome_align2)

                    if (edit_dist < curr_edit_dist):
                        new_cigar = "%s:%s:%d-%d:%s" % (mate1_chrom, genomic_target_strand, min(all_genome_positions)+1, max(all_genome_positions)+1, cigarstring)
                        candidate_merge = MergeMatesResult(targetID, curr_readID, S, Q, cigarstring, edit_dist)
                        alt_genomic_merges.append(candidate_merge)
                        #bed12_line = createBed12Line(mate1_chrom, genomic_target_strand, list(all_genome_positions))
                else:
                    print >> sys.stderr, "WARNING: non-overlapping genomic alignment found for %s." % curr_readID
                    edit_distance = -1    # TODO: this is nM tag value + #S,I,D. Be careful not to double-count the S,I,D
                    if (candidate_merge != None and edit_dist < curr_edit_dist):
                        S = "%s,%s" % (mate1_seq, mate2_seq)
                        Q = "%s,%s" % (mate1_qual, mate2_qual)
                        cigarstring = "%s:%s:%d-%d:%s,%s" % (mate1_chrom, genomic_target_strand, min(all_genome_positions)+1, max(all_genome_positions)+1,
                                                             genome_align1.cigarstring, genome_align2.cigarstring)
                        candidate_merge = MergeMatesResult(targetID, curr_readID, S, Q, cigarstring, edit_distance)
                        alt_genomic_merges.append(candidate_merge)

    return alt_genomic_merges


def resolveIncompleteAlignmentsToTarget(curr_readID, target_align1, target_align2, target, incompl_merged_to_target, genome_seq):
    '''Merges the two mates into one sequence as best it can using target.'''
    merge_result = None
    
    assert (not target_align1.is_reverse and target_align2.is_reverse), "Mates 1/2 don't align in F/R orientation to target sequence." % curr_readID

    target_align1_positions = set(target_align1.positions)
    target_align2_positions = set(target_align2.positions)
    pairwise_alignments_overlap = not target_align1_positions.isdisjoint(target_align2_positions)
    
    # 1) Compare incompl_merged_to_target alignment to target_align1/2.
    #    - If merged has same alignment positions as target_align1/2, confirm/correct nucleotide sequence in merged based on quals. Use merged's CIGAR
    #           - Flag cases with soft-clip. Investigate.
    #           - This case could also have indels. May need todo alignment of mates back to merged.
    #    - If merged has more matches, align with mates and the target. Check that merged sequence handles indels and nucleotide/quals properly.
    #           - Flag cases with soft-clip. Investigate.
    #           - If muscle doesn't properly align all four sequences:
    #                - Recreate merged to target alignment.
    #                - Align mates to merged sequences.
    #                - Use muscle to align the two profiles.
    #    - If merged has fewer matches, investigate.
    if (incompl_merged_to_target.has_key( (curr_readID, target.ID) )):
        merged_align = incompl_merged_to_target[(curr_readID, target.ID)]
        assert (len(merged_align) == 1), "Multiple alignments of merged sequence to target sequence for %s" % curr_readID
        merged_align = merged_align[0]
        #merged_align_positions = set(merged_align.positions)
        #target_align_positions = target_align1_positions.union(target_align2_positions)

        merged_num_match_position, merged_num_mismatches = 0,0
        for q,t in filter(lambda x: x[0]!='-' and x[1]!='-', zip(merged_align[0], merged_align[1])):
            if (q != t):
                merged_num_mismatches += 1
            merged_num_match_position += 1

        mismatch_frac = float(merged_num_mismatches)/float(merged_num_match_position)
        
        mate1_seq = target_align1.query_sequence
        mate2_seq = target_align2.query_sequence
        target_seq = target.seq

        #merged_is_CandC = ('-' not in merged_align[0] and '-' not in merged_align[1] )
        #if (len(target_align1.cigartuples)==1 and len(target_align2.cigartuples)==1 and merged_is_CandC): # len(merged_align.cigartuples)==1
        #    #assert (merged_align_positions == target_align_positions), "Alignment of merged mates sequence and individual mates to target not the same, for %s" % curr_readID
        #    fasta_lines = ''
        #    for name, mate_seq, pw_alignment in [("mate1",mate1_seq,target_align1), ("mate2",mate2_seq,target_align2)]: #, ("target",merged_align)]:
        #        S = ['-'] * len(target_seq)
        #        for (query_pos, target_pos) in pw_alignment.get_aligned_pairs():
        #            S[query_pos] = mate_seq[query_pos]
        #        fasta_lines += ">%s\n%s\n" % (name, ''.join(S))
        #    fasta_lines += ">target\n%s\n" % target_seq
        #    pdb.set_trace()
        #    msa = AlignIO.read(StringIO(fasta_lines), "fasta")
        #else:
        if (mismatch_frac <= 0.1): # Do I need the above piece? What was the thinking behind it?
            # How to handle: Align mates to merged and correct merged.
            # If merged sequence was corrected, realign merged back to target with usearch to get new cigar, etc. Not implemented.
            # OR
            merged_seq = merged_align[0].replace('-','') #merged_align.query_sequence
            # 1) Recreate merged to target alignment.
            #merged_align_str, target_align_str = '',''
            #for query_pos, target_pos in merged_align.get_aligned_pairs():
            #    merged_align_str += '-' if (query_pos == None) else merged_seq[query_pos]
            #    target_align_str += '-' if (target_pos == None) else target_seq[target_pos]
            a1 = NamedTemporaryFile(prefix='tmp', dir=None, delete=False)
            a1.write(">%s\n%s\n>%s\n%s\n" % ("Merged1", merged_align[0], "target", merged_align[1]))
            a1.close()
            
            # 2) muscle align the mates to the merged sequences
            three_seqs = ">Merged2\n%s\n>mate1\n%s\n>mate2\n%s\n" % (merged_seq, mate1_seq, mate2_seq)
            child = subprocess.Popen("muscle -quiet", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            child.stdin.write(three_seqs)
            child.stdin.close()
            align = AlignIO.read(child.stdout, "fasta")
            align.sort() # Plased the alignments in the order Merged2, mate1, mate2
            
            # Check the alignment of the mates back to the merged sequence that was formed from them. If okay, proceed.
            if (not any(map(lambda c: align[1:3,c]=="--", range(align.get_alignment_length())))):
                a2 = NamedTemporaryFile(prefix='tmp', dir=None, delete=False)
                a2.write(">%s\n%s\n>%s\n%s\n>%s\n%s\n" % (align[0].name, align[0].seq, align[1].name, align[1].seq, align[2].name, align[2].seq))
                a2.close()

                # 3) muscle profile_align the merged+target msa to the merged+mates msa
                cmd_line = "muscle -quiet -profile -in1 %s -in2 %s" % (a1.name, a2.name)
                child = subprocess.Popen(cmd_line, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                stdout_data, stderr_data = child.communicate()
                msa = AlignIO.read(StringIO(stdout_data), "fasta")
                unlink(a2.name)
                
                # 4) remove all but mate1, mate2, and target and reconstruct S, Q, C, and edit_distance as already implemented.
                msa.sort() # Places the alignments in the order Merged1, Merged2, mate1, mate2, target
                msa = msa[2:]

                merge_result = reconstructSequenceFromMSA(target.ID, curr_readID, msa, target_align1.qual, target_align2.qual, True)

            unlink(a1.name)


    elif (pairwise_alignments_overlap):
        # This case should just be for mates from very long targets for which mates wouldn't overlap.
        print "WARNING: Found case where mate alignments to target overlap, but they could not be merged: %s" % curr_readID
        # How to handle. parameterize PEAR with ~number of overlap positions and/or product length, and with very loose P-value.
        #pdb.set_trace()
    else:
        # Confirm that mates should not overlap, given their combined length and length of the target.
        print >> sys.stderr, "WARNING: mates for %s did not overlap and were not merged" % curr_readID
        #pdb.set_trace()

    return merge_result


#    # For the two mates, with soft-clipped positions removed
#    target_align1_seq, target_align1_qual = "", ""
#    target_align2_seq, target_align2_qual = "", ""
#
#    # Attempt merging using pairwise alignments from BWA
#    # For Mate 1
#    align_w_cigar1 = zip(target_align1.get_aligned_pairs(), list(''.join(map(lambda x: conv_cigar[str(x[0])] * x[1], target_align1.cigartuples))))
#    alignment1_query = ['-'] * len(align_w_cigar1)
#    alignment1_qual = [None] * len(align_w_cigar1)
#    alignment1_targ = ['-'] * len(align_w_cigar1)
#
#    for c, ((query_pos, target_pos), cigar_letter) in enumerate(align_w_cigar1):
#        if (query_pos!=None):
#            alignment1_query[c] = target_align1.query_sequence[query_pos]
#            alignment1_qual[c] = target_align1.qual[query_pos]
#        if (target_pos!=None):
#            alignment1_targ[c] = target.seq[target_pos]
#
#    assert (target_align1.cigartuples[0][0] in [0,4]), "Mate1 of %s, alignment problem to 5p"
#    #trim_lindex = target_align1.cigartuples[0][1] if (target_align1.cigartuples[0][0] == 4) else 0
#    # In the upstream standardization step mate 1 was deemed to align to, at least, the 5p---even if it had mismatches. Those mismatches could be
#    # soft-clipped here, which I want to override.
#    trim_lindex = 0
#    trim_rindex = len(align_w_cigar1) - target_align1.cigartuples[-1][1] if (target_align1.cigartuples[-1][0] == 4) else len(align_w_cigar1)
#
#    # Alignments with soft-clipped trimmed
#    alignment1_query = ''.join(alignment1_query[trim_lindex:trim_rindex])
#    alignment1_qual = alignment1_qual[trim_lindex:trim_rindex]
#    alignment1_targ = ''.join(alignment1_targ[trim_lindex:trim_rindex])
#
#    # Soft-clip trimmed sequences
#    target_align1_seq = alignment1_query.replace('-','')
#    target_align1_qual = filter(lambda x: x!=None, alignment1_qual)
#    target_align1_target = alignment1_targ.replace('-','')
#
#
#    # For Mate 2
#    align_w_cigar2 = zip(target_align2.get_aligned_pairs(), list(''.join(map(lambda x: conv_cigar[str(x[0])] * x[1], target_align2.cigartuples))))
#    alignment2_query = ['-'] * len(align_w_cigar2)
#    alignment2_qual = [None] * len(align_w_cigar2)
#    alignment2_targ = ['-'] * len(align_w_cigar2)
#
#    for c, ((query_pos, target_pos), cigar_letter) in enumerate(align_w_cigar2):
#        if (query_pos!=None):
#            alignment2_query[c] = target_align2.query_sequence[query_pos]
#            alignment2_qual[c] = target_align2.qual[query_pos]
#        if (target_pos!=None):
#            alignment2_targ[c] = target.seq[target_pos]
#
#    trim_lindex = target_align2.cigartuples[0][1] if (target_align2.cigartuples[0][0] == 4) else 0
#    assert (target_align2.cigartuples[0][0] in [0,4]), "Mate2 of %s, alignment problem to 3p"
#    #trim_rindex = len(align_w_cigar2) - target_align2.cigartuples[-1][1] if (target_align2.cigartuples[-1][0] == 4) else len(align_w_cigar2)
#    # In the upstream standardization step mate 2 was deemed to align to, at least, the 3p---even if it had mismatches. Those mismatches could be
#    # soft-clipped here, which I want to override.
#    trim_rindex = len(align_w_cigar2)
#
#    # Alignments with soft-clipped trimmed
#    alignment2_query = ''.join(alignment2_query[trim_lindex:trim_rindex])
#    alignment2_qual = alignment2_qual[trim_lindex:trim_rindex]
#    alignment2_targ = ''.join(alignment2_targ[trim_lindex:trim_rindex])
#
#    # Soft-clip trimmed sequences
#    target_align2_seq = alignment2_query.replace('-','')
#    target_align2_qual = filter(lambda x:x!=None, alignment2_qual)
#    target_align2_target = alignment2_targ.replace('-','')


#    alignments1 = [(alignment1_query, alignment1_targ, None, None, None)]
#    alignments2 = [(alignment2_query, alignment2_targ, None, None, None)]
#
#    olap_indel_positions1, olap_indel_positions2 = set(), set()
#    if ((len(alignments1[0][0]) == len(alignments2[0][0]))):
#        olap_cols = map(lambda w: w[0], filter(lambda y: y[1][0]!='-' and y[1][1]!='-', enumerate(zip(alignments1[0][0], alignments2[0][0]))))
#            
#        if (len(olap_cols) > 0):
#            olap_start, olap_stop = (olap_cols[0], olap_cols[-1]) 
#            olap_indel_positions1.update( [pos for pos, char in enumerate(alignments1[0][0][olap_start:olap_stop+1]) if char == '-'] )
#            olap_indel_positions1.update( [pos for pos, char in enumerate(alignments1[0][1][olap_start:olap_stop+1]) if char == '-'] )
#            olap_indel_positions2.update( [pos for pos, char in enumerate(alignments2[0][0][olap_start:olap_stop+1]) if char == '-'] )
#            olap_indel_positions2.update( [pos for pos, char in enumerate(alignments2[0][1][olap_start:olap_stop+1]) if char == '-'] )
#
#    # If the pairwise alignments differ in lengths or on the columns where indels occur, align the mate sequences and the target sequence altogether
#    if ((len(alignments1[0][0]) != len(alignments2[0][0])) or (len(olap_indel_positions1 ^ olap_indel_positions2)>0)):
#        method_used = "MSA"
#        three_seqs = ">target\n%s\n>mate1\n%s\n>mate2\n%s" % (target.seq, target_align1_seq, target_align2_seq)
#        child = subprocess.Popen("muscle -quiet", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#        child.stdin.write(three_seqs)
#        child.stdin.close()
#        align = AlignIO.read(child.stdout, "fasta")
#        align.sort() # Places the alignments in the order mate1, mate2, target
#        try:
#            merge_result = reconstructSequenceFromMSA(target.ID, curr_readID, align, target_align1_qual, target_align2_qual, True)
#        except AssertionError:
#            pdb.set_trace()
#    else:
#        method_used = "PW"
#        assert (alignments1[0][1] == alignments2[0][1]) 
#        x = np.array([alignments1[0][0], alignments2[0][0], alignments1[0][1]], dtype=str)
#        align = x.view('S1').reshape((x.size,-1))
#        try:
#            merge_result = reconstructSequenceFromMSA(target.ID, curr_readID, align, alignment1_qual, alignment2_qual, True)
#        except AssertionError:
#            pdb.set_trace()
#    
#    if (False and merge_result != None and method_used == "MSA"): # 
#        should_be_merged = len(target_align1_seq) + len(target_align2_seq) > len(target.seq)
#        print >> sys.stderr, "--------------------------------------------------------"
#        print >> sys.stderr, "Should be merged =", should_be_merged
#        print >> sys.stderr, "Method used       =", method_used
#        print >> sys.stderr, "### MSA"
#        print >> sys.stderr, str(align[0].seq), "mate1", target_align1.cigarstring
#        print >> sys.stderr, str(align[1].seq), "mate2", target_align2.cigarstring
#        print >> sys.stderr, str(align[2].seq), "target"
#        print >> sys.stderr, merge_result.seq, "manual_merge_seq"
#
#        print >> sys.stderr, "\n### PW"
#        alignments1 = pairwise2.align.localms(target_align1_seq, target.seq, 1,-2,-6,-1)
#        print >> sys.stderr, alignments1[0][0], "mate1", target_align1.cigarstring
#        print >> sys.stderr, alignments1[0][1], "target"
#        alignments2 = pairwise2.align.localms(target_align2_seq, target.seq, 1,-2,-6,-1)
#        print >> sys.stderr, alignments2[0][0], "mate2", target_align2.cigarstring
#        print >> sys.stderr, alignments2[0][1], "target"
#        print >> sys.stderr, "\n"
#    
#    elif (False and merge_result != None):
#        should_be_merged = len(target_align1_seq) + len(target_align2_seq) > len(target.seq)
#        print >> sys.stderr, "Should be merged =", should_be_merged
#        print >> sys.stderr, "Method used      =", method_used
#        print >> sys.stderr, "### PW"
#        alignments1 = pairwise2.align.localms(target_align1_seq, target.seq, 1,-2,-6,-1)
#        print >> sys.stderr, alignments1[0][0], "mate1", target_align1.cigarstring
#        print >> sys.stderr, alignments1[0][1], "target"
#        alignments2 = pairwise2.align.localms(target_align2_seq, target.seq, 1,-2,-6,-1)
#        print >> sys.stderr, alignments2[0][0], "mate2", target_align2.cigarstring
#        print >> sys.stderr, alignments2[0][1], "target"
#        print >> sys.stderr, merge_result.seq, "manual_merge_seq"
#        print >> sys.stderr, "\n"
#
#    if (merge_result == None):
#        # TODO: implement tryHarder() to try more involved procedures to merge the mates.
#        #    Try both usearch -merge_fastq (maybe w/ lenient olap restriction), usearch_chimera,
#        #    muscle -profile to aligh the the mate1+target and mate2+target pairwise alignments, etc
#        print >> sys.stderr, "WARNING: could not merge mates for %s" % curr_readID
#        #pdb.set_trace()
#
#    return merge_result


def performReadPairQC(curr_readID, target_seq_data, bam_handle, target_align1, target_align2,
                      alt_compl_alignments, compl_merged_to_genome,
                      alt_incompl_alignments, incompl_merged_to_genome):

    '''The first standardized and trimmed read aligns, at a minimum, to the
    beginning primer sequence of the target and the second read to the end
    primer sequence of the target. That was a requirement to their inclusion in the
    standardized reads input file.'''

    one_mate_unaligned = target_align1.is_unmapped ^ target_align2.is_unmapped
    both_mates_unaligned_to_target = target_align1.is_unmapped and target_align2.is_unmapped
    mates_aligned_different_targets = target_align1.tid != target_align2.tid

    if (mates_aligned_different_targets or both_mates_unaligned_to_target):
        target_seq_ID = None
        mates_not_aligned_to_primers = True
    else:
        target_seq_ID = bam_handle.references[target_align1.tid]
        target_seq_data_tuple = target_seq_data[target_seq_ID]

        fwd_primer_len = target_seq_data_tuple[5]
        rev_primer_len = target_seq_data_tuple[6]

        # Confirm that each read at least matches its primer 
        mates_not_aligned_to_primers = not((target_align1.cigartuples[0][0]==0 and target_align1.cigartuples[0][1] >= fwd_primer_len) and \
          (target_align2.cigartuples[-1][0]==0 and target_align2.cigartuples[-1][1] >= rev_primer_len))
    
    # Do both reads align Completely and Contiguously?
    both_aligned_CandC = (not mates_aligned_different_targets and
                          len(target_align1.cigartuples)==1 and target_align1.cigartuples[0][0]==0 and 
                          len(target_align2.cigartuples)==1 and target_align2.cigartuples[0][0]==0)

    target_align_has_internal_softclipping = target_align1.cigartuples[-1][-1]==4 or target_align2.cigartuples[0][0]==4
    has_merged_genome_aln = incompl_merged_to_genome.has_key(curr_readID) or compl_merged_to_genome.has_key(curr_readID)
    has_paired_genome_aln = alt_incompl_alignments.has_key(curr_readID) or incompl_merged_to_genome.has_key(curr_readID)
    mates_to_target_align_clipped_and_unsupported = target_align_has_internal_softclipping and not (has_merged_genome_aln or has_paired_genome_aln)
        
    standardized_reads_aligned_ok = not(both_mates_unaligned_to_target or mates_aligned_different_targets or
                                        mates_not_aligned_to_primers or mates_to_target_align_clipped_and_unsupported)

    qc_report = ReadPairQC(both_mates_unaligned_to_target, mates_aligned_different_targets, one_mate_unaligned,
                           mates_not_aligned_to_primers, both_aligned_CandC, mates_to_target_align_clipped_and_unsupported,
                           standardized_reads_aligned_ok)

    return (target_seq_ID, qc_report)


def judgeBestMerge(target_based_merge, alt_genomic_merges):
    '''Genome alignment better if it has smaller edit distance.'''
    candidates = [] # tuples of (merge result, edit_dist, is alt genomic)
    
    if (target_based_merge != None):
        candidates.append( (target_based_merge, target_based_merge.edit_dist, False) )

    for alt in alt_genomic_merges:
        candidates.append( (alt, alt.edit_dist, True) )

    if (len(candidates) > 0):
        candidates = sorted(candidates, key=itemgetter(1))
        smallest_edit_dist = candidates[0][1]
        best_merges = filter(lambda x: x[1] <= smallest_edit_dist, candidates)

        if (len(best_merges)>1):
            if (any(map(lambda x: x[2], best_merges)) and any(map(lambda x: not x[2], best_merges))):
                pass
            elif (all(map(lambda x: not x[2], best_merges))):
                print >> sys.stderr, "WARNING: enountered case of multiple equally good merges to target. Shouldn't happen. %s" % target_based_merge.matesID
                pdb.set_trace()
    else:
        best_merges = []

    return best_merges


def mergeStandardizedReads(BWA_bam, target_seq_data, alt_incompl_alignments, alt_compl_alignments,
                           compl_merged_to_genome, incompl_merged_to_genome, incompl_merged_to_target,
                           genome_seq, clustered_seqs_output):

    filtered_out_counts = defaultdict(int)
    merge_stats = ProcessingStats()

    ip_bam = pysam.AlignmentFile(BWA_bam, 'rb')
    op = gzip.open(clustered_seqs_output, 'wb')
    #op_debug = open("align_diff_targets.txt", 'w')
    
    try:
        prev_readID = None
        #counter = 0
        while (True):
            #counter += 1
            #if (counter == 100):
            #    raise StopIteration

            target_align1 = ip_bam.next()
            while (target_align1.is_secondary or target_align1.is_supplementary):
                target_align1 = ip_bam.next()

            target_align2 = ip_bam.next()
            while (target_align2.is_secondary or target_align2.is_supplementary):
                target_align2 = ip_bam.next()

            assert (target_align1.qname == target_align2.qname), "Names differ for the two mates: %s and %s" % (target_align1.qname, target_align2.qname)
            curr_readID = target_align1.qname
            
            # Expecting just one target match per read pair, and in standardized ordering
            assert (curr_readID != prev_readID), "In mergeStandardizedReads(), curr_readID == prev_readID == %s" % curr_readID
            assert (target_align1.is_read1 and target_align2.is_read2), "SAM line pair not in standardized order for %s" % curr_readID

            target_seq_ID, qc_report = performReadPairQC(curr_readID, target_seq_data, ip_bam, target_align1, target_align2,
                                                         alt_compl_alignments, compl_merged_to_genome,
                                                         alt_incompl_alignments, incompl_merged_to_genome)

            target_based_merge, alt_genomic_merges = None, []
            if (qc_report.standardized_reads_aligned_ok):
                target = MergeTarget(target_seq_ID, target_seq_data[target_seq_ID][0])

                if (qc_report.both_aligned_CandC):
                    target_based_merge = mergeCompleteTargetAlignments(curr_readID, target_align1, target_align2, target, genome_seq)
                    target_edit_dist = target_based_merge.edit_dist if (target_based_merge != None) else 1e10 
                    if (alt_compl_alignments.has_key(curr_readID) or compl_merged_to_genome.has_key(curr_readID)):
                        alt_genomic_merges = altGenomicAlignments2(target.ID, curr_readID, target_align1, target_align2, target_edit_dist,
                                                                   alt_compl_alignments, compl_merged_to_genome, genome_seq)
                else:
                    target_based_merge = resolveIncompleteAlignmentsToTarget(curr_readID, target_align1, target_align2, target, incompl_merged_to_target, genome_seq)
                    target_edit_dist = target_based_merge.edit_dist if (target_based_merge != None) else 1e10
                    if (alt_incompl_alignments.has_key(curr_readID) or incompl_merged_to_genome.has_key(curr_readID)):
                        alt_genomic_merges = altGenomicAlignments2(target.ID, curr_readID, target_align1, target_align2, target_edit_dist,
                                                                   alt_incompl_alignments, incompl_merged_to_genome, genome_seq)

            best_merges = judgeBestMerge(target_based_merge, alt_genomic_merges)
            if (len(best_merges) != 0):
                writeSequenceOutput(op, best_merges)
                    
            # merge_result != None, best_is_genomic
            merge_stats.updateForMergeAttempt(target_seq_ID, best_merges, qc_report)

            prev_readID = curr_readID
            target_align1, target_align2 = None, None

    except StopIteration:
        # This works because the input file should only have primary and properly-paired reads
        if (target_align1 != None or target_align2 != None):
            pdb.set_trace()
        assert (target_align1 == target_align2 == None), "ERROR: input parsing completed with only one mate read"
    except AssertionError, ae:
        print >> sys.stderr, "Assertion Error:", ae.message
        print >> sys.stderr, "Current read ID is %s" % curr_readID
    finally:
        ip_bam.close()
        op.close()
        #op_debug.close()
        
    return merge_stats


def writeSequenceOutput(op, merge_results):
    for merge_result in map(itemgetter(0), merge_results):
        readID, target_and_umi = merge_result.matesID.split(',',1)
        target_part, umi_part = target_and_umi.rsplit('_',1)
        try:
            assert (target_part == merge_result.targetID), "ERROR: mate's for %s not matched to expected target %s" % (merge_result.matesID, merge_result.targetID)
        except AssertionError, ae:
            print >> sys.stderr, ae.message
            #pdb.set_trace()

        op.write("%s\t%s\t%s\t%d\t%s\n" % (target_part, umi_part, merge_result.cigar, merge_result.edit_dist, merge_result.seq))
    

if (__name__ == "__main__"):
    genome_fasta, BWA_bam, incompl_HISAT_bam, incompl_STAR_bam, compl_HISAT_bam, compl_STAR_bam, \
    compl_merged_to_genome_HISAT_bam, compl_merged_to_genome_STAR_bam, \
    incompl_merged_to_genome_HISAT_bam, incompl_merged_to_genome_STAR_bam, \
    incompl_merged_to_target_fasta, \
    primers_fasta, products_fasta, output_log, clustered_seqs_output, failed_merge_counts = sys.argv[1:]

    genome_seq = Fasta(genome_fasta)
    target_seq_data = parseTargetSeqData(primers_fasta, products_fasta)

    print >> sys.stderr, "INFO: reading complete and incomplete merged-sequence genome/target alignments:"
    print >> sys.stderr, "\t", compl_merged_to_genome_HISAT_bam, compl_merged_to_genome_STAR_bam
    compl_merged_to_genome = readSingletonAlignmentsBam([compl_merged_to_genome_HISAT_bam, compl_merged_to_genome_STAR_bam])
    print >> sys.stderr, "\t", incompl_merged_to_genome_HISAT_bam, incompl_merged_to_genome_STAR_bam
    incompl_merged_to_genome = readSingletonAlignmentsBam([incompl_merged_to_genome_HISAT_bam, incompl_merged_to_genome_STAR_bam])
    print >> sys.stderr, "\t", incompl_merged_to_target_fasta
    incompl_merged_to_target = readSingletonAlignmentsFasta([incompl_merged_to_target_fasta])
    
    print >> sys.stderr, "INFO: reading complete and incomplete paired-end genome alignments:"
    print >> sys.stderr, "\t", incompl_HISAT_bam, incompl_STAR_bam
    alt_incompl_alignments = readGenomicReadPairAlignments([incompl_HISAT_bam, incompl_STAR_bam], target_seq_data)
    print >> sys.stderr, "\t", compl_HISAT_bam, compl_STAR_bam
    alt_compl_alignments = readGenomicReadPairAlignments([compl_HISAT_bam, compl_STAR_bam], target_seq_data)

    print >> sys.stderr, "INFO: merging read pairs"
    merge_stats = mergeStandardizedReads(BWA_bam, target_seq_data, alt_incompl_alignments, alt_compl_alignments,
                                         compl_merged_to_genome, incompl_merged_to_genome, incompl_merged_to_target,
                                         genome_seq, clustered_seqs_output)
    merge_stats.checkAndReport(output_log)

    print >> sys.stderr, "INFO: Done"
    sys.exit(0)
