#!/usr/bin/env python

from operator import itemgetter
from collections import Counter, defaultdict
import sys


def readSubpools(pp_subpools_file):
    subpool_per_pp = {}
    with open(pp_subpools_file, 'r') as ip:
        for line in ip:
            subpool, pp_name = line.strip().split("\t")
            subpool_per_pp[pp_name] = subpool
    return subpool_per_pp

    
def readNamesExpectedPrimerPairs(pp_names_file):
    pp_full_names = set()

    with open(pp_names_file, 'r') as ip:
        for line in ip:
            pp_name = line.strip()
            pp_full_names.add(pp_name)

    return pp_full_names            


def standardizePrimerName(full_primer_name):
    if ("CTRL" in full_primer_name or full_primer_name.startswith("UDTD") or full_primer_name.startswith("ERCC")):
        standardized_name = full_primer_name
    else:
        elems = full_primer_name.split('_')

        F, R = elems[3].split('R')

        R = "R" + R
        if (elems[-1] == "5p"):
            standardized_name = "%s_%s_%s_%s" % (elems[0], elems[1], elems[2], F)
        else:
            standardized_name = "%s_%s_%s_%s" % (elems[0], elems[1], elems[2], R)

    return standardized_name



def processLinesForCurrRead(curr_read, curr_read_lines, pp_full_names):
    processed_read_pair = None
    assert (len(curr_read_lines)>0)

    R1_seq, R2_seq = None, None
    R1_start, R1_stop = None, None
    R2_start, R2_stop = None, None
    R1_sense_UMI, R2_sense_UMI = 'None', 'None'
    R1_sense_primer_name, R1_antisense_primer_name = 'None', 'None'
    R2_sense_primer_name, R2_antisense_primer_name = 'None', 'None'
    R1_primer_signature, R2_primer_signature = 'None', 'None'
        
    has_swapped_sense_primers = False

    R1_lines = filter(lambda x:x[-1]=="R1", curr_read_lines)
    sense_primer_src = antisense_primer_src = None
    sense_equiv_region = antisense_equiv_region = None
    R1_sense_primer_name, R1_antisense_primer_name = 'None', 'None'
    R1_multiple_starts = set()
    for R1_line_elems in R1_lines:
        if (len(R1_line_elems)==2):
            assert (R1_seq == None), "R1_seq already set"
            R1_seq = R1_line_elems[0]

        elif (R1_line_elems[0] == '+'):
            if (R1_start != None):
                R1_multiple_starts.add(curr_read)
                #print("%s has multiple R1 starts" % curr_read)
                    
            R1_sense_primer_name = R1_line_elems[4]
            if (R1_sense_primer_name.startswith('UDTD')):
                R1_start = int(R2_line_elems[2])
            else:
                sense_primer_src, sense_equiv_region = R1_line_elems[4].rsplit('_',1)
                if (sense_equiv_region[0] != 'F'):
                    has_swapped_sense_primers = True
                    #print("%s has R1 with sense rev primer" % curr_read)

                if (R1_start == None or int(R1_line_elems[1]) < R1_start):
                    R1_start = int(R1_line_elems[1])

            if (R1_stop == None):
                R1_stop = int(R1_line_elems[3])

        elif (R1_line_elems[0] == '-' and (R1_antisense_primer_name == "None" or int(R1_line_elems[1]) < R1_stop)):
            R1_antisense_primer_name = R1_line_elems[4]
            if (R1_antisense_primer_name.startswith('UDTD')):
                R1_stop = int(R1_line_elems[1])
                #antisense_primer_src, antisense_equiv_region = (R1_antisense_primer_name, R1_antisense_primer_name)
            else:
                antisense_primer_src, antisense_equiv_region = R1_line_elems[4].rsplit('_',1)
                R1_stop = int(R1_line_elems[2])
                
    if (R1_start != None):
        R1_sense_UMI = R1_seq[0:R1_start]
        R1_seq = R1_seq[R1_start-1:R1_stop]


    if (sense_primer_src != None and antisense_primer_src == None):
        R1_primer_signature = (sense_primer_src, sense_equiv_region, None, None)
    elif (sense_primer_src == None and antisense_primer_src != None):
        R1_primer_signature = (None, None, antisense_primer_src, antisense_equiv_region)
    elif (R1_sense_primer_name == R1_antisense_primer_name == 'None'):
        R1_primer_signature = (None, None, None, None)
    else:
        R1_primer_signature = (sense_primer_src, sense_equiv_region, antisense_primer_src, antisense_equiv_region)


    R2_lines = filter(lambda x:x[-1]=="R2", curr_read_lines)
    sense_primer_src = antisense_primer_src = None
    sense_equiv_region = antisense_equiv_region = None
    R2_antisense_primer_name, R2_sense_primer_name = 'None', 'None'
    R2_multiple_starts = set()
    for R2_line_elems in R2_lines:
        if (len(R2_line_elems)==2):
            assert (R2_seq == None), "R2_seq already set"
            R2_seq = R2_line_elems[0]

        elif (R2_line_elems[0] == '+'):
            if (R2_start != None):
                R2_multiple_starts.add(curr_read)
                #print("%s has multiple R2 starts" % curr_read)

            R2_sense_primer_name = R2_line_elems[4]
            if (R2_sense_primer_name.startswith('UDTD')):
                R2_start = int(R2_line_elems[2])
            else:
                sense_primer_src, sense_equiv_region = R2_line_elems[4].rsplit('_',1)
                if (sense_equiv_region[0] != 'R'):
                    has_swapped_sense_primers = True
                    #print("%s has R2 with sense fwd primer" % curr_read)

                if (R2_start == None or int(R2_line_elems[1]) < R2_start):
                    R2_start = int(R2_line_elems[1])

            if (R2_stop == None):
                R2_stop = int(R2_line_elems[3])

        elif (R2_line_elems[0] == '-' and (R2_antisense_primer_name == "None" or int(R2_line_elems[1]) < R2_stop)):
            R2_antisense_primer_name = R2_line_elems[4]
            if (R2_antisense_primer_name.startswith('UDTD')):
                #antisense_primer_src, antisense_equiv_region = (R2_antisense_primer_name, R2_antisense_primer_name)
                R2_stop = int(R2_line_elems[1])
            else:
                antisense_primer_src, antisense_equiv_region = R2_line_elems[4].rsplit('_',1)
                R2_stop = int(R2_line_elems[2])

    if (R2_start != None):
        R2_sense_UMI = R2_seq[0:R2_start]
        R2_seq = R2_seq[R2_start-1:R2_stop]

    if (sense_primer_src != None and antisense_primer_src == None):
        R2_primer_signature = (None, None, sense_primer_src, sense_equiv_region)
    elif (sense_primer_src == None and antisense_primer_src != None):
        R2_primer_signature = (antisense_primer_src, antisense_equiv_region, None, None)
    elif (R2_antisense_primer_name == R2_sense_primer_name == 'None'):
        R2_primer_signature = (None, None, None, None)
    else:
        R2_primer_signature = (antisense_primer_src, antisense_equiv_region, sense_primer_src, sense_equiv_region)


    primer_signature, triage_result = makePrimerSignatureForReadPair(R1_primer_signature, R2_primer_signature, pp_full_names)
    #print("%s\t%s" % (primer_signature, triage_result))
    processed_read_pair = (primer_signature, R1_sense_UMI, R2_sense_UMI, R1_seq, R2_seq, triage_result)

    return (processed_read_pair, R1_multiple_starts, R2_multiple_starts, has_swapped_sense_primers)


# For each read, collect
#  - primer beginning at ~7
#  - most downstream negative strand-matching primer (non-UDTD)
#
# Then, handle the cases
# Get insert and UMI(s)
def processLinesForCurrRead2(curr_read, curr_read_lines, pp_full_names):
    processed_read_pair = None
    assert (len(curr_read_lines)>0)
    orig_R1_seq, orig_R2_seq = None, None

    # Parse R1 
    R1_UDTD_start = None
    R1_pos_primer = None
    R1_neg_primer = None
    R1_neg_primer_candidates = []
    R1_seq = None
    R1_matches = filter(lambda x:x[-1]=="R1", curr_read_lines)
    R1_matches.sort(key=itemgetter(1), reverse=True)
    for R1_line_elems in R1_matches:
        if (len(R1_line_elems)==3):
            R1_seq, R1_err = R1_line_elems[0:2]
            orig_R1_seq = R1_seq[:]
            
        elif (R1_line_elems[0] == '+' and not R1_line_elems[4].startswith('UDTD')):
            if (R1_line_elems[1] == 7):
                if (R1_pos_primer != None): # If this fails, need to select primer with highest % identity
                    R1_pos_primer = R1_line_elems if (R1_line_elems[10] > R1_pos_primer[10]) else R1_pos_primer
                else:
                    R1_pos_primer = R1_line_elems

            elif (R1_line_elems[1] == 6 or R1_line_elems[1] == 8): 
                if (R1_pos_primer != None): # Select primer with highest % identity
                    R1_pos_primer = R1_line_elems if (R1_line_elems[10] > R1_pos_primer[10]) else R1_pos_primer
                else:
                    R1_pos_primer = R1_line_elems

        elif (R1_line_elems[0] == '-'):
            if (R1_line_elems[4].startswith('UDTD')):
                R1_UDTD_start = R1_line_elems[1] # TODO: subtract unmatched length
            else:
                R1_neg_primer_candidates.append( R1_line_elems )

    if (R1_UDTD_start != None and len(R1_neg_primer_candidates)>0):
        R1_neg_primer_candidates = list(filter(lambda R1_line_elems: R1_line_elems[2] < R1_UDTD_start, R1_neg_primer_candidates))

    if (R1_UDTD_start != None):
        R1_seq = R1_seq[0:R1_UDTD_start]


    # Parse R2
    R2_UDTD_start = None
    R2_pos_primer = None
    R2_neg_primer = None
    R2_neg_primer_candidates = []
    R2_seq = None
    R2_matches = filter(lambda x:x[-1]=="R2", curr_read_lines)
    R2_matches.sort(key=itemgetter(1), reverse=True)
    for R2_line_elems in R2_matches:
        if (len(R2_line_elems)==3):
            R2_seq, R2_err = R2_line_elems[0:2]
            orig_R2_seq = R2_seq[:]

        elif (R2_line_elems[0] == '+' and not R2_line_elems[4].startswith('UDTD')):
            if (R2_line_elems[1] == 7):
                if (R2_pos_primer != None): # Select primer with highest % identity
                    R2_pos_primer = R2_line_elems if (R2_line_elems[10] > R2_pos_primer[10]) else R2_pos_primer
                else:
                    R2_pos_primer = R2_line_elems

            elif (R2_line_elems[1] == 6 or R2_line_elems[1] == 8):
                if (R2_pos_primer != None): # Select primer with highest % identity
                    R2_pos_primer = R2_line_elems if (R2_line_elems[10] > R2_pos_primer[10]) else R2_pos_primer
                else:
                    R2_pos_primer = R2_line_elems

        elif (R2_line_elems[0] == '-'):
            if (R2_line_elems[4].startswith('UDTD')):
                R2_UDTD_start = R2_line_elems[1]
            else:
                R2_neg_primer_candidates.append( R2_line_elems )

    if (R2_UDTD_start != None and len(R2_neg_primer_candidates)>0):
        R2_neg_primer_candidates = list(filter(lambda R2_line_elems: R2_line_elems[2] < R2_UDTD_start, R2_neg_primer_candidates))

    if (R2_UDTD_start != None):
        R2_seq = R2_seq[0:R2_UDTD_start]

    neg_R2_in_R1_read = R2_pos_primer != None and any(map(lambda R1_neg_cand: R2_pos_primer[4] == R1_neg_cand[4], R1_neg_primer_candidates))
    neg_R1_in_R2_read = R1_pos_primer != None and any(map(lambda R2_neg_cand: R1_pos_primer[4] == R2_neg_cand[4], R2_neg_primer_candidates))
    
    if (neg_R2_in_R1_read and neg_R1_in_R2_read):
        R1_neg_primers = list(filter(lambda R1_neg_cand: R2_pos_primer[4] == R1_neg_cand[4], R1_neg_primer_candidates))
        R1_neg_primers.sort(key=itemgetter(2))
        R1_neg_primer = R1_neg_primers[-1]

        R2_neg_primers = list(filter(lambda R2_neg_cand: R1_pos_primer[4] == R2_neg_cand[4], R2_neg_primer_candidates))
        R2_neg_primers.sort(key=itemgetter(2))
        R2_neg_primer = R2_neg_primers[-1]


    if (R1_pos_primer != None):
        R1_sense_primer_src, R1_sense_equiv_region = R1_pos_primer[4].rsplit('_',1)
    else:
        R1_sense_primer_src, R1_sense_equiv_region = None, None
        
    if (R1_neg_primer != None):
        R1_antisense_primer_src, R1_antisense_equiv_region = R1_neg_primer[4].rsplit('_',1)
    else:
        R1_antisense_primer_src, R1_antisense_equiv_region = None, None


    if (R2_pos_primer != None):
        try:
            R2_sense_primer_src, R2_sense_equiv_region = R2_pos_primer[4].rsplit('_',1)
        except ValueError:
            pdb.set_trace()
            print(curr_read)
    else:
        R2_sense_primer_src, R2_sense_equiv_region = None, None
        
    if (R2_neg_primer != None):
        R2_antisense_primer_src, R2_antisense_equiv_region = R2_neg_primer[4].rsplit('_',1)
    else:
        R2_antisense_primer_src, R2_antisense_equiv_region = None, None


    # Categorize the primer matches
    R1_sense_primer_name = 'None' if (R1_pos_primer == None) else R1_pos_primer[4]
    R1_antisense_primer_name = 'None' if (R1_neg_primer == None) else R1_neg_primer[4]
    R2_sense_primer_name = 'None' if (R2_pos_primer == None) else R2_pos_primer[4]
    R2_antisense_primer_name = 'None' if (R2_neg_primer == None) else R2_neg_primer[4]

    primer_signature = "%s-%s" % (R1_sense_primer_name, R2_sense_primer_name)

    if (R1_sense_primer_src == R2_sense_primer_src != None):
        pp_full_name = "%s_%s%s" % (R1_sense_primer_src, R1_sense_equiv_region, R2_sense_equiv_region)
        triage_result = "ExpectedTarget" if (pp_full_name in pp_full_names) else "UnexpectedTarget"
        R1_sense_UMI = R1_seq[0:R1_pos_primer[1]-1]
        R2_sense_UMI = R2_seq[0:R2_pos_primer[1]-1]

        if (R1_antisense_primer_name == R2_sense_primer_name):
            R1_seq = R1_seq[R1_pos_primer[1]-1:R1_neg_primer[2]]
        else:
            R1_seq = R1_seq[R1_pos_primer[1]-1:]

        if (R2_antisense_primer_name == R1_sense_primer_name):
            R2_seq = R2_seq[R2_pos_primer[1]-1:R2_neg_primer[2]]
        else:
            R2_seq = R2_seq[R2_pos_primer[1]-1:]

    elif (R1_sense_primer_src != None and R2_sense_primer_src != None):
        triage_result = "Misprime"
        R1_sense_UMI = R1_seq[0:R1_pos_primer[1]-1]
        R2_sense_UMI = R2_seq[0:R2_pos_primer[1]-1]

        if (R1_antisense_primer_name == R2_sense_primer_name):
            R1_seq = R1_seq[R1_pos_primer[1]-1:R1_neg_primer[2]]
        else:
            R1_seq = R1_seq[R1_pos_primer[1]-1:]

        if (R2_antisense_primer_name == R1_sense_primer_name):
            R2_seq = R2_seq[R2_pos_primer[1]-1:R2_neg_primer[2]]
        else:
            R2_seq = R2_seq[R2_pos_primer[1]-1:]

    elif (R1_sense_primer_src == None and R2_sense_primer_src != None):
        triage_result = "Misprime"
        R1_sense_UMI = "None"
        R2_sense_UMI = R2_seq[0:R2_pos_primer[1]-1]
        R2_seq = R2_seq[R2_pos_primer[1]-1:]

    elif (R1_sense_primer_src != None and R2_sense_primer_src == None):
        triage_result = "Misprime"
        R1_sense_UMI = R1_seq[0:R1_pos_primer[1]-1]
        R2_sense_UMI = "None"
        R1_seq = R1_seq[R1_pos_primer[1]-1:]

    else:
        triage_result = "Misprime"
        R1_sense_UMI = "None"
        R2_sense_UMI = "None"

    return (primer_signature, R1_sense_UMI, R2_sense_UMI, R1_err, R2_err, R1_seq, R2_seq, triage_result)


def makePrimerSignatureForReadPair(R1_sig, R2_sig, pp_full_names):
    primer_signature = None
    triage_result = ""
    
    if (all(R1_sig) and all(R2_sig)):
        # Both primers in both reads. Determine if the primers should be paired and "on target" (expected from the subpool or not) or if they are a mispriming event.
        if (R1_sig == R2_sig):
            primer_signature = "%s_%s-%s_%s" % (R1_sig[0], R1_sig[1], R2_sig[2], R2_sig[3])
            if (R1_sig[0] == R1_sig[2]):
                pp_full_name = "%s_%s%s" % (R1_sig[0], R1_sig[1], R1_sig[3])
                triage_result = "ExpectedTarget" if (pp_full_name in pp_full_names) else "Misprime"
            else:
                triage_result = "Misprime"
        else:
            pdb.set_trace()
            triage_result = ''


    elif ((R1_sig[0:2] != (None, None) and R2_sig[2:4] != (None, None)) and (R1_sig[2:4] == R2_sig[0:2] == (None, None))):
        # Each read begins with a primer. No other primer matches found.
        # Determine if the primers should be paired and "on target" (expected from the subpool or not) or if they are a mispriming event.
        primer_signature = "%s_%s-%s_%s" % (R1_sig[0], R1_sig[1], R2_sig[2], R2_sig[3])
        if (R1_sig[0] == R2_sig[2]):
            pp_full_name = "%s_%s%s" % (R1_sig[0], R1_sig[1], R2_sig[3])
            triage_result = "ExpectedTarget" if (pp_full_name in pp_full_names) else "Misprime"
        else:
            triage_result = "Misprime"


    elif (all(R1_sig) and not any(R2_sig)):
        # Saw fwd & rev primers in R1, neither in R2. Assume R2 is bad.
        primer_signature = "%s_%s-%s_%s" % R1_sig
        triage_result = "Misprime"
    

    elif (not any(R1_sig) and all(R2_sig)):
        # Saw fwd & rev primers in R2, neither in R1. Assume R1 is bad.
        primer_signature = "%s_%s-%s_%s" % R2_sig
        triage_result = "Misprime"


    elif (all(R1_sig[0:2]) and not any(R1_sig[2:4]) and not any(R2_sig)):
        # Only saw fwd primer in R1 and no other primers.
        primer_signature = "%s_%s-None" % R1_sig[0:2]
        triage_result = "Misprime"
        

    elif (all(R2_sig[0:2]) and not any(R2_sig[2:4]) and not any(R1_sig)):
        # Only saw fwd primer in R2 and no other primers.
        primer_signature = "%s_%s-None" % R2_sig[0:2]
        triage_result = "Misprime"


    elif (not any(R1_sig[0:2]) and all(R1_sig[2:4]) and not any(R2_sig)):
        # rev primer only in R1 and no other primers.
        primer_signature = "None-%s_%s" % R1_sig[2:4]
        triage_result = "Misprime"
        

    elif (not any(R2_sig[0:2]) and all(R2_sig[2:4]) and not any(R1_sig)):
        # rev primer only in R2 and no other primers.
        primer_signature = "None-%s_%s" % R2_sig[2:4]
        triage_result = "Misprime"


    elif ((R1_sig[0:2] == R2_sig[0:2] == (None, None)) and R1_sig[2:4] != (None, None) and R2_sig[2:4] != (None, None)):
        #  rev primer in R1 & R2, no fwd primers
        if (R1_sig[2:4] == R2_sig[2:4]):
            primer_signature = "None-%s_%s" % R1_sig[2:4]
        else:
            primer_signature = "%s_%s-%s_%s" % (R2_sig[2], R2_sig[3], R1_sig[0], R1_sig[1])
        triage_result = "Misprime"


    elif ((R1_sig[0:2] == (None, None) or R2_sig[0:2] == (None, None)) and (R1_sig[2:4] == R2_sig[2:4] != (None, None))):
        #  rev primer in R1 & R2, fwd primer only in one
        if (R1_sig[0:2] == (None, None)):
            primer_signature = "%s_%s-%s_%s" % R2_sig
            if (R2_sig[0] != R1_sig[2]):
                triage_result = "Misprime"
            else:
                pp_full_name = "%s_%s%s" % (R2_sig[0], R2_sig[1], R2_sig[3])
                triage_result = "ExpectedTarget" if (pp_full_name in pp_full_names) else "Misprime"
        else:
            primer_signature = "%s_%s-%s_%s" % R1_sig
            if (R1_sig[0] != R2_sig[2]):
                triage_result = "Misprime"
            else:
                pp_full_name = "%s_%s%s" % (R1_sig[0], R1_sig[1], R1_sig[3])
                triage_result = "ExpectedTarget" if (pp_full_name in pp_full_names) else "Misprime"


    elif (R1_sig[0:2] != (None, None) and R2_sig[0:2] != (None, None) and (R1_sig[2:4] == R2_sig[2:4] == (None, None))):
        #  fwd primer in R1 & R2, no rev primers
        if (R1_sig[0:2] == R2_sig[0:2]):
            primer_signature = "%s_%s-None" % R1_sig[0:2]
        else:
            primer_signature = "%s_%s-%s_%s" % (R1_sig[0], R1_sig[1], R2_sig[2], R2_sig[3])
        triage_result = "Misprime"


    elif ((R1_sig[2:4] == (None, None) or R2_sig[2:4] == (None, None)) and (R1_sig[0:2] == R2_sig[0:2] != (None, None))):
        #  fwd primer in R1 & R2, rev primer only in one
        if (R1_sig[2:4] == (None, None)):
            primer_signature = "%s_%s-%s_%s" % R2_sig
            if (R1_sig[0] != R2_sig[2]):
                triage_result = "Misprime"
            else:
                pp_full_name = "%s_%s%s" % (R2_sig[0], R2_sig[1], R2_sig[3])
                triage_result = "ExpectedTarget" if (pp_full_name in pp_full_names) else "Misprime"
        else:
            primer_signature = "%s_%s-%s_%s" % R1_sig
            if (R1_sig[0] != R1_sig[2]):
                triage_result = "Misprime"
            else:
                pp_full_name = "%s_%s%s" % (R1_sig[0], R1_sig[1], R1_sig[3])
                triage_result = "ExpectedTarget" if (pp_full_name in pp_full_names) else "Misprime"


    elif (R1_sig == R2_sig == (None, None, None, None)):
        # No primers matches in either read
        primer_signature = "None-None"
        triage_result = "Misprime"
    

    else:
        primer_signature = "Unaccounted-Unaccounted"
        triage_result = "Misprime"


    #assert (primer_signature != None and triage_result != "")
    
    return (primer_signature, triage_result)


def groupReads(merged_sorted_data, pp_full_names):
    all_nonunique_umi_instances_per_signature = defaultdict(list)
    grouped_reads = defaultdict(list)
    triage_results = {}

    ip = open(merged_sorted_data, 'r')

    curr_read = None
    curr_read_lines = []
    line = ip.readline()

    #all_R1_multiple_starts = set()
    #all_R2_multiple_starts = set()
    #swapped_sense_primers = set()

    None_count, nonNone_count = 0, 0
    #counter = 0
    #can_continue = False
    while (line != ''):
        fields = line.strip().split("\t")
        read_ID = fields[0].split()[0]

        #if (can_continue or read_ID=="M01048:156:000000000-BKB6V:1:1105:21593:14125"):
        #    can_continue = True
        #else:
        #    line = ip.readline()
        #    continue

        if (curr_read == None):
            curr_read = read_ID

        if (read_ID == curr_read):
            if (len(fields)==4):
                curr_read_lines.append( fields[1:] )
            elif (len(fields) > 1):
                fields = map(lambda x: x[0](x[1]), zip((str,str,int,int,int,str,str,int,int,int,int,float,str),fields))
                curr_read_lines.append( fields[1:] )
        else:
            primer_signature, R1_sense_UMI, R2_sense_UMI, R1_err, R2_err, R1_insert, R2_insert, triage_result = processLinesForCurrRead2(curr_read, curr_read_lines, pp_full_names)
            #processed_read_pair, R1_multiple_starts, R2_multiple_starts, has_swapped_sense_primers
            #all_R1_multiple_starts.update(R1_multiple_starts)
            #all_R2_multiple_starts.update(R2_multiple_starts)
            #if (has_swapped_sense_primers):
            #    swapped_sense_primers.add(curr_read)

            #if (processed_read_pair != None):
                #primer_signature, R1_umi, R2_umi, R1_insert, R2_insert, triage_result = processed_read_pair
                #if (primer_signature == 'chr7_OS5_0_F39-chr3_OS114_0_R38' and triage_result == 'ExpectedTarget'):
                #    pdb.set_trace()
                    
            if (primer_signature in triage_results):
                try:
                    assert (triage_results[primer_signature] == triage_result)
                except AssertionError:
                    pdb.set_trace()
            else:
                triage_results[primer_signature] = triage_result

            grouped_reads[primer_signature].append( (primer_signature, R1_sense_UMI, R2_sense_UMI, R1_err, R2_err, R1_insert, R2_insert) )
            all_nonunique_umi_instances_per_signature[primer_signature].append( (R1_sense_UMI, R2_sense_UMI) )
            #nonNone_count += 1
            #else:
            #    pdb.set_trace()
            #    None_count += 1

            curr_read = read_ID
            if (len(fields[1:])!=0): # Assert should always be len 0
                curr_read_lines = [ fields[1:] ]
            else:
                curr_read_lines = []

            #counter += 1
            #if (counter == 10000):
            #    break

        line = ip.readline()

    ip.close()

    primer_signature, R1_sense_UMI, R2_sense_UMI, R1_err, R2_err, R1_insert, R2_insert, triage_result = processLinesForCurrRead2(curr_read, curr_read_lines, pp_full_names)
    #processed_read_pair, R1_multiple_starts, R2_multiple_starts, has_swapped_sense_primers = processLinesForCurrRead(curr_read, curr_read_lines, pp_full_names)
    #all_R1_multiple_starts.update(R1_multiple_starts)
    #all_R2_multiple_starts.update(R2_multiple_starts)
    #if (has_swapped_sense_primers):
    #    swapped_sense_primers.add(curr_read)

    #if (processed_read_pair != None):
    #    primer_signature, R1_umi, R2_umi, R1_insert, R2_insert = processed_read_pair
    if (primer_signature in triage_results):
        assert (triage_results[primer_signature] == triage_result)
    else:
        triage_results[primer_signature] = triage_result

    grouped_reads[primer_signature].append( (primer_signature, R1_sense_UMI, R2_sense_UMI, R1_err, R2_err, R1_insert, R2_insert) )
    all_nonunique_umi_instances_per_signature[primer_signature].append( (R1_sense_UMI, R2_sense_UMI) )
    #    nonNone_count += 1
    #else:
    #    None_count += 1

    #print("Number of R1 with multiple starts = %s" % len(all_R1_multiple_starts))
    #print("Number of R2 with multiple starts = %s" % len(all_R2_multiple_starts))
    #print("Number of read pairs with swapped_sense_primers = %d" % len(swapped_sense_primers))

    return (grouped_reads, all_nonunique_umi_instances_per_signature, triage_results) #, None_count, nonNone_count)


def writeSummaryAndGroupedReads(summary_file, category_summary_file, dir_for_grouping, grouped_reads,
                                all_nonunique_umi_instances_per_signature, triage_results, pp_full_names, subpool_per_pp, target_subpools):
    for primer_signature, tuples_list in grouped_reads.items():
        filename = "%s/%s.tsv" % (dir_for_grouping,primer_signature)
        with open(filename, 'w') as op:
            for tup in tuples_list:
                op.write("%s\t%s\t%s\t%s\t%s\t%s\n" % tup[1:])

    signature_counts = []
    signature_category_counts = defaultdict(lambda: (0,0))
    for primer_signature, umis in all_nonunique_umi_instances_per_signature.items():
        triage_result = triage_results[primer_signature]
        unique_umis = set(umis)
        signature_counts.append( (primer_signature, triage_result, len(umis), len(unique_umis)) )

        if (primer_signature == "None-None"):
            sig_cat = "None-None\tMisprime"
        elif (primer_signature.startswith("None")):
            sig_cat = "None-Rev\tMisPrime"
        elif (primer_signature.endswith("None")):
            sig_cat = "Fwd-None\tMisPrime"
        elif (triage_result == "Misprime"):
            sig_cat = "Fwd-Rev\tMisPrime"
        elif (triage_result == "ExpectedTarget"):
            sig_cat = "Fwd-Rev\tExpectedTarget"
        elif (triage_result == "UnexpectedTarget"):
            sig_cat = "Fwd-Rev\tUnexpectedTarget"
        
        Numi, Nunique_umi = signature_category_counts[sig_cat]
        signature_category_counts[sig_cat] = (Numi+len(umis), Nunique_umi+len(unique_umis))

    Numi_total = sum(map(itemgetter(0), signature_category_counts.values()))
    Nunique_umi_total = sum(map(itemgetter(1), signature_category_counts.values()))

    op = open(category_summary_file, 'w')
    for sig_cat in ["Fwd-Rev\tExpectedTarget", "Fwd-Rev\tMisPrime", "Fwd-None\tMisPrime",
                    "None-Rev\tMisPrime", "None-None\tMisprime", "Fwd-Rev\tUnexpectedTarget"]:
        Numi, Nunique_umi = signature_category_counts[sig_cat]
        perc_Numi = 100.0 * float(Numi)/float(Numi_total)
        perc_Nunique_umi = 100.0 * float(Nunique_umi)/float(Nunique_umi_total)  
        op.write("%s\t%3.1f%%\t%3.1f%%\n" % (sig_cat, perc_Numi, perc_Nunique_umi))
    op.close()

    signature_counts.sort(key=itemgetter(2,3), reverse=True)
    op = open(summary_file, 'w')
    for tup in signature_counts:
        op.write("%s\t%s\t%d\t%d\n" % tup)
    op.close()
    

def writeGroupedReadsAndPrimingEventsReport(grouped_reads, triage_results, all_nonunique_umi_instances_per_signature, dir_for_grouping, priming_events_output):
    tups_to_write = []
    for primer_signature in grouped_reads.keys():
        triage_result = triage_results[primer_signature]
        count = len(all_nonunique_umi_instances_per_signature[primer_signature])
        
        if (primer_signature == "None-None"):
            sig_cat = "None-None\tMisprime"
        elif (primer_signature.startswith("None")):
            sig_cat = "None-Rev\tMisPrime"
        elif (primer_signature.endswith("None")):
            sig_cat = "Fwd-None\tMisPrime"
        elif (triage_result == "Misprime"):
            sig_cat = "Fwd-Rev\tMisPrime"
        elif (triage_result == "ExpectedTarget"):
            sig_cat = "Fwd-Rev\tExpectedTarget"
        elif (triage_result == "UnexpectedTarget"):
            sig_cat = "Fwd-Rev\tUnexpectedTarget"

        tups_to_write.append( (primer_signature, triage_result, count) )


    tups_to_write.sort(key=itemgetter(2), reverse=True)
    with open(priming_events_output, 'w') as op:
        for tup in tups_to_write:
            op.write("%s\t%s\t%d\n" % tup)


    for primer_signature, tuples_list in grouped_reads.items():
        filename = "%s/%s.tsv" % (dir_for_grouping,primer_signature)
        with open(filename, 'w') as op:
            for tup in tuples_list:
                op.write("%s\t%s\t%s\t%s\t%s\t%s\n" % tup[1:])


if (__name__ == "__main__"):
    target_subpools, merged_sorted_data, pp_names_file, pp_subpools_file, dir_for_grouping, priming_events_output = sys.argv[1:]
    target_subpools = set(target_subpools.split(','))

    if (pp_subpools_file != 'None'):
        subpool_per_pp = readSubpools(pp_subpools_file)
    else:
        subpool_per_pp = None
    pp_full_names = readNamesExpectedPrimerPairs(pp_names_file)
    grouped_reads, all_nonunique_umi_instances_per_signature, triage_results = groupReads(merged_sorted_data, pp_full_names)
    writeGroupedReadsAndPrimingEventsReport(grouped_reads, triage_results, all_nonunique_umi_instances_per_signature, dir_for_grouping, priming_events_output)

    sys.exit(0)
