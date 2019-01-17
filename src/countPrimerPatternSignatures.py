#!/usr/bin/env python

from operator import itemgetter
from collections import Counter, defaultdict
import sys


def standardizePrimerName(full_primer_name):
    if ("CTRL" in full_primer_name or full_primer_name.startswith("UDTD") or full_primer_name.startswith("ERCC") or
        full_primer_name.startswith("P5") or full_primer_name.startswith("P7")):
        standardized_name = full_primer_name
    else:
        elems = full_primer_name.split('_')

        try:
            F, R = elems[3].split('R')
        except IndexError as ie:
            pdb.set_trace()
            
        R = "R" + R
        if (elems[-1] == "5p"):
            standardized_name = "%s_%s_%s_%s" % (elems[0], elems[1], elems[2], F)
        else:
            standardized_name = "%s_%s_%s_%s" % (elems[0], elems[1], elems[2], R)

    return standardized_name



def processLinesForCurrRead(curr_read_ID, curr_read_lines, pattern_counts):
    ppr_5p_occurrences = set()
    ppr_3p_occurrences = set()

    if (len(curr_read_lines)==0):
        pattern_counts.update(["no primer matches"])
    else:
        R1_lines = filter(lambda x:x[-1]=="R1", curr_read_lines)
        R1_matches = []
        for R1_line_elems in R1_lines:
            primer_name = R1_line_elems[0]
            standardized_primer_name = standardizePrimerName(primer_name)
            if (primer_name.endswith("_5p")):
                ppr_5p_occurrences.add( primer_name[0:-3] )
            elif (primer_name.endswith("_3p")):
                ppr_3p_occurrences.add( primer_name[0:-3] )

            
            new_match = (standardized_primer_name, R1_line_elems[1], int(R1_line_elems[3]), int(R1_line_elems[4]))
            if (new_match not in R1_matches):
                R1_matches.append( new_match )
        R1_matches.sort(key=itemgetter(2))
        R1_matches = map(itemgetter(0,1), R1_matches)    # Removes positions
        R1_pattern = " ".join(map(str, R1_matches))

        R2_lines = filter(lambda x:x[-1]=="R2", curr_read_lines)
        R2_matches = []
        for R2_line_elems in R2_lines:
            primer_name = R2_line_elems[0]
            if (primer_name.endswith("_5p")):
                ppr_5p_occurrences.add( primer_name[0:-3] )
            elif (primer_name.endswith("_3p")):
                ppr_3p_occurrences.add( primer_name[0:-3] )

            standardized_primer_name = standardizePrimerName(primer_name)
            new_match = (standardized_primer_name, R2_line_elems[1], int(R2_line_elems[3]), int(R2_line_elems[4]))
            if (new_match not in R2_matches):
                R2_matches.append( new_match )
        R2_matches.sort(key=itemgetter(2))
        R2_matches = map(itemgetter(0,1), R2_matches)    # Removes positions
        R2_pattern = " ".join(map(str, R2_matches))

        common_5p_3p_ppr = ppr_5p_occurrences & ppr_3p_occurrences
        if (len(common_5p_3p_ppr) == 1):
            common_ppr = list(common_5p_3p_ppr)[0]
        elif (len(common_5p_3p_ppr) > 1):
            print >> sys.stderr, "More than one ppr primed, read %s" % curr_read_ID
            common_5p_3p_ppr = sorted(common_5p_3p_ppr)
            common_ppr = ",".join(common_5p_3p_ppr)
        else:
            common_ppr = "No common ppr"

        R1_R2_pattern = "%s | %s\t%s" % (R1_pattern, R2_pattern, common_ppr)
        pattern_counts.update([R1_R2_pattern])


def processLinesForCurrRead2(curr_read_ID, curr_read_lines, pattern_counts, umi_counts):
    fwd_matches = set()
    rev_matches = set()
    R1_seq = None
    R2_seq = None
    min_R1_pos = 1e10
    min_R2_pos = 1e10

    if (len(curr_read_lines)==0):
        pattern_counts.update(["no primer matches"])
    else:
        for R1_line_elems in filter(lambda x:x[-1]=="R1", curr_read_lines):
            if (len(R1_line_elems)==2):
                R1_seq = R1_line_elems[0]
            else:
                standardized_primer_name = standardizePrimerName(R1_line_elems[4])
                if ("CTRL" in standardized_primer_name or standardized_primer_name.startswith("chr")):
                    min_R1_pos = min(min_R1_pos, int(R1_line_elems[1]))
                assert (R1_line_elems[0] == '+' or R1_line_elems[0] == '-')
                if (R1_line_elems[0] == '+'):
                    fwd_matches.add(standardized_primer_name)
                else:
                    rev_matches.add(standardized_primer_name)

        for R2_line_elems in filter(lambda x:x[-1]=="R2", curr_read_lines):
            if (len(R2_line_elems)==2):
                R2_seq = R2_line_elems[0]
            else:
                standardized_primer_name = standardizePrimerName(R2_line_elems[4])
                if ("CTRL" in standardized_primer_name or standardized_primer_name.startswith("chr")):
                    min_R2_pos = min(min_R2_pos, int(R2_line_elems[1]))
                assert (R2_line_elems[0] == '+' or R2_line_elems[0] == '-')
                if (R2_line_elems[0] == '+'):
                    rev_matches.add(standardized_primer_name)
                else:
                    fwd_matches.add(standardized_primer_name)

        
        umi_1 = R1_seq[0:min_R1_pos] if (min_R1_pos < 1e10) else ''
        umi_2 = R2_seq[0:min_R2_pos] if (min_R2_pos < 1e10) else ''
        umi_12 = "%s,%s" % (umi_1,umi_2)

        fwd_pattern = ' '.join(sorted(fwd_matches))
        rev_pattern = ' '.join(sorted(rev_matches))

        fwd_rev_pattern = "%s || %s" % (fwd_pattern, rev_pattern)
        pattern_counts.update([fwd_rev_pattern])
        umi_counts[fwd_rev_pattern].append(umi_12)
        

def compilePatternCounts(merged_sorted_data):
    pattern_counts = Counter()
    umi_counts = defaultdict(list)

    ip = open(merged_sorted_data, 'r')

    curr_read_ID = None
    curr_read_lines = []

    for line in ip:
        fields = line.strip().split("\t")
        read_ID = fields[0].split()[0]

        if (curr_read_ID == None):
            curr_read_ID = read_ID

        if (read_ID == curr_read_ID):
            if (len(fields) != 1):
                curr_read_lines.append( fields[1:] )
        else:
            processLinesForCurrRead2(curr_read_ID, curr_read_lines, pattern_counts, umi_counts)
            curr_read_ID = read_ID
            if (len(fields) > 1):
                curr_read_lines = [ fields[1:] ]
            else:
                curr_read_lines = []

    ip.close()

    processLinesForCurrRead2(curr_read_ID, curr_read_lines, pattern_counts, umi_counts)

    return (pattern_counts, umi_counts)


def writePatternCounts(pattern_counts, umi_counts, output_file):
    op = open(output_file, 'w')
    for pattern, pattern_count in pattern_counts.most_common():
        if ("RBM17_CTRL_5p || RBM17_CTRL_3p" in pattern):
            pdb.set_trace()
        try:
            umi_counter = Counter(umi_counts[pattern])
            most_common_umi, umi_count = umi_counter.most_common(1)[0]
            perc = 100.0 * float(umi_count)/float(pattern_count)
            op.write("%s\t%d\t%5.3f\t%s\n" % (pattern, pattern_count, perc, most_common_umi))
        except TypeError:
            pdb.set_trace()
            op.write("%s\t%d\t%-\n" % (pattern, pattern_count))

    op.close()


if (__name__ == "__main__"):
    merged_sorted_data, output_file = sys.argv[1:]

    pattern_counts, umi_counts = compilePatternCounts(merged_sorted_data)
    writePatternCounts(pattern_counts, umi_counts, output_file)

    sys.exit(0)
