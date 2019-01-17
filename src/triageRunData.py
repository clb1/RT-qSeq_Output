#!/usr/bin/env python

from operator import itemgetter
from collections import Counter
import sys


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



def processLinesForCurrRead(curr_read, curr_read_lines, pattern_counts):
    if (len(curr_read_lines)==0):
        pattern_counts.update(["no primer matches"])
    else:
        R1_lines = filter(lambda x:x[-1]=="R1", curr_read_lines)
        R1_matches = []
        for R1_line_elems in R1_lines:
            standardized_name = standardizePrimerName(R1_line_elems[0])
            new_match = (standardized_name, R1_line_elems[1], int(R1_line_elems[3]), int(R1_line_elems[4]))
            if (new_match not in R1_matches):
                R1_matches.append( new_match )
        R1_matches.sort(key=itemgetter(2))
        R1_matches = map(itemgetter(0,1), R1_matches)    # Removes positions
        R1_pattern = " ".join(map(str, R1_matches))

        R2_lines = filter(lambda x:x[-1]=="R2", curr_read_lines)
        R2_matches = []
        for R2_line_elems in R2_lines:
            standardized_name = standardizePrimerName(R2_line_elems[0])
            new_match = (standardized_name, R2_line_elems[1], int(R2_line_elems[3]), int(R2_line_elems[4]))
            if (new_match not in R2_matches):
                R2_matches.append( new_match )
        R2_matches.sort(key=itemgetter(2))
        R2_matches = map(itemgetter(0,1), R2_matches)    # Removes positions
        R2_pattern = " ".join(map(str, R2_matches))

        if (R1_pattern == ''):
            pattern_counts.update([R2_pattern])
            #print >> sys.stderr, R2_pattern
        elif (R2_pattern == ''):
            pattern_counts.update([R1_pattern])
            #print >> sys.stderr, R1_pattern
        else:
            R1_R2_pattern = "%s | %s" % (R1_pattern, R2_pattern)
            R2_R1_pattern = "%s | %s" % (R2_pattern, R1_pattern)
            if (pattern_counts.has_key(R1_R2_pattern)):
                pattern_counts.update([R1_R2_pattern])
                #print >> sys.stderr, R1_R2_pattern
            elif (pattern_counts.has_key(R2_R1_pattern)):
                pattern_counts.update([R2_R1_pattern])
                #print >> sys.stderr, R2_R1_pattern
            else:
                pattern_counts.update([R1_R2_pattern])
                #print >> sys.stderr, R1_R2_pattern


def compilePatternCounts(merged_sorted_data):
    pattern_counts = Counter()

    ip = open(merged_sorted_data, 'r')

    curr_read = None
    curr_read_lines = []
    line = ip.readline()

    while (line != ''):
        fields = line.strip().split("\t")
        read_ID = fields[0].split()[0]

        if (curr_read == None):
            curr_read = read_ID

        if (read_ID == curr_read):
            if (len(fields[1:])!=0):
                curr_read_lines.append( fields[1:] )
        else:
            processLinesForCurrRead(curr_read, curr_read_lines, pattern_counts)
            curr_read = read_ID
            if (len(fields[1:])!=0): # Assert should always be len 0
                curr_read_lines = [ fields[1:] ]
            else:
                curr_read_lines = []

        line = ip.readline()

    ip.close()

    processLinesForCurrRead(curr_read, curr_read_lines, pattern_counts)

    return pattern_counts


def writePatternCounts(pattern_counts, output_file):
    op = open(output_file, 'w')
    for pattern, count in pattern_counts.most_common():
        op.write("%s\t%d\n" % (pattern, count))
    op.close()


if (__name__ == "__main__"):
    merged_sorted_data, output_file = sys.argv[1:]

    pattern_counts = compilePatternCounts(merged_sorted_data)
    writePatternCounts(pattern_counts, output_file)

    sys.exit(0)
