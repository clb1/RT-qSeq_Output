#!/usr/bin/env python

from collections import defaultdict
from operator import itemgetter
import sys


def compileUMIsPerTarget():
    UMIs_per_target = defaultdict(list)

    ip = open(sys.stdin, 'r')
    for line in ip:
        fields = line.strip().split("\t")
        assert (fields[0] == fields[3]), "Read names do not match"
        assert (fields[1] == fields[4]), "Target names do not match"
        UMIs_per_target[fields[1]].append( (fields[2], fields[5]) )
    ip.close()

    return UMIs_per_target


def writeCounts(UMIs_per_target, output_file):
    L = []
    for target_label, counts in UMIs_per_target.items():
        nonunique_counts = len(counts)
        unique_counts = len(set(counts))
        L.append( (target_label, nonunique_counts, unique_counts) )

    L.sort(key=itemgetter(1,2), reverse=True)

    op = open(output_file, 'w')
    for tup in L:
        op.write("%s\t%d\t%d\n" % tup)
    op.close()

    
if (__name__ == "__main__"):
    output_file = sys.argv[1]

    UMIs_per_target = compileUMIsPerTarget()
    writeCounts(UMIs_per_target, output_file)

    sys.exit(0)
