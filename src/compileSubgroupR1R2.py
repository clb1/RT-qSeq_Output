#!/usr/bin/env python3

import os
import re
import sys


def getPrimerSignaturesMetadata(primer_signatures_summary_file):
    primer_signatures_meta = {}
    with open(primer_signatures_summary_file, 'r') as ip:
        header = ip.readline()
        for line in ip:
            primer_signature, triage, num_subgroups, num_nonunique_reads, num_unique_reads = line.strip().split("\t")
            if (triage == "ExpectedAtLocus"):
                fwd_primer_signature, rev_primer_signature = primer_signature.split('-')
                has_correct_primer_w_adapter = "_F" in fwd_primer_signature and "_R" in rev_primer_signature
                if (primer_signature != "None-None" and has_correct_primer_w_adapter):
                    primer_signatures_meta[primer_signature] = (triage, int(num_subgroups), int(num_nonunique_reads), int(num_unique_reads))
                elif ("None" not in primer_signature):
                    print("WARNING: skipping, adapter-primer mismatch %s" % primer_signature, file=sys.stderr)
    return primer_signatures_meta


def collectR1R2(dir_for_grouping, primer_signatures_meta, read_length):
    all_R1, all_R2 = [], []
    reSubgroupHeader = re.compile("^(\S+)\s+(\S+)\s+(\d+) subgroups\s+(\d+) nonunique UMIs\s+(\d+) unique UMIs$")

    for primer_signature in primer_signatures_meta.keys():
        subgroups_file = "%s/%s_subgroups.tsv" % (dir_for_grouping, primer_signature)
        assert (os.path.exists(subgroups_file)), "Subgroups file %s does not exist" % subgroups_file
        with open(subgroups_file, 'r') as ip:
            header = ip.readline()
            header_primer_signature, triage, num_subgroups, total_num_nonunique_reads, total_num_unique_reads = reSubgroupHeader.match(header).groups()
            assert (header_primer_signature == primer_signature)
            for line in ip:
                subgroup_index, R1, R2, num_nonunique_reads, num_unique_reads = line.strip().split("\t")
                subgroup_index = int(subgroup_index)
                all_R1.append( ">%s-sg%d\n%s\n" % (primer_signature, subgroup_index, R1[0:read_length]) )
                all_R2.append( ">%s-sg%d\n%s\n" % (primer_signature, subgroup_index, R2[0:read_length]) )
            assert (subgroup_index == int(num_subgroups)), "Did not find expected number of subgroups for %s" % subgroup

    return (all_R1, all_R2)


def writeReads(R, R_fasta):
    with open(R_fasta, 'w') as op:
        for line in R:
            op.write(line)


if (__name__ == "__main__"):
    read_length, dir_for_grouping, primer_signatures_summary_file, R1_fasta, R2_fasta = sys.argv[1:]
    read_length = int(read_length)
    
    primer_signatures_meta = getPrimerSignaturesMetadata(primer_signatures_summary_file)
    all_R1, all_R2 = collectR1R2(dir_for_grouping, primer_signatures_meta, read_length)
    writeReads(all_R1, R1_fasta)
    writeReads(all_R2, R2_fasta)

    sys.exit(0)
    
