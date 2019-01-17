#!/usr/bin/env python3

from collections import Counter, defaultdict
import difflib
from itertools import product
from operator import itemgetter
import os
import subprocess
from subprocess import check_output, CalledProcessError, Popen
import sys
import tempfile

import pickle


def subgroupReadsFromEachPrimingEvent(scratch_dir, grouping_dir, priming_events_file):
    sig_and_triage = []
    
    with open(priming_events_file, 'r') as ip:
        for line in ip:
            primers_signature, triage_result, num_read_pairs = line.strip().split("\t")
            sig_and_triage.append( (primers_signature, triage_result, int(num_read_pairs)) )
        
    all_subgroup_summaries = []
    for count, (primers_signature, triage_result, num_read_pairs) in enumerate(sig_and_triage,1):
        print("%d of %d\t%s\tsize = %d" % (count, len(sig_and_triage), primers_signature, num_read_pairs), file=sys.stderr, flush=True)
        grouped_reads_file = "%s/%s.tsv" % (grouping_dir, primers_signature)
        output_file = "%s/%s_subgroups.tsv" % (grouping_dir, primers_signature)

        if (os.path.exists(output_file)):
            continue
        
        print("\tReading UMIs and read pairs", file=sys.stderr, flush=True)
        all_lines, umi_len = readUMIsAndReadPairs(grouped_reads_file)
        print("\tCorrecing UMIs", file=sys.stderr, flush=True)
        umi_corrected_lines, nonunique_umi_class_counts = correctUMIs(all_lines, umi_len)

        print("\tMaking consensus read pair for each UMI pair", file=sys.stderr, flush=True)
        read_pair_err, read_pair_seqs, read_pair_IDs_per_cluster, orig_read_pair_IDs_per_cluster, umi_map_cluster = makeConsensusReadPairsForEachUMI(umi_corrected_lines, scratch_dir)

        print("\tComputing read consensi", file=sys.stderr, flush=True)
        R1_R2_consensi_map_umis = computeReadConsensiIndividually(umi_map_cluster, orig_read_pair_IDs_per_cluster, read_pair_seqs, nonunique_umi_class_counts, scratch_dir)

        print("\tWriting subgroups", file=sys.stderr, flush=True)
        subgroup_summary = writeSubgroupResults(R1_R2_consensi_map_umis, primers_signature, triage_result, output_file) # expected_amplicons.bed/tsv
        all_subgroup_summaries.append(subgroup_summary)

    return all_subgroup_summaries


def readUMIsAndReadPairs(grouped_reads_file):
    all_lines = []
    umi_lens = []
    with open(grouped_reads_file,'r') as ip:
        for line in ip:
            fields = line.strip().split("\t")
            all_lines.append(fields)
            umi_lens.append(len(fields[0]))
            umi_lens.append(len(fields[1]))
    C = Counter(umi_lens)
    most_freq_umi_len = C.most_common()[0][0]
    return (all_lines, most_freq_umi_len)


def correctUMIs(all_lines, umi_len):
    umi_map_class = {}

    unclassified_umis = list(filter(lambda x:len(x[0])==len(x[1])==umi_len, map(itemgetter(0,1), all_lines)))
    unclassified_nonstandard_umis = list(filter(lambda x:len(x[0])!=umi_len or len(x[1])!=umi_len, map(itemgetter(0,1), all_lines)))

    # Classify UMIs that are of the expected UMI size
    while (len(unclassified_umis)>0):
        umi_counts = Counter(unclassified_umis)

        # The most frequent UMI pair is the next UMI class
        (R1_UMI, R2_UMI), umi_count = umi_counts.most_common(1)[0]

        # Add all exact and similar unclassified UMI pairs to the new UMI class
        to_remove = set()
        for (R1_umi, R2_umi) in unclassified_umis:
            if ((R1_umi,R2_umi) not in umi_map_class):  # TODO: is this 'if' necessary? No (R1_umi,R2_umi) should be in umi_map_class. Should have been removed.
                num_R1_diffs = sum(map(lambda i: R1_umi[-i]!=R1_UMI[-i], range(umi_len)))
                num_R2_diffs = sum(map(lambda i: R2_umi[-i]!=R2_UMI[-i], range(umi_len)))
                num_diffs = num_R1_diffs + num_R2_diffs
                if (num_diffs <= 2):
                    umi_map_class[(R1_umi,R2_umi)] = (R1_UMI, R2_UMI)
                    to_remove.add( (R1_umi,R2_umi) )

        # Filter out all of the UMI pairs that have now been classified
        unclassified_umis = list(filter(lambda a: a not in to_remove, unclassified_umis))


    # Classify UMIs that are of an unexpected UMI size
    sm = difflib.SequenceMatcher()
    umi_classes_R1 = set(map(itemgetter(0), umi_map_class.values()))
    umi_classes_R2 = set(map(itemgetter(1), umi_map_class.values()))

    while (len(unclassified_nonstandard_umis)>0):
        #print(len(unclassified_nonstandard_umis))
        umi_counts = Counter(unclassified_nonstandard_umis)

        # Get the most frequent nonstandard UMI pair, which may or may not form a new UMI class
        (R1_umi, R2_umi), umi_count = umi_counts.most_common(1)[0]

        class_matches_R1 = list(enumerate(difflib.get_close_matches(R1_umi, umi_classes_R1, 10, 0.66)))
        class_matches_R2 = list(enumerate(difflib.get_close_matches(R2_umi, umi_classes_R2, 10, 0.66)))

        was_put_in_existing_class = False
        if (len(class_matches_R1)>0 and len(class_matches_R2)>0):
            ranked_R1R2_matches = list(map(lambda x: (x[0][1],x[1][1],x[0][0]+x[1][0]), product(class_matches_R1, class_matches_R2)))
            ranked_R1R2_matches.sort(key=itemgetter(2))
            for (R1_UMI, R2_UMI, rank_sum) in ranked_R1R2_matches:
                if ((R1_UMI, R2_UMI) in umi_map_class.values()):
                    if (len(R1_umi) == len(R1_UMI)):
                        assert( len(R1_umi)==umi_len )
                        num_R1_diffs = sum(map(lambda i: R1_umi[-i]!=R1_UMI[-i], range(umi_len)))
                    else:
                        sm.set_seqs(R1_umi,R1_UMI)
                        num_R1_matches = int(round(sm.ratio() * 0.5 * float(len(R1_umi)+len(R1_UMI))))
                        num_R1_diffs = len(R1_umi)-num_R1_matches if (len(R1_umi)<len(R1_UMI)) else len(R1_UMI)-num_R1_matches

                    if (len(R2_umi) == len(R2_UMI)):
                        assert( len(R2_umi)==umi_len )
                        num_R2_diffs = sum(map(lambda i: R2_umi[-i]!=R2_UMI[-i], range(umi_len)))
                    else:
                        sm.set_seqs(R2_umi,R2_UMI)
                        num_R2_matches = int(round(sm.ratio() * 0.5 * float(len(R2_umi)+len(R2_UMI))))
                        num_R2_diffs = len(R2_umi)-num_R2_matches if (len(R2_umi)<len(R2_UMI)) else len(R2_UMI)-num_R2_matches

                    if (num_R1_diffs + num_R2_diffs <= 2):
                        #print("Added to class %s %s" % (R1_UMI, R2_UMI))
                        # The UMI pair is similar enough to an existing UMI class that it likely belongs to it
                        umi_map_class[(R1_umi,R2_umi)] = (R1_UMI, R2_UMI)
                        was_put_in_existing_class = True
                        break

        if (not was_put_in_existing_class):
            umi_map_class[(R1_umi,R2_umi)] = (R1_umi, R2_umi)

        unclassified_nonstandard_umis = list(filter(lambda a: a != (R1_umi,R2_umi), unclassified_nonstandard_umis))

    umi_corrected_lines = []
    for (R1_umi, R2_umi, R1_err, R2_err, R1_seq, R2_seq) in all_lines:
        corrected_R1_umi, corrected_R2_umi = umi_map_class[(R1_umi,R2_umi)]
        umi_corrected_lines.append( (corrected_R1_umi, corrected_R2_umi, float(R1_err), float(R2_err), R1_seq, R2_seq) )

    nonunique_umi_class_counts = Counter(map(itemgetter(0,1), umi_corrected_lines))

    return (umi_corrected_lines, nonunique_umi_class_counts)


def makeConsensusReadPairsForEachUMI(umi_corrected_lines, scratch_dir):
    read_pair_err = {}

    seqs_per_umi = defaultdict(list)
    for (R1_umi, R2_umi, R1_err, R2_err, R1_seq, R2_seq) in umi_corrected_lines:
        seqs_per_umi[(R1_umi, R2_umi)].append( (R1_err+R2_err, R1_seq, R2_seq) )
    
    # Write out all read pairs in order of lower expected number of errors
    read_pairs_to_write = []
    read_pair_seqs = {}
    for (R1_umi,R2_umi), L in seqs_per_umi.items():
        L.sort(key=itemgetter(0))
        for counter, (sum_err, R1_seq, R2_seq) in enumerate(L,1):
            read_pair_ID = "%s_%s_%d" % (R1_umi,R2_umi,counter)
            read_pair_err[read_pair_ID] = sum_err
            read_pairs_to_write.append( (sum_err, ">%s\n%s%s\n" % (read_pair_ID,R1_seq,R2_seq)) )
            read_pair_seqs[read_pair_ID] = (R1_seq,R2_seq)

    read_pairs_to_write.sort(key=itemgetter(0))
    read_pairs_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=scratch_dir, delete=True)
    for sum_err, fasta_entry in read_pairs_to_write:
        read_pairs_fasta.write(fasta_entry)
    read_pairs_fasta.flush()

    with tempfile.TemporaryDirectory(dir=scratch_dir) as temp_dir:
        consensus_fasta = "%s/consensus.fa" % temp_dir
        cluster_cmd = "usearch -cluster_fast %s -id 0.97 -strand plus -maxaccepts 0 -maxrejects 0 -rightjust -consout %s -msaout %s/Cluster" % \
                      (read_pairs_fasta.name, consensus_fasta, temp_dir)
        try:
            cluster_output = check_output(cluster_cmd, stderr=subprocess.STDOUT, shell=True)
        except CalledProcessError as cpe:
            pdb.set_trace()
    
        consensus_seqs, ordered_cluster_names = readConsensusSeqs(consensus_fasta)
        umi_classes_per_cluster, read_pair_IDs_per_cluster, orig_read_pair_IDs_per_cluster = getClusterMembers(temp_dir, ordered_cluster_names)
        umi_map_cluster = assignUMIToCluster(consensus_seqs, ordered_cluster_names, umi_classes_per_cluster)
        attemptToMergeSmallAndHighErrorClusters(umi_map_cluster, read_pair_IDs_per_cluster, umi_classes_per_cluster)

    read_pairs_fasta.close()

    return (read_pair_err, read_pair_seqs, read_pair_IDs_per_cluster, orig_read_pair_IDs_per_cluster, umi_map_cluster)


def readConsensusSeqs(consensus_fasta):
    consensus_seqs = {}
    ordered_cluster_names = []   #TODO: currently ordered by 'Cluster#'. Maybe better to be by # seqs in cluster
    
    with open(consensus_fasta, 'r') as ip:
        curr_ID = None
        curr_seq = ""
        for line in ip:
            if (line[0] == '>'):
                if (curr_ID is not None):
                    consensus_seqs[curr_ID] = curr_seq.upper()
                    curr_seq = ""
                curr_ID = line[1:].strip()
                cluster_num = int(curr_ID[7:])
                ordered_cluster_names.append( (cluster_num, curr_ID) )
            else:
                curr_seq += line.strip()

        consensus_seqs[curr_ID] = curr_seq.upper()

    ordered_cluster_names.sort(key=itemgetter(0))
    ordered_cluster_names = list(map(itemgetter(1), ordered_cluster_names))

    return (consensus_seqs, ordered_cluster_names)


def getClusterMembers(temp_dir, cluster_names):
    umi_classes_per_cluster = defaultdict(set)
    read_pair_IDs_per_cluster = defaultdict(set)

    for cluster_name in cluster_names:
        with open("%s/%s" % (temp_dir, cluster_name), 'r') as ip:
            for line in ip:
                if (line[0] == '>'):
                    read_pair_ID = line[1:].strip()
                    read_pair_IDs_per_cluster[cluster_name].add(read_pair_ID)
                    umi_class = tuple(read_pair_ID.split('_')[0:2])
                    umi_classes_per_cluster[cluster_name].add( umi_class )

    orig_read_pair_IDs_per_cluster = {}
    for cluster_name, read_pair_IDs in read_pair_IDs_per_cluster.items():
        orig_read_pair_IDs_per_cluster[cluster_name] = list(read_pair_IDs)

    return (umi_classes_per_cluster, read_pair_IDs_per_cluster, orig_read_pair_IDs_per_cluster)


def assignUMIToCluster(consensus_seqs, ordered_cluster_names, umi_classes_per_cluster):
    umi_classes_already_assigned = set()
    umi_map_cluster = {}

    for cluster_name in ordered_cluster_names:
        # If no UMIs in cluster have already been assigned to another cluster, then assign them to this cluster
        # Else assign all unassigned UMIs to the most-represented cluster of the assigned UMIs

        umis_already_assigned = umi_classes_already_assigned & umi_classes_per_cluster[cluster_name]
        if (len(umis_already_assigned) == 0):
            umi_classes_already_assigned.update(umi_classes_per_cluster[cluster_name])
            for umi in umi_classes_per_cluster[cluster_name]:
                umi_map_cluster[umi] = cluster_name
        else:
            repr_cluster_counts = Counter(map(lambda k: umi_map_cluster[k], umis_already_assigned))
            most_repr_cluster = repr_cluster_counts.most_common(1)[0][0]
            for umi in umi_classes_per_cluster[cluster_name] - umi_classes_already_assigned:
                umi_classes_already_assigned.add(umi)
                umi_map_cluster[umi] = most_repr_cluster

    return umi_map_cluster


def attemptToMergeSmallAndHighErrorClusters(umi_map_cluster, read_pair_IDs_per_cluster, umi_classes_per_cluster):
    cluster_sizes = list(map(lambda k: (k, len(read_pair_IDs_per_cluster[k])), umi_map_cluster.values()))
        
    cluster_sizes.sort(key=itemgetter(1), reverse=True)
    largest_cluster = cluster_sizes[0][0]

    clusters_to_remove = []
    for (cluster, cluster_size) in cluster_sizes[1:]:
        if (cluster_size <= 2):
            read_pair_IDs_per_cluster[largest_cluster] |= read_pair_IDs_per_cluster[cluster]
            umi_classes_per_cluster[largest_cluster] |= umi_classes_per_cluster[cluster]
            for umi_to_remap in map(itemgetter(0), filter(lambda x: x[1]==cluster, umi_map_cluster.items())):
                umi_map_cluster[umi_to_remap] = largest_cluster
            
    for cluster in clusters_to_remove:
        del read_pair_IDs_per_cluster[cluster]
        del umi_classes_per_cluster[cluster]


def computeReadConsensiIndividually(umi_map_cluster, orig_read_pair_IDs_per_cluster, read_pair_seqs, nonunique_umi_class_counts, scratch_dir):
    R1_R2_consensi_map_umis = {}
    
    umis_per_cluster = defaultdict(set)
    for (umi,cluster) in umi_map_cluster.items():
        umis_per_cluster[cluster].add(umi)

    with tempfile.TemporaryDirectory(dir=scratch_dir) as temp_dir:
        for cluster in set(umi_map_cluster.values()):
            #print("%s\t%d" % (cluster, len(umis_per_cluster[cluster])))
            if (len(orig_read_pair_IDs_per_cluster[cluster]) == 1):
                read_pair_ID = orig_read_pair_IDs_per_cluster[cluster][0]
                R1_consensus_seq, R2_consensus_seq = read_pair_seqs[read_pair_ID]
            else:
                R1_fa = open("%s/R1.fa" % temp_dir, 'w')
                R2_fa = open("%s/R2.fa" % temp_dir, 'w')
                R1_consensus = "%s/R1_consensus.fa" % temp_dir
                R2_consensus = "%s/R2_consensus.fa" % temp_dir

                for read_pair_ID in orig_read_pair_IDs_per_cluster[cluster]:
                    R1_seq, R2_seq = read_pair_seqs[read_pair_ID]
                    R1_fa.write(">%s\n%s\n" % (read_pair_ID, R1_seq))
                    R2_fa.write(">%s\n%s\n" % (read_pair_ID, R2_seq))
                R1_fa.close()
                R2_fa.close()

                # Clustering was done on read *pairs*, but now I need to cluster the individual reads. R2 typically has higher error rates,
                # which would be compensated by lower error rates in R1. Here, use low identity threshold to force all reads to be put into one cluster.
                R1_cluster_cmd = "usearch -cluster_fast %s -id 0.75 -strand plus -maxaccepts 0 -maxrejects 0 -consout %s" % (R1_fa.name, R1_consensus)
                R2_cluster_cmd = "usearch -cluster_fast %s -id 0.75 -strand plus -maxaccepts 0 -maxrejects 0 -consout %s" % (R2_fa.name, R2_consensus)

                try:
                    R1_cluster_output = check_output(R1_cluster_cmd, stderr=subprocess.STDOUT, shell=True)
                    R2_cluster_output = check_output(R2_cluster_cmd, stderr=subprocess.STDOUT, shell=True)
                except CalledProcessError as cpe:
                    pdb.set_trace()

                R1_consensus_seqs, R1_ordered_cluster_names = readConsensusSeqs(R1_consensus)
                R2_consensus_seqs, R2_ordered_cluster_names = readConsensusSeqs(R2_consensus)
                
                R1_consensus_seq = R1_consensus_seqs['Cluster0']
                R2_consensus_seq = R2_consensus_seqs['Cluster0']

            sum_nonunique_umi_class_counts = sum(map(lambda umi: nonunique_umi_class_counts[umi], umis_per_cluster[cluster]))
            R1_R2_consensi_map_umis[(R1_consensus_seq, R2_consensus_seq)] = (sum_nonunique_umi_class_counts, len(umis_per_cluster[cluster]))

    return R1_R2_consensi_map_umis
    
        
def printFasta(id_and_seqs, filename):
    with open(filename,'w') as op:
        for seqid, seq in id_and_seqs.items():
            op.write(">%s\n%s\n" % (seqid,seq))


def writeSubgroupResults(R1_R2_consensi_map_umis, primers_signature, triage_result, output_file): # expected_amplicons.bed/tsv
    total_nonunique, total_unique = 0, 0
    lines = []
    for (R1_consensus_seq, R2_consensus_seq), (sum_nonunique_umi_class_counts, sum_unique_umi_class_counts) in R1_R2_consensi_map_umis.items():
        total_nonunique += sum_nonunique_umi_class_counts
        total_unique += sum_unique_umi_class_counts
        lines.append( (R1_consensus_seq, R2_consensus_seq, sum_nonunique_umi_class_counts, sum_unique_umi_class_counts) )
    
    lines.sort(key=itemgetter(2,3), reverse=True)

    op = open(output_file, 'w')
    op.write("%s\t%s\t%d subgroups\t%d nonunique UMIs\t%d unique UMIs\n" % (primers_signature, triage_result, len(lines), total_nonunique, total_unique))
    for subgroup_number, tup in enumerate(lines,1):
        #labeled_tup = (primers_signature, triage_result) + tup
        #op.write("%s\t%s\t%s\t%s\t%d\t%d\n" % labeled_tup)
        labeled_tup = (subgroup_number,) + tup
        op.write("%d\t%s\t%s\t%d\t%d\n" % labeled_tup)
    op.close()

    return (primers_signature, triage_result, len(lines), total_nonunique, total_unique)


def writeSubgroupSummaries(all_subgroup_summaries, subgroup_summary_file):
    with open(subgroup_summary_file, 'w') as op:
        op.write("#PrimersSignature\tTriage\tNumSubgroups\tTotalNonuniqueReads\tTotalUniqueReads\n")
        for tup in all_subgroup_summaries:
            op.write("%s\t%s\t%d\t%d\t%d\n" % tup)


if (__name__ == "__main__"):
    root_scratch_dir, grouping_dir, priming_events_file, subgroup_summary_file = sys.argv[1:]  # TODO: add expected_amplicons.bed/tsv

    with tempfile.TemporaryDirectory(dir=root_scratch_dir) as scratch_dir:
        all_subgroup_summaries = subgroupReadsFromEachPrimingEvent(scratch_dir, grouping_dir, priming_events_file)

    writeSubgroupSummaries(all_subgroup_summaries, subgroup_summary_file)
    sys.exit(0)
