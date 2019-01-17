#!/usr/bin/env python3

from collections import defaultdict
from itertools import groupby
from operator import itemgetter
import sys
import tempfile

from groupReadsByPrimerPairUsingAuxData import BED, PrimerAuxData, readPrimerAuxData

import subprocess
from subprocess import check_output, CalledProcessError


def getBestConsensiMatchToExpected(tempdir, pp_aux_data, subgroup_consensi_R1, subgroup_consensi_R2):
    all_R1 = defaultdict(list)
    all_R2 = defaultdict(list)
    best_R1 = defaultdict(lambda: None)
    best_R2 = defaultdict(lambda: None)

    expected_R1_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=tempdir, delete=True)
    expected_R2_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=tempdir, delete=True)
    for primer_signature_and_isoform_IDs, aux_data in pp_aux_data.items():
        isoform_IDs = aux_data.getIsoformIDs()
        R1_seq, R2_seq = aux_data.getReadPair()
        expected_R1_fasta.write(">%s\n%s\n" % (primer_signature_and_isoform_IDs, R1_seq))
        expected_R2_fasta.write(">%s\n%s\n" % (primer_signature_and_isoform_IDs, R2_seq))
    expected_R1_fasta.flush()
    expected_R1_fasta.flush()

    ip_usearch_R1 = tempfile.NamedTemporaryFile(mode='w+t', suffix=".usearch", dir=tempdir, delete=True)
    usearch_R1_cmd = ['usearch', '-search_global', subgroup_consensi_R1, '-db', expected_R1_fasta.name, '-strand', 'plus', '-id', '0.85', '-leftjust', '-rightjust',
                      '-userout', ip_usearch_R1.name, '-userfields', 'query+qlo+qhi+ql+qstrand+target+tlo+thi+tl+tstrand+alnlen+mism+ids+opens+exts+diffs+id+aln']

    ip_usearch_R2 = tempfile.NamedTemporaryFile(mode='w+t', suffix=".usearch", dir=tempdir, delete=True)    
    usearch_R2_cmd = ['usearch', '-search_global', subgroup_consensi_R2, '-db', expected_R2_fasta.name, '-strand', 'plus', '-id', '0.85', '-leftjust', '-rightjust',
                      '-userout', ip_usearch_R2.name, '-userfields', 'query+qlo+qhi+ql+qstrand+target+tlo+thi+tl+tstrand+alnlen+mism+ids+opens+exts+diffs+id+aln']

    try:
        usearch_R1_output = check_output(usearch_R1_cmd, stderr=subprocess.STDOUT)
        ip_usearch_R1.flush()
        ip_usearch_R1.seek(0)
    except CalledProcessError as cpe:
        pdb.set_trace()

    for line in ip_usearch_R1:
        query, qlo, qhi, ql, qstrand, template, tlo, thi, tl, tstrand, alnlen, mism, ids, opens, exts, diffs, frac_id, cigar = \
            map(lambda z: z[0](z[1]), zip((str,int,int,int,str,str,int,int,int,str,int,int,int,int,int, int,float,str),line.strip().split('\t')))

        num_gap_pos = opens + exts
        IorD_cigar = list(filter(lambda x:x[0]!='M', [(k,len(list(g))) for k, g in groupby(cigar)]))
        max_IorD_length = 0 if (len(IorD_cigar)==0) else max(map(itemgetter(1), IorD_cigar))
        edit_dist = num_gap_pos + mism
        if (max_IorD_length <= 2 and mism <= int(tl/100)*7 and num_gap_pos <= 4):
            all_R1[query].append( (template, edit_dist, frac_id, ids, mism, max_IorD_length, opens, num_gap_pos) )                

    try:
        usearch_R2_output = check_output(usearch_R2_cmd, stderr=subprocess.STDOUT)
        ip_usearch_R2.flush()
        ip_usearch_R2.seek(0)
    except CalledProcessError as cpe:
        pdb.set_trace()

    for line in ip_usearch_R2:
        query, qlo, qhi, ql, qstrand, template, tlo, thi, tl, tstrand, alnlen, mism, ids, opens, exts, diffs, frac_id, cigar = \
            map(lambda z: z[0](z[1]), zip((str,int,int,int,str,str,int,int,int,str,int,int,int,int,int,int,float,str),line.strip().split('\t')))

        num_gap_pos = opens + exts
        IorD_cigar = list(filter(lambda x:x[0]!='M', [(k,len(list(g))) for k, g in groupby(cigar)]))
        max_IorD_length = 0 if (len(IorD_cigar)==0) else max(map(itemgetter(1), IorD_cigar))
        edit_dist = num_gap_pos + mism
        if (max_IorD_length <= 2 and mism <= int(tl/100)*9 and num_gap_pos <= 5): # mism and gap criteria larger for R2 because R2 is noisier
            all_R2[query].append( (template, edit_dist, frac_id, ids, mism, max_IorD_length, opens, num_gap_pos) )

    ip_usearch_R1.close()
    ip_usearch_R2.close()

    expected_R1_fasta.close()
    expected_R2_fasta.close()

    # Retain the subgroups whose best R1 & R2 expected-template matches are the same
    best_sg_consensus_to_template_matches = defaultdict(list)
    for sg in all_R1.keys() & all_R2.keys():
        all_R1_for_sg = all_R1[sg]
        all_R1_for_sg.sort(key=itemgetter(1))
        smallest_edit_dist = all_R1_for_sg[0][1]
        best_R1_for_sg = list(filter(lambda x: x[1]==smallest_edit_dist, all_R1_for_sg))

        all_R2_for_sg = all_R2[sg]
        all_R2_for_sg.sort(key=itemgetter(1))
        smallest_edit_dist = all_R2_for_sg[0][1]
        best_R2_for_sg = list(filter(lambda x: x[1]==smallest_edit_dist, all_R2_for_sg))

        for common_best_template in set(map(itemgetter(0), best_R1_for_sg)) & set(map(itemgetter(0), best_R2_for_sg)):
            R1_template_matches = list(filter(lambda x:x[0]==common_best_template, all_R1[sg]))
            assert (len(R1_template_matches)==1)
            R1_template_match = list(R1_template_matches)[0]
            R2_template_matches = list(filter(lambda x:x[0]==common_best_template, all_R2[sg]))
            assert (len(R2_template_matches)==1)
            R2_template_match = list(R2_template_matches)[0]
            best_sg_consensus_to_template_matches[sg].append( (R1_template_match, R2_template_match) )

    return best_sg_consensus_to_template_matches

                
def writeConsensusMatches(best_sg_consensus_to_template_matches, matched_concordant_output):
    with open(matched_concordant_output, 'w') as op:
        for subgroup_consensus, best_R12_template_matches in best_sg_consensus_to_template_matches.items():
            for (best_R1_template_match, best_R2_template_match) in best_R12_template_matches:
                (target, edit_dist, frac_id, ids, mism, max_IorD_length, opens, num_gap_pos) = best_R1_template_match # TODO: for now
                primer_signature, target_isoform_group = target.split('+')
                op.write("%s\t%s\t%s\t%d\t%3.2f\t%d\t%d\t%d\t%d\t%d\n" % \
                         (subgroup_consensus, primer_signature, target_isoform_group, edit_dist, frac_id, ids, mism, max_IorD_length, opens, num_gap_pos))


if (__name__ == "__main__"):
    root_scratch_dir, pp_aux_data_file, subgroup_consensi_R1, subgroup_consensi_R2, matched_concordant_output = sys.argv[1:]

    pp_aux_data = readPrimerAuxData(pp_aux_data_file)

    with tempfile.TemporaryDirectory(dir=root_scratch_dir) as tempdir:
        best_sg_consensus_to_template_matches = getBestConsensiMatchToExpected(tempdir, pp_aux_data, subgroup_consensi_R1, subgroup_consensi_R2)

    writeConsensusMatches(best_sg_consensus_to_template_matches, matched_concordant_output)

    sys.exit(0)
