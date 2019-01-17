#!/usr/bin/env python3

from collections import defaultdict
from operator import itemgetter
import os
import re
import sys
import tempfile

from groupReadsByPrimerPairUsingAuxData import BED, PrimerAuxData, readPrimerAuxData
from compileSubgroupR1R2 import getPrimerSignaturesMetadata
from matchConcordantSubgroupConsensiToExpected import getBestConsensiMatchToExpected

import pickle


def collectBEDForSubgroupMatchToBestTemplate(subgroup_template_intersect_genome_aligns, pp_aux_data, best_sg_template_matches):
    sg_bed_of_best_template_matches = defaultdict(set)
    template_overlapping_sg_genome_aligns = defaultdict(set)
    
    sg_intersects = defaultdict(list)
    with open(subgroup_template_intersect_genome_aligns, 'r') as ip:
        for line in ip.readlines():
            fields = line.strip().split("\t")
            sg = fields[3][0:-2]
            sg_intersects[sg].append(fields)

    for sg, all_intersects in sg_intersects.items():
        R1_intersect = {}
        R2_intersect = {}
        for fields in all_intersects:
            template, R_number = fields[15].split('/')
            try:
                if (R_number == '1'):
                    sg_R1_BED = BED(fields[0], fields[1], fields[2], fields[3], fields[5], fields[9], fields[10], fields[11])
                    template_R1_BED = BED(fields[12], fields[13], fields[14], fields[15], fields[17], fields[21], fields[22], fields[23])
                    sg_R1_score = int(fields[4])
                    template_R1_score = int(fields[16])
                    R1_overlap_amount = int(fields[-1])

                    if (pp_aux_data[template].hasSameR1BED(template_R1_BED)):
                        if (template not in R1_intersect):
                            R1_intersect[template] = (sg_R1_BED, template_R1_BED, sg_R1_score, template_R1_score, R1_overlap_amount)
                        elif (not sg_R1_BED.isSameAs(R1_intersect[template][0]) or not template_R1_BED.isSameAs(R1_intersect[template][1]) ): # Input intersect file has duplicates
                            assert (R1_intersect[template][4] > R1_overlap_amount or
                                    R1_intersect[template][3] < template_R1_score or
                                    R1_intersect[template][2] < sg_R1_score), "Should intersect less"
                else:
                    sg_R2_BED = BED(fields[0], fields[1], fields[2], fields[3], fields[5], fields[9], fields[10], fields[11])
                    template_R2_BED = BED(fields[12], fields[13], fields[14], fields[15], fields[17], fields[21], fields[22], fields[23])
                    sg_R2_score = int(fields[4])
                    template_R2_score = int(fields[16])
                    R2_overlap_amount = int(fields[-1])

                    if (pp_aux_data[template].hasSameR2BED(template_R2_BED)):
                        if (template not in R2_intersect):
                            R2_intersect[template] = (sg_R2_BED, template_R2_BED, sg_R2_score, template_R2_score, R2_overlap_amount)
                        elif (not sg_R2_BED.isSameAs(R2_intersect[template][0]) or not template_R2_BED.isSameAs(R2_intersect[template][1]) ): # Input intersect file has duplicates
                            assert (R2_intersect[template][4] > R2_overlap_amount or
                                    R2_intersect[template][3] < template_R2_score or
                                    R2_intersect[template][2] < sg_R2_score), "Should intersect less"
            except AssertionError as ae:
                pdb.set_trace()
                print(ae, file=sys.stderr)

                
        common_templates = R1_intersect.keys() & R2_intersect.keys()
        #expected_templates = set(map(lambda t: t[0][0], best_sg_template_matches[sg]))
        #try:
        #    assert (expected_templates.issubset(common_templates))
        #except AssertionError as ae:
        #    pdb.set_trace()

        for template in common_templates:
            sg_R1_BED, template_R1_BED, sg_R1_score, template_R1_score, R1_overlap_amount = R1_intersect[template]
            sg_R2_BED, template_R2_BED, sg_R2_score, template_R2_score, R2_overlap_amount = R2_intersect[template]
                
            sg_R1_label = sg_R1_BED.getLabel()
            sg_R2_label = sg_R2_BED.getLabel()
            template_R1_label = template_R1_BED.getLabel()
            template_R2_label = template_R2_BED.getLabel()
            
            try:
                assert (sg_R1_label[0:-2] == sg_R2_label[0:-2]), "Different subgroup consensus R1 R2 read as a pair"
                assert (sg_R1_label.endswith("/1") and sg_R2_label.endswith("/2")), "Subgroup consensus pair alignments not for both R1 and R2"
                assert (sg == sg_R1_label[0:-2] == sg_R2_label[0:-2])
                assert(template_R1_label[0:-2] == template_R2_label[0:-2]), "Different template R1 R2 read as a pair"
                assert (template_R1_label.endswith("/1") and template_R2_label.endswith("/2")), "Template pair alignments not for both R1 and R2"
            except AssertionError as ae:
                pdb.set_trace()
                print(ae, file=sys.stderr)
                
            sum_sg_scores = sg_R1_score + sg_R2_score
            template_overlapping_sg_genome_aligns[sg].add( (sg_R1_BED, sg_R2_BED, sum_sg_scores, template) )

            if (pp_aux_data[template].hasSameBED(template_R1_BED, template_R2_BED)):
                sg_bed_of_best_template_matches[sg].add( (sg_R1_BED, sg_R2_BED, sum_sg_scores, template) )
                
    return (template_overlapping_sg_genome_aligns, sg_bed_of_best_template_matches)


def collectNewSubgroupConsensusGenomeAlignments(subgroups_genome_aligns_bed, template_overlapping_sg_genome_aligns):
    '''Collect subgroup alignments to genome that are new genomic_alignment_tuples. Add to D[subgroup], with template == "NEW" '''

    nontemplate_overlapping_sg_genome_aligns = defaultdict(set)

    with open(subgroups_genome_aligns_bed, 'r') as ip:
        all_lines = ip.readlines()
    assert (len(all_lines)%2 == 0), "Read an odd number of read pair alignment lines"

    for i in range(0,len(all_lines),2):
        R1_fields = all_lines[i].strip().split("\t")
        R2_fields = all_lines[i+1].strip().split("\t")
        assert(R1_fields[3][0:-2] == R2_fields[3][0:-2]), "Different subgroup consensus R1 R2 read as a pair"
        assert (R1_fields[3].endswith("/1") and R2_fields[3].endswith("/2")), "Pair alignments not for both R1 and R2"
        sg = R1_fields[3][0:-2]

        sg_R1_BED = BED(R1_fields[0], R1_fields[1], R1_fields[2], R1_fields[3], R1_fields[5], R1_fields[9], R1_fields[10], R1_fields[11])
        sg_R2_BED = BED(R2_fields[0], R2_fields[1], R2_fields[2], R2_fields[3], R2_fields[5], R2_fields[9], R2_fields[10], R2_fields[11])
        sum_sg_scores = int(R1_fields[4]) + int(R2_fields[4])

        if (sg in template_overlapping_sg_genome_aligns):
            is_new_alignment = True
            for known_sg_R1R2_BED in template_overlapping_sg_genome_aligns[sg]:
                if (sg_R1_BED.isSameAs(known_sg_R1R2_BED[0]) and sg_R2_BED.isSameAs(known_sg_R1R2_BED[1])):
                    is_new_alignment = False
                    break
            if (is_new_alignment):
                nontemplate_overlapping_sg_genome_aligns[sg].add( (sg_R1_BED, sg_R2_BED, sum_sg_scores, "New") )
        else:
            nontemplate_overlapping_sg_genome_aligns[sg].add( (sg_R1_BED, sg_R2_BED, sum_sg_scores, "New") )

    return nontemplate_overlapping_sg_genome_aligns


def selectGenomeAlignmentForEachSubgroupConsensus(subgroup_counts, template_overlapping_sg_genome_aligns, nontemplate_overlapping_sg_genome_aligns, sg_bed_of_best_template_match):
    align_for_sg = {}
    all_sg = template_overlapping_sg_genome_aligns.keys() | nontemplate_overlapping_sg_genome_aligns.keys()
    for (primer_signature, subgroup_index), (num_nonunique_reads, num_unique_reads) in subgroup_counts.items():
        sg = "%s-%s" % (primer_signature, subgroup_index)
        print(sg, file=sys.stderr)
        if (sg in all_sg):
            assert (len(template_overlapping_sg_genome_aligns[sg])>0 or len(nontemplate_overlapping_sg_genome_aligns[sg])>0)
            template_sg_aligns = list(template_overlapping_sg_genome_aligns[sg])
            nontemplate_sg_aligns = list(nontemplate_overlapping_sg_genome_aligns[sg])

            template_sg_aligns.sort(key=itemgetter(2))
            nontemplate_sg_aligns.sort(key=itemgetter(2))

            if (len(template_sg_aligns)==0):
                best_nontemplate_sg_align = nontemplate_sg_aligns[0]
                align_for_sg[sg] = ("new",) + best_nontemplate_sg_align[0:2]

            elif (len(nontemplate_sg_aligns)==0):
                best_template_sg_align = template_sg_aligns[0]
                if (sg in sg_bed_of_best_template_match):
                    is_new_alignment = True
                    for sg_R1R2_BED in sg_bed_of_best_template_match[sg]:
                        if (best_template_sg_align[0].isSameAs(sg_R1R2_BED[0]) and best_template_sg_align[1].isSameAs(sg_R1R2_BED[1])):
                            template = sg_R1R2_BED[3]
                            align_for_sg[sg] = (template,) + best_template_sg_align[0:2]
                            is_new_alignment = False
                            break
                    if (is_new_alignment):
                        pdb.set_trace()
                        print("A1", file=sys.stderr)
                else:
                    pdb.set_trace()
                    print("A2", file=sys.stderr)
            else:
                best_template_sg_align = template_sg_aligns[0]
                best_nontemplate_sg_align = nontemplate_sg_aligns[0]
                if (best_template_sg_align[2] <= best_nontemplate_sg_align[2]):
                    if (sg in sg_bed_of_best_template_match):
                        is_new_alignment = True
                        for sg_R1R2_BED in sg_bed_of_best_template_match[sg]:
                            if (best_template_sg_align[0].isSameAs(sg_R1R2_BED[0]) and  best_template_sg_align[1].isSameAs(sg_R1R2_BED[1])):
                                template = sg_R1R2_BED[3]
                                align_for_sg[sg] = (template,) + best_template_sg_align[0:2]
                                is_new_alignment = False
                                break
                        if (is_new_alignment):
                            pdb.set_trace()
                            print("B1", file=sys.stderr)
                    else:
                        pdb.set_trace()
                        print("B2", file=sys.stderr)
                else:
                    template = best_nontemplate_sg_align[3]  # Should it be 0, not 3?
                    align_for_sg[sg] = (template,) + best_nontemplate_sg_align[0:2]
        else:
            print("WARNING: no template or genome alignment for %s (%d unique counts), " % (sg, num_unique_reads), file=sys.stderr)

    return align_for_sg


def consolidateSubgroupsAndCounts(align_for_sg, subgroup_counts):
    new_align_for_sg, new_subgroup_counts = {}, {}

    # Identify subgroups that have same genome alignment
    alignment_map_sg = defaultdict(list)
    for sg, alignment in align_for_sg.items():
        alignment_map_sg[alignment].append(sg)

    for alignment, subgroups in alignment_map_sg.items():
        if (len(subgroups) > 1):
            # Use lowest-index subgroup
            index_and_subgroup = list(map(lambda sg: (int(sg.split('-')[1][2:]), sg), subgroups))
            index_and_subgroup.sort(key=itemgetter(0))
            sg_to_use = index_and_subgroup[0][1]

            new_align_for_sg[sg_to_use] = alignment
            
            sum_nonunique_counts = sum(map(lambda sg: subgroup_counts[sg][0], subgroups))
            sum_unique_counts = sum(map(lambda sg: subgroup_counts[sg][1], subgroups))  # This could be inaccurate. Ideally, would need to check UMIs
            new_subgroup_counts[sg_to_use] = (sum_nonunique_counts, sum_unique_counts)
            
        else:
            new_align_for_sg[subgroups[0]] = alignment
            new_subgroup_counts[sg] = subgroup_counts[sg]

    return (new_align_for_sg, new_subgroup_counts)


def standardizeReadPair_DEPRACATED(R1_fields, R2_fields):
    R1_primer_signature, R1_subgroup = R1_fields[3][0:-2].rsplit('-',1)
    R2_primer_signature, R2_subgroup = R2_fields[3][0:-2].rsplit('-',1)
    assert (R1_primer_signature == R2_primer_signature), "R1 and R2 have different primer signatures"
    assert (R1_subgroup == R2_subgroup), "R1 and R2 from different subgroups"
    assert (R1_fields[5] != R2_fields[5]), "R1 and R2 not from opposite strands"
        
    primer_signature = R1_primer_signature
    subgroup = R1_subgroup

    R1_chrom = R1_fields[0]
    R1_start, R1_stop = list(map(int, R1_fields[1:3]))
    R2_chrom = R2_fields[0]
    R2_start, R2_stop = list(map(int, R2_fields[1:3]))
    R1_start += 1  # Convert from BED 0-based to 1-based position
    R2_start += 1

    assert (R1_chrom == R2_chrom)
    chromosome = R1_chrom
    
    strand = R1_fields[5]

    if (strand == '+'):
        assert (R1_start <= R2_start and R2_stop >= R1_stop)
    else:
        assert (strand == '-')
        assert (R1_stop >= R2_stop and R2_start <= R1_start)

    R1_score = int(R1_fields[4])
    R2_score = int(R2_fields[4])
    total_score = R1_score + R2_score

    R1_num_exons = int(R1_fields[9])
    R2_num_exons = int(R2_fields[9])

    R1_exon_lengths = list(map(int, R1_fields[10].split(',')))
    R1_exon_start_offsets = list(map(int, R1_fields[11].split(',')))
    R1_exon_ends = []
    for j in range(R1_num_exons):
        R1_exon_ends.append(R1_start + R1_exon_start_offsets[j])
        R1_exon_ends.append(R1_start + R1_exon_start_offsets[j] + R1_exon_lengths[j] - 1)
    if (strand == "+"):
        R1_exon_ends.pop()
    else:
        R1_exon_ends.pop(0)

    R2_exon_lengths = list(map(int, R2_fields[10].split(',')))
    R2_exon_start_offsets = list(map(int, R2_fields[11].split(',')))
    R2_exon_ends = []
    for j in range(R2_num_exons):
        R2_exon_ends.append(R2_start + R2_exon_start_offsets[j])
        R2_exon_ends.append(R2_start + R2_exon_start_offsets[j] + R2_exon_lengths[j] - 1)
    if (strand == "+"):
        R2_exon_ends.pop(0)
    else:
        R2_exon_ends.pop()
        
    alignment_exon_ends = list(set(R1_exon_ends + R2_exon_ends))
    alignment_exon_ends.sort()

    if (strand == '-'):
        assert (R1_stop in alignment_exon_ends and R2_start in alignment_exon_ends)
    else:
        assert (R1_start in alignment_exon_ends and R2_stop in alignment_exon_ends)

    return (primer_signature, subgroup, chromosome, tuple(alignment_exon_ends), total_score)


def readAndStandardizeAlignments_DEPRACATED(temp_input):
    pair_alignments = defaultdict(list)

    with open(temp_input, 'r') as op:
        all_lines = op.readlines()
    assert (len(all_lines)%2 == 0), "Read an odd number of read pair alignment lines"

    for i in range(0,len(all_lines),2):
        R1_fields = all_lines[i].strip().split("\t")
        R2_fields = all_lines[i+1].strip().split("\t")
        assert (R1_fields[3].endswith("/1") and R2_fields[3].endswith("/2")), "Pair alignments not for both R1 and R2"
        primer_signature, subgroup, chromosome, alignment_exon_ends, total_score = standardizeReadPair(R1_fields, R2_fields)
        pair_alignments[(primer_signature, subgroup)].append( (chromosome, alignment_exon_ends, total_score) )

    return pair_alignments


def collectSubgroupData(dir_for_grouping, primer_signatures_meta):
    subgroup_counts = {}
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
                sg_tup = (primer_signature, subgroup_index)
                subgroup_counts[sg_tup] = (int(num_nonunique_reads), int(num_unique_reads))
            assert (subgroup_index == int(num_subgroups)), "Did not find expected number of subgroups for %s" % subgroup

    return subgroup_counts


def reportAlignmentsAndCounts(align_for_sg, subgroup_counts, output_file):
    write_tuples = []
    for sg, (template, R1_bed, R2_bed) in align_for_sg.items():
        primer_signature, subgroup = sg.rsplit('-',1)
        subgroup_index = subgroup[2:]
        num_nonunique_reads, num_unique_reads = subgroup_counts[(primer_signature, subgroup_index)]
        write_tuples.append( (primer_signature, subgroup, template, str(R1_bed), str(R2_bed), num_nonunique_reads, num_unique_reads) )

    write_tuples.sort(key=itemgetter(6), reverse=True)

    op = open(output_file, 'w')
    for tup in write_tuples:
        op.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\n" % tup)
    op.close()
    

if (__name__ == "__main__"):
    root_scratch_dir, dir_for_grouping, pp_aux_data_file, primer_signatures_summary_file, subgroup_consensi_R1, subgroup_consensi_R2, \
        subgroups_genome_aligns_bed, subgroup_template_intersect_genome_aligns, output_file = sys.argv[1:]

    pdb.set_trace()
    pp_aux_data = readPrimerAuxData(pp_aux_data_file)
    
    with tempfile.TemporaryDirectory(dir=root_scratch_dir) as tempdir:
        best_sg_template_matches = getBestConsensiMatchToExpected(tempdir, pp_aux_data, subgroup_consensi_R1, subgroup_consensi_R2)

    template_overlapping_sg_genome_aligns, sg_bed_of_best_template_matches = \
                collectBEDForSubgroupMatchToBestTemplate(subgroup_template_intersect_genome_aligns, pp_aux_data, best_sg_template_matches)
    nontemplate_overlapping_sg_genome_aligns = collectNewSubgroupConsensusGenomeAlignments(subgroups_genome_aligns_bed, template_overlapping_sg_genome_aligns)

    # TODO: 1) Determine the best BED (template or genomic) for each subgroup
    #       2) Combine subgroups that have the same best BED
    primer_signatures_meta = getPrimerSignaturesMetadata(primer_signatures_summary_file)
    subgroup_counts = collectSubgroupData(dir_for_grouping, primer_signatures_meta)

    #pickle.dump((subgroup_counts, template_overlapping_sg_genome_aligns, nontemplate_overlapping_sg_genome_aligns, sg_bed_of_best_template_matches), open("data.pkl", 'wb'))
    #subgroup_counts, template_overlapping_sg_genome_aligns, nontemplate_overlapping_sg_genome_aligns, sg_bed_of_best_template_matches = pickle.load(open("data.pkl", 'rb'))
    align_for_sg = selectGenomeAlignmentForEachSubgroupConsensus(subgroup_counts, template_overlapping_sg_genome_aligns, nontemplate_overlapping_sg_genome_aligns,
                                                                 sg_bed_of_best_template_matches)

    new_align_for_sg, new_subgroup_counts = consolidateSubgroupsAndCounts(align_for_sg, subgroup_counts)

    # Select best matching template or genomic alignment for each subgroup (use /raid1/projects/NewIsoforms/scripts/selectBestAlignments.py)
    # Combine subgroups with same best, merging counts
    # Report for each subgroup
    reportAlignmentsAndCounts(new_align_for_sg, new_subgroup_counts, output_file)

    sys.exit(0)
