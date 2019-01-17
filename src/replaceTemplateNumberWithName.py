#!/usr/bin/env python3

import sys

if (__name__ == "__main__"):
    template_R1_fasta, template_R2_fasta, input_bed, output_bed = sys.argv[1:]

    number_to_name = {}
    with open(template_R1_fasta, 'r') as ip:
        for line in ip:
            if (line[0] == '>'):
                num, name = line[1:].strip().split()
                number_to_name[num] = name

    with open(template_R1_fasta, 'r') as ip:
        for line in ip:
            if (line[0] == '>'):
                num, name = line[1:].strip().split()
                assert (num in number_to_name and number_to_name[num] == name)

    with open(output_bed,'w') as op:
        with open(input_bed,'r') as ip:
            for line in ip:
                fields = line.strip().split()
                num, R_num = fields[3].split('/')
                assert (num in number_to_name)
                fields[3] = "%s/%s" % (number_to_name[num], R_num)
                op_line = "\t".join(fields)
                op.write("%s\n" % op_line)

    sys.exit(0)
    
