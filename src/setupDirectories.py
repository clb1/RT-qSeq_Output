#!/usr/bin/env python3

import os
import sys

if (__name__ == "__main__"):
    setup_data_file, root_directory = sys.argv[1:]

    ip = open(setup_data_file, 'r')

    for line in ip:
        if (line[0] == '#'):
            continue
        
        dest_dir, R1, R2, fasta, aux = line.strip().split()
        assert("R1" in R1 and "R2" not in R1)
        assert("R2" in R2 and "R1" not in R2)
        assert(fasta.endswith(".fa"))
        assert(aux.endswith("_aux.tsv"))

        full_dest_dir = "%s/%s" % (root_directory, dest_dir)
        assert(os.path.exists(full_dest_dir))
        makefile_path = "%s/Makefile" % full_dest_dir

        with open(makefile_path, 'w') as op:
            op.write("R1_FASTQ=%s\n" % R1)
            op.write("R2_FASTQ=%s\n" % R2)
            op.write("SUBPOOL_NAMES = 1\n")
            op.write("PRIMERS_FASTA = %s\n" % fasta)
            op.write("PRIMERS_AUX_DATA = %s\n" % aux)
            op.write("include ../Makefile.common\n")
        
    ip.close()

    sys.exit(0)
