import os
import argparse
import pysam
import time


def separate_reads(bam_file, gfa_dict, prefix):
    """
    separate reads by haplotype
    """
    bam = pysam.AlignmentFile(bam_file, threads=8, check_sq=False)
    pat_reads = pysam.AlignmentFile(prefix + ".pat.bam", "wb", template=bam)
    mat_reads = pysam.AlignmentFile(prefix + ".mat.bam", "wb", template=bam)
    amb_reads = pysam.AlignmentFile(prefix + ".amb.bam", "wb", template=bam)
    un_reads = pysam.AlignmentFile(prefix + ".un.bam", "wb", template=bam)
    for read in bam.fetch(until_eof=True):
        if read.query_name in gfa_dict['p']:
            pat_reads.write(read)
        elif read.query_name in gfa_dict['m']:
            mat_reads.write(read)
        elif read.query_name in gfa_dict['a']:
            amb_reads.write(read)
        else:
            un_reads.write(read)


def main():
    parser = argparse.ArgumentParser(description='Separate reads from a bam file by haplotype using a hifiasm gfa file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file with Mm and Ml tags')
    parser.add_argument('-g', '--gfa', type=str, required=True,
                        help='a hifiasm gfa file containing the haplotype information')


    args = parser.parse_args()
    prefix = args.bam[:-4]
    bam_file = os.path.abspath(args.bam)
    gfa_file = os.path.abspath(args.gfa)
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    if not os.path.exists(gfa_file):
        raise ValueError("--gfa file does not exist!")
    start_time = time.time()
    gfa_dict = {}
    with open(gfa_file, 'r') as f:  # read the gfa file
        for line in f.readlines():
            line_list = line.strip().split('\t')
            if line_list[0] == 'A':
                hap = line_list[8][-1:]
                if hap not in gfa_dict:
                    gfa_dict[hap] = [line_list[4]]
                else:
                    gfa_dict[hap].append(line_list[4])

    separate_reads(bam_file, gfa_dict, prefix)

    end_time = time.time()
    print("--- %s hours ---" % ((end_time - start_time)/3600))

if __name__ == '__main__':
    main()