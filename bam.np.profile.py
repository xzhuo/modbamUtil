import os
import argparse
import pysam
import time

def make_table(bam_file, out_file):
    """
    calculate the methylation average of the inserted regions in the bam file
    """

    bam = pysam.AlignmentFile(bam_file, threads = 8, check_sq=False)
    out_list = []
    for read in bam.fetch(until_eof=True):
        if read.is_supplementary or read.is_secondary or read.is_unmapped:
            continue
        query_name = read.query_name
        query_length = read.query_length
        if read.has_tag('np') and read.has_tag('rq') and read.has_tag('ec'):
            np = read.get_tag('np')
            rq = read.get_tag('rq')
            ec = read.get_tag('ec')
            out_list.append([query_name, query_length, np, ec, rq])

    bam.close()
    with open(out_file, "w") as out:
        for line in out_list:
            out.write("{:s}\t{:d}\t{:d}\t{:.4f}\t{:.4f}\n".format(line[0], line[1], line[2], line[3], line[4]))

def main():
    parser = argparse.ArgumentParser(description='calculate CpG methylation average of inserted regions in the bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input revio bam file')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output table file to store np, length, quality')

    args = parser.parse_args()
    bam_file = os.path.abspath(args.bam)
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    start_time = time.time()
    make_table(bam_file, args.out)
    end_time = time.time()
    print("--- %s seconds ---" % (end_time - start_time))


if __name__ == '__main__':
    main()
