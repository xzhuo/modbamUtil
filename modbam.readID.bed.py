import os
from collections import defaultdict
from concurrent.futures import process
from modbampy import ModBam
import argparse
import time
from mpire import WorkerPool
from itertools import repeat

def process_bam(bam_file, window_dict):
    chrom = window_dict["chrom"]
    start = window_dict["start"]
    end = window_dict["end"]
    outfile_list = []
    with ModBam(bam_file) as bam:
        nested_dict = lambda: defaultdict(nested_dict)
        output = nested_dict()
        for read in bam.reads(chrom, start, end):
            # print(read.mod_sites)
            for pos_mod in read.mod_sites:
                """read_id,
                reference position,
                query (read) position,
                reference strand (+ or -),
                modification strand (0 or 1, as defined in the HTSlib tag specification. This is invariable 0),
                canonical base associated with modification,
                modified base,
                modified-base score (scaled to 0-255)."""
                pos = pos_mod[1]
                strand = pos_mod[3]
                if pos not in output[chrom]:
                    output[chrom][pos] = {"strand": strand, "methylated": [], "unmethylated": []}

                if pos_mod[1] > start and pos_mod[1] <= end:
                    if pos_mod[7]/255 > 0.5:
                        output[chrom][pos]["methylated"].append(pos_mod[0])
                    else:
                        output[chrom][pos]["unmethylated"].append(pos_mod[0])

    for chrom in output:
        for pos in output[chrom]:
            outfile_list.append(chrom, pos, output[chrom][pos]["strand"], ",".join(output[chrom][pos]["methylated"]), ",".join(output[chrom][pos]["unmethylated"]))

    return outfile_list

def process_chromsize(chromsize_file):
    size_list = []
    with open(chromsize_file, 'r') as f:  # read the chrom size file
        for line in f.readlines():
            chrom_list = line.strip().split()
            chrom_dict = {"chrom": chrom_list[0], "start": 0, "end": int(chrom_list[1])}
            size_list.append(chrom_dict)
    return size_list

def main():
    parser = argparse.ArgumentParser(description='the bam file')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='multi-threading')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file with Mm and Ml tags')
    parser.add_argument('-c', '--chrom', type=str, required=True,
                        help='input chromsize file')
    parser.add_argument('-w', '--window', type=int, default=100000000,
                        help='processing window size')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output bed like txt file storing the methylation data in the defined regions')

    args = parser.parse_args()
    bam_file = os.path.abspath(args.bam)
    chromsize_file = os.path.abspath(args.chrom)
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    if not os.path.exists(chromsize_file):
        raise ValueError("--chromsize file does not exist!")

    start_time = time.time()
    size_list = process_chromsize(chromsize_file)

    outputs = []
    if args.threads == 1:
        for i in size_list:
            outputs.append(process_bam(bam_file, i))
    else:
        with WorkerPool(n_jobs=args.threads) as pool:
            outputs = pool.imap(process_bam, zip(repeat(bam_file), size_list), iterable_len=len(size_list), progress_bar=True)

    outputs.sort(key=lambda x: (x[0],x[1]))
    with open(args.out, "w") as out:
        for line in outputs:
                out.write("{:s}\t{:d}\t{:s}\t{:s}\t{:s}\n".format(
                    line[0], line[1], line[2], line[3], line[4]))

    end_time = time.time()
    print("--- %s hours ---" % ((end_time - start_time)/3600))
