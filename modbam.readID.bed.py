import os
from collections import defaultdict
from concurrent.futures import process
import pysam
from modbampy import ModBam
import argparse
import time
from mpire import WorkerPool
from itertools import repeat

def process_bam(bam_file, window_dict, merge):
    chrom = window_dict["chrom"]
    start = window_dict["start"]
    end = window_dict["end"]
    outfile_list = []
    with ModBam(bam_file) as bam:
        nested_dict = lambda: defaultdict(nested_dict)
        output = nested_dict()
        for read in bam.reads(chrom, start, end):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
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
                if pos > start and pos <= end:
                    strand = pos_mod[3]
                    if merge:
                        pos = pos if strand == "+" else pos - 1
                        strand = "."
                    if pos not in output[chrom]:
                        output[chrom][pos] = {"strand": strand, "methylated": [], "unmethylated": []}

                    if pos_mod[7]/255 > 0.5:
                        output[chrom][pos]["methylated"].append(pos_mod[0])
                    else:
                        output[chrom][pos]["unmethylated"].append(pos_mod[0])

    for chrom in output:
        for pos in output[chrom]:
            outfile_list.append([chrom, pos, output[chrom][pos]["strand"], ",".join(output[chrom][pos]["methylated"]), ",".join(output[chrom][pos]["unmethylated"])])

    return outfile_list

def process_chromsize(bam_file, window):
    size_list = []
    bam = pysam.AlignmentFile(bam_file)
    chrom_sizes = zip(bam.references, bam.lengths)
    for chrom, size in chrom_sizes:
        for i in range(0,size,window):
            if i > 0:
                size_list.append({"chrom": chrom, "start": last_i, "end": i})
            last_i = i
        size_list.append({"chrom": chrom, "start": last_i, "end": size})
    # with open(chromsize_file, 'r') as f:  # read the chrom size file
    #     for line in f.readlines():
    #         chrom_list = line.strip().split()
    #         for i in range(0,int(chrom_list[1]),window):
    #             if i > 0:
    #                 size_list.append({"chrom": chrom_list[0], "start": last_i, "end": i})
    #             last_i = i
    #         size_list.append({"chrom": chrom_list[0], "start": last_i, "end": int(chrom_list[1])})
    return size_list

def main():
    parser = argparse.ArgumentParser(description='parser the bam file and summarize which reads are methylated/unmethylated at each CpG site')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='multi-threading')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file with Mm and Ml tags')
    parser.add_argument('-w', '--window', type=int, default=10000000,
                        help='processing window size')
    parser.add_argument('-m', '--merge', action=argparse.BooleanOptionalAction, default=True,
                        help='merge both strand or not')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output bed like txt file storing the methylation data in the defined regions')

    args = parser.parse_args()
    bam_file = os.path.abspath(args.bam)
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")

    start_time = time.time()
    size_list = process_chromsize(bam_file, args.window)

    outputs = []
    if args.threads == 1:
        for i in size_list:
            outputs.append(process_bam(bam_file, i, args.merge))
    else:
        with WorkerPool(n_jobs=args.threads) as pool:
            outputs = pool.map(process_bam, zip(repeat(bam_file), size_list, repeat(args.merge)), iterable_len=len(size_list), progress_bar=True)

    # merge the adjacent windows if they are of the same coordinate:
    for i, batch in enumerate(outputs):
        if i > 0 and len(batch) > 0 and len(outputs[i-1]) > 0:
            if outputs[i-1][-1][0] == batch[0][0] and outputs[i-1][-1][1] == batch[0][1]:
                outputs[i-1][-1][3] = outputs[i-1][-1][3] + "," + batch[0][3]
                outputs[i-1][-1][4] = outputs[i-1][-1][4] + "," + batch[0][4]
                outputs[i] = outputs[i][1:]

    flat_outputs = [item for batch in outputs for item in batch]
    flat_outputs.sort(key=lambda x: (x[0],x[1]))
    with open(args.out, "w") as out:
        for line in flat_outputs:
                if line[1] >= 0:
                    out.write("{:s}\t{:d}\t{:s}\t{:s}\t{:s}\n".format(
                        line[0], line[1], line[2], line[3], line[4]))

    end_time = time.time()
    print("--- %s hours ---" % ((end_time - start_time)/3600))

if __name__ == '__main__':
    main()