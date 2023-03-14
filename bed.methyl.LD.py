import os
import argparse
import time
from scipy import stats
from mpire import WorkerPool
from itertools import repeat

def process_line(line):
    line_list = line.strip().split("\t")
    chrom=line_list[0]
    pos=line_list[1]
    try:
        methylation_set=set(line_list[2].split(","))
    except IndexError:
        methylation_set=set()
    try:
        unmethylation_set=set(line_list[3].split(","))
    except IndexError:
        unmethylation_set=set()

    return chrom, pos, methylation_set, unmethylation_set

def process_input(input_file):
    out_list = []
    with open(input_file, 'r') as f:  # read the chrom size file
        last_line = ""
        for line in f.readlines():
            if last_line != "":
                chrom, pos, met, unmet = process_line(line)
                last_chrom, last_pos, last_met, last_unmet = process_line(last_line)
                if chrom == last_chrom:
                    distance = int(pos) - int(last_pos)
                    both_met = len(met.intersection(last_met))
                    both_unmet = len(unmet.intersection(last_unmet))
                    unmet_met = len(met.intersection(last_unmet))
                    met_unmet = len(unmet.intersection(last_met))
                    fisher_ratio, fisher_p = stats.fisher_exact(table=[[both_met,unmet_met],[met_unmet,both_unmet]], alternative="greater")
                    out_list.append([chrom, pos, distance, fisher_ratio, fisher_p])
            last_line = line

    return out_list


def main():
    parser = argparse.ArgumentParser(description='parser the pos bed like file to calculate methylation epiallele using LD or fisher exact test')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='multi-threading')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='input bed like file with pos and read ID separated by methylation calls')
    parser.add_argument('-w', '--window', type=int, default=10000000,
                        help='processing window chunk size')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output txt file storing the methylation epiallele linkage')

    args = parser.parse_args()
    input_file = os.path.abspath(args.input)
    if not os.path.exists(input_file):
        raise ValueError("--input bed like file does not exist!")


    start_time = time.time()

    out_list = process_input(input_file)
    with open(args.out, "w") as out:
        for line in out_list:
            out.write("{:s}\t{:d}\t{:d}\t{:s}\t{:s}\n".format(
                line[0], line[1], line[2], line[3], line[4]))
    end_time = time.time()
    print("--- %s hours ---" % ((end_time - start_time)/3600))

if __name__ == '__main__':
    main()