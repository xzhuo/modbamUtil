import os
import argparse
import time
from scipy import stats
import math
from mpire import WorkerPool
from itertools import repeat

def process_line(line):
    line_list = line.strip().split("\t")
    chrom=line_list[0]
    pos=line_list[1]
    try:
        methylation_set=set(line_list[3].split(","))
    except IndexError:
        methylation_set=set()
    try:
        unmethylation_set=set(line_list[4].split(","))
    except IndexError:
        unmethylation_set=set()

    return chrom, pos, methylation_set, unmethylation_set

# def process_input(input_file):
#     out_list = []
#     with open(input_file, 'r') as f:  # read the chrom size file
#         last_line = ""
#         for line in f.readlines():
#             if last_line != "":
#                 chrom, pos, met, unmet = process_line(line)
#                 last_chrom, last_pos, last_met, last_unmet = process_line(last_line)
#                 if chrom == last_chrom:
#                     distance = int(pos) - int(last_pos)
#                     both_met = len(met.intersection(last_met))
#                     both_unmet = len(unmet.intersection(last_unmet))
#                     unmet_met = len(met.intersection(last_unmet))
#                     met_unmet = len(unmet.intersection(last_met))
#                     depth = len(met) + len(unmet)
#                     met_perc = len(met)/depth if depth > 0 else float("NaN")
#                     fisher_ratio, fisher_p = stats.fisher_exact(table=[[both_met,unmet_met],[met_unmet,both_unmet]], alternative="greater")
#                     log_fisher_p = -math.log10(fisher_p) + 0 if fisher_p > 0 else float("inf")  # get the log10 p value and add 0 to get rid of the -0.0 or return inf if p value is 0.
#                     out_list.append([chrom, int(pos), distance, depth, met_perc, fisher_ratio, log_fisher_p])
#             last_line = line

#     return out_list

# def process_lines(last_line, line):
#     out_list = []
#     chrom, pos, met, unmet = process_line(line)
#     last_chrom, last_pos, last_met, last_unmet = process_line(last_line)
#     if chrom == last_chrom:
#         distance = int(pos) - int(last_pos)
#         both_met = len(met.intersection(last_met))
#         both_unmet = len(unmet.intersection(last_unmet))
#         unmet_met = len(met.intersection(last_unmet))
#         met_unmet = len(unmet.intersection(last_met))
#         depth = len(met) + len(unmet)
#         met_perc = len(met)/depth if depth > 0 else float("NaN")
#         fisher_ratio, fisher_p = stats.fisher_exact(table=[[both_met,unmet_met],[met_unmet,both_unmet]], alternative="greater")
#         log_fisher_p = -math.log10(fisher_p) + 0 if fisher_p > 0 else float("inf")  # get the log10 p value and add 0 to get rid of the -0.0 or return inf if p value is 0.
#         out_list = [chrom, int(pos), distance, depth, met_perc, fisher_ratio, log_fisher_p]

#     return out_list

def process_arrays(array, i):
    out_list = []
    line = array[i]
    last_line = array[i-1]
    chrom, pos, met, unmet = process_line(line)
    last_chrom, last_pos, last_met, last_unmet = process_line(last_line)
    if chrom == last_chrom:
        distance = int(pos) - int(last_pos)
        both_met = len(met.intersection(last_met))
        both_unmet = len(unmet.intersection(last_unmet))
        unmet_met = len(met.intersection(last_unmet))
        met_unmet = len(unmet.intersection(last_met))
        depth = len(met) + len(unmet)
        met_perc = len(met)/depth if depth > 0 else float("NaN")
        fisher_ratio, fisher_p = stats.fisher_exact(table=[[both_met,unmet_met],[met_unmet,both_unmet]], alternative="greater")
        log_fisher_p = -math.log10(fisher_p) + 0 if fisher_p > 0 else float("inf")  # get the log10 p value and add 0 to get rid of the -0.0 or return inf if p value is 0.
        out_list = [chrom, int(pos), distance, depth, met_perc, fisher_ratio, log_fisher_p]

    return out_list

def main():
    parser = argparse.ArgumentParser(description='parser the pos bed like file to calculate methylation epiallele using LD or fisher exact test')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='multi-threading')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='input bed like file with pos and read ID separated by methylation calls')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output txt file storing the methylation epiallele linkage')

    args = parser.parse_args()
    input_file = os.path.abspath(args.input)
    if not os.path.exists(input_file):
        raise ValueError("--input bed like file does not exist!")


    start_time = time.time()
    infile_array = []
    out_list = []
    with open(input_file, 'r') as f:  # read the bed file
        for line in f.readlines():
            infile_array.append(line.strip())

    # out_list = process_input(input_file)

    with WorkerPool(n_jobs=args.threads, shared_objects=infile_array, start_method='fork') as pool:
        out_list = pool.imap(process_arrays, range(1, len(infile_array)), iterable_len=len(infile_array) - 1, progress_bar=True)
    with open(args.out, "w") as out:
        for line in out_list:
            if line == []:
                continue
            out.write("{:s}\t{:d}\t{:d}\t{:d}\t{:0.2f}\t{:0.2f}\t{:0.4f}\n".format(
                line[0], line[1], line[2], line[3], line[4], line[5], line[6]))
    end_time = time.time()
    print("--- %s hours ---" % ((end_time - start_time)/3600))

if __name__ == '__main__':
    main()