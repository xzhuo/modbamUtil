import os
import sys
import argparse
import time
from mpire import WorkerPool
from itertools import repeat

class Locus:
    def __init__(self, chr, start, end):
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.name = chr + ":" + start + "-" + end
        self.reads = {}

    def add_read(self, read):
        self.reads[read.name] = read

    def remove_reads(self):
        self.reads = {}

    def aggregate_methylation(self, start, length):
        empty_b_cpgs = []
        empty_a_cpgs = []
        insertion_b_cpgs = []
        insertion_a_cpgs = []
        i_cpgs = []
        total_flanking_reads = 0
        total_insertion_reads = 0
        for i in self.reads:
            read = self.reads[i]
            read.methylation_list(start, length)
            if read.readtype == "a,b":
                empty_b_cpgs.extend(read.b_cpgs)
                empty_a_cpgs.extend(read.a_cpgs)
                total_flanking_reads += 1
            elif read.readtype == "a,b,i":
                i_cpgs.extend(read.i_cpgs)
                insertion_b_cpgs.extend(read.b_cpgs)
                insertion_a_cpgs.extend(read.a_cpgs)
                total_flanking_reads += 1
                total_insertion_reads += 1

        self.empty_b_hyper = len([i for i in empty_b_cpgs if i > 0.5])
        self.empty_b_all = len(empty_b_cpgs)
        self.empty_a_hyper = len([i for i in empty_a_cpgs if i > 0.5])
        self.empty_a_all = len(empty_a_cpgs)
        self.with_b_hyper = len([i for i in insertion_b_cpgs if i > 0.5])
        self.with_b_all = len(insertion_b_cpgs)
        self.with_a_hyper = len([i for i in insertion_a_cpgs if i > 0.5])
        self.with_a_all = len(insertion_a_cpgs)
        self.insertion_hyper = len([i for i in i_cpgs if i > 0.5])
        self.insertion_all = len(i_cpgs)
        self.total_flanking_reads = total_flanking_reads
        self.total_insertion_reads = total_insertion_reads
        self.remove_reads()

    def get_read(self, read_name):
        return self.reads[read_name]

    def number_of_reads(self):
        return len(self.reads)

class Read:
    def __init__(self, read, CGins, readtype):
        self.name = read
        self.CGins = CGins
        self.readtype = readtype
        self.cpg = []

    def add_cpg(self, cpg):
        self.cpg.append(cpg)

    def remove_cpg(self):
        self.cpg = []

    def methylation_list(self, start, length):
        end = start + length
        self.b_cpgs = [i.methylation for i in self.cpg if i.rel_pos < start * -1 and i.rel_pos >= end * -1 and i.pos > 0]
        self.a_cpgs = [i.methylation for i in self.cpg if i.rel_pos > start and i.rel_pos <= end and i.pos > 0]
        self.i_cpgs = [i.methylation for i in self.cpg if i.type == "i"]
        # mean_methylation = sum(methylation_list)/len(methylation_list) if len(methylation_list)>0 else -1

class CpG:
    def __init__(self, pos, rel_pos, methylation, type):
        self.pos = int(pos)
        self.rel_pos = int(rel_pos)
        self.methylation = float(methylation)
        self.type = type

def multi_process_aggregate_func(loci_chr, threads, start, length):
    outputs = []
    if threads == 1:
        for i in loci_chr.keys():
            outputs.append(aggregate_func(loci_chr[i], start, length))
    else:
        with WorkerPool(n_jobs=threads) as pool:
            outputs = pool.map(aggregate_func, zip(loci_chr.values(), repeat(start), repeat(length)), iterable_len=len(loci_chr.values()), progress_bar=True)
    return outputs

def aggregate_func(locus, start, length):
    locus.aggregate_methylation(start, length)
    output = [locus.chr, locus.start, locus.end, locus.empty_b_hyper, locus.empty_b_all, locus.empty_a_hyper, locus.empty_a_all, locus.with_b_hyper, locus.with_b_all, locus.with_a_hyper, locus.with_a_all, locus.insertion_hyper, locus.insertion_all, locus.total_flanking_reads, locus.total_insertion_reads]
    return output

def main():
    parser = argparse.ArgumentParser(description='Calculate mean methylation of both flanking and insertion regions of all reads')
    parser.add_argument('-p', '--position_file', type=str, required=True, help='input file with each CpG in all reads')
    parser.add_argument('-g', '--groupby_file', type=str, required=True, help='input file of all reads with read type annotation')
    parser.add_argument('-o', '--output', type=str, required=True, help='output file')
    parser.add_argument('-s', '--start', type=int, default=0, help='position outside of insertion to start averaging methylations')
    parser.add_argument('-l', '--len', type=int, default=500, help='length of flanking regions to average methylation')
    parser.add_argument('-c', '--locus', action=argparse.BooleanOptionalAction, default=False, help='if true, process locus by locus instead of chr by chr')
    parser.add_argument('-t', '--threads', type=int, default=1, help='multi-threading')
    args = parser.parse_args()
    position_file = os.path.abspath(args.position_file)
    groupby_file = os.path.abspath(args.groupby_file)
    if not os.path.exists(position_file):
        raise ValueError("--position file does not exist!")
    if not os.path.exists(groupby_file):
        raise ValueError("--groupby file does not exist!")
    start_time = time.time()
    threads = 1 if args.locus else args.threads  # if locus by locus, no need to use multi-threading
    locus_dict = {}
    """
    groupby file:
    chr1	10861	10862	m64136_200710_174522/11207845/ccs	12	a,b
    chr1	10861	10862	m64136_200710_174522/39323705/ccs	9	a,b
    chr1	10861	10862	m64136_200711_235843/89718828/ccs	82	a,b
    """
    with open(groupby_file, 'r') as f:  # read the region file
        for line in f.readlines():
            chr, start, end, read, CGins, readtype = line.split()
            locus_name = chr + ":" + start + "-" + end
            read_item = Read(read, CGins, readtype)
            if chr not in locus_dict:
                locus_dict[chr] = {}
                locus_dict[chr][locus_name] = Locus(chr, start, end)
                locus_dict[chr][locus_name].add_read(read_item)
            elif locus_name not in locus_dict[chr]:
                locus_dict[chr][locus_name] = Locus(chr, start, end)
                locus_dict[chr][locus_name].add_read(read_item)
            else:
                locus_dict[chr][locus_name].add_read(read_item)
    groupby_time = time.time()
    print("--- It took %s hours to read the groupby file ---" % ((groupby_time - start_time)/3600))

    outputs = []

    """
    position file:
    chr1	10861	10862	10082	m64136_200710_174522/11207845/ccs	-796	0.96	+	b
    chr1	10861	10862	10207	m64136_200710_174522/11207845/ccs	-670	0.98	+	b
    chr1	10861	10862	10468	m64136_200710_174522/11207845/ccs	-418	1.00	+	b
    """
    with open(position_file, 'r') as f:  # read the sorted region file. The file is too big for the memory, so we have to sort it first and then process it locus by locus.
        last_chr = ""
        last_locus = ""
        for line in f.readlines():
            chr, start, end, pos, read, rel_pos, methylation, strand, type = line.split()
            cpg_item = CpG(pos, rel_pos, methylation, type)
            locus = chr + ":" + start + "-" + end
            if chr == last_chr:
                if locus in locus_dict[chr] and read in locus_dict[chr][locus].reads:
                    if args.locus and last_locus != "" and last_locus != locus:
                        print("--- Processing %s ---" % (last_locus))
                        outputs.append(aggregate_func(locus_dict[chr][last_locus], args.start, args.len))
                        del locus_dict[chr][last_locus]
                        last_locus = locus
                    read = locus_dict[chr][locus].get_read(read)
                    read.add_cpg(cpg_item)
            else:
                if last_chr != "":
                    print("--- Processing last locus for %s ---" % (last_chr))
                    outputs.extend(multi_process_aggregate_func(locus_dict[last_chr], threads, args.start, args.len))
                    print("--- Delete %s ---" % (last_chr))
                    del locus_dict[last_chr]
                if locus in locus_dict[chr] and read in locus_dict[chr][locus].reads:
                    read = locus_dict[chr][locus].get_read(read)
                    read.add_cpg(cpg_item)
                last_chr = chr
                last_locus = locus
        print("--- Processing %s ---" % (last_locus))
        outputs.extend(multi_process_aggregate_func(locus_dict[last_chr], threads, args.start, args.len))
        print("--- Delete last chr %s ---" % (last_chr))
        del locus_dict[last_chr]

    # if args.threads == 1: 
    #     for i in locus_dict.keys():
    #         outputs.append(aggregate_func(locus_dict[i], args.len, args.aggregation))
    # else:
    #     with WorkerPool(n_jobs=args.threads) as pool:
    #         outputs = pool.map(aggregate_func, zip(locus_dict.values(), repeat(args.len), repeat(args.aggregation)), iterable_len=len(locus_dict.values()), progress_bar=True)

    with open(args.output, "w") as out:
        out.write("chr\tstart\tend\tempty_before_hyper\tempty_before_all\tempty_after_hyper\tempty_after_all\twith_before_hyper\twith_before_all\twith_after_hyper\twith_after_all\tinsertion_hyper\tinsertion_all\ttotal_flanking_reads\ttotal_insertion_reads\n")
        for i in outputs:
            if len(i) > 0:
                out.write("{:s}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(
                    i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10], i[11], i[12], i[13], i[14]))

    end_time = time.time()
    print("--- In total it took %s hours ---" % ((end_time - start_time)/3600))

if __name__ == '__main__':
    main()