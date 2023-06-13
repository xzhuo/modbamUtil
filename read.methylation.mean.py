import os
import sys
import argparse
import time
from mpire import WorkerPool

class Locus:
    def __init__(self, chr, start, end):
        self.chr = chr
        self.start = start
        self.end = end
        self.name = chr + ":" + start + "-" + end
        self.reads = {}

    def add_read(self, read):
        self.reads[read.name] = read

    def aggregate_methylation(self, length, aggregate_type = "mean"):
        empty_b_cpgs = []
        empty_a_cpgs = []
        insertion_b_cpgs = []
        insertion_a_cpgs = []
        i_cpgs = []
        total_flanking_reads = 0
        total_insertion_reads = 0
        for read in self.reads:
            read.methylation_list(length)
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

        if aggregate_type == "mean":
            if len(empty_b_cpgs) > 0 and len(empty_a_cpgs) > 0 and len(insertion_b_cpgs) > 0 and len(insertion_a_cpgs) > 0 and len(i_cpgs) > 0:
                self.empty_b_methylation = sum(empty_b_cpgs)/len(empty_b_cpgs)
                self.empty_a_methylation = sum(empty_a_cpgs)/len(empty_a_cpgs)
                self.insertion_b_methylation = sum(insertion_b_cpgs)/len(insertion_b_cpgs)
                self.insertion_a_methylation = sum(insertion_a_cpgs)/len(insertion_a_cpgs)
                self.i_methylation = sum(i_cpgs)/len(i_cpgs)
        elif aggregate_type == "count":
            if total_flanking_reads > 0 and total_insertion_reads > 0:
                self.empty_b_methylation = len([i for i in empty_b_cpgs if i > 0.5])/total_flanking_reads
                self.empty_a_methylation = len([i for i in empty_a_cpgs if i > 0.5])/total_flanking_reads
                self.insertion_b_methylation = len([i for i in insertion_b_cpgs if i > 0.5])/total_insertion_reads
                self.insertion_a_methylation = len([i for i in insertion_a_cpgs if i > 0.5])/total_insertion_reads
                self.i_methylation = sum(i_cpgs)/len(i_cpgs)  # always return mean methylation for insertion region.
        else:
            raise ValueError("aggregate_type can only be mean or count")

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

    def methylation_list(self,length):
        self.b_cpgs = [i["methylation"] for i in self.cpg if i["rel_pos"]<0 and i["rel_pos"]>length * -1 and i["pos"] > 0]
        self.a_cpgs = [i["methylation"] for i in self.cpg if i["rel_pos"]>0 and i["rel_pos"]<length and i["pos"] > 0]
        self.i_cpgs = [i["methylation"] for i in self.cpg if i["type"] == "i"]
        # mean_methylation = sum(methylation_list)/len(methylation_list) if len(methylation_list)>0 else -1

class CpG:
    def __init__(self, pos, rel_pos, methylation, type):
        self.pos = pos
        self.rel_pos = rel_pos
        self.methylation = methylation
        self.type = type

def main():
    parser = argparse.ArgumentParser(description='Calculate mean methylation of both flanking and insertion regions of all reads')
    parser.add_argument('-p', '--position_file', type=str, required=True, help='input file with each CpG in all reads')
    parser.add_argument('-g', '--groupby_file', type=str, required=True, help='input file of all reads with read type annotation')
    parser.add_argument('-o', '--output', type=str, required=True, help='output file')
    parser.add_argument('-l', '--len', type=int, default=500, help='length of flanking regions to average methylation')
    parser.add_argument('-a', '--aggregation', choices=['mean', 'count'], default= 'mean', help='how to aggregate the methylation of flanking regions. can be either mean or count')
    parser.add_argument('-t', '--threads', type=int, default=1, help='multi-threading')
    args = parser.parse_args()
    position_file = os.path.abspath(args.position_file)
    groupby_file = os.path.abspath(args.groupby_file)
    if not os.path.exists(position_file):
        raise ValueError("--position file does not exist!")
    if not os.path.exists(groupby_file):
        raise ValueError("--groupby file does not exist!")
    start_time = time.time()

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
            if locus_name not in locus_dict:
                locus_dict[locus_name] = Locus(chr, start, end)
                locus_dict[locus_name].add_read(read_item)
            else:
                locus_dict[locus_name].add_read(read_item)

    """
    position file:
    chr1	10861	10862	10082	m64136_200710_174522/11207845/ccs	-796	0.96	+	b
    chr1	10861	10862	10207	m64136_200710_174522/11207845/ccs	-670	0.98	+	b
    chr1	10861	10862	10468	m64136_200710_174522/11207845/ccs	-418	1.00	+	b
    """
    with open(position_file, 'r') as f:  # read the region file
        for line in f.readlines():
            chr, start, end, pos, read, rel_pos, methylation, strand, type = line.split()
            cpg_item = CpG(pos, rel_pos, methylation, type)
            locus = chr + ":" + start + "-" + end
            if locus in locus_dict and read in locus_dict[locus].reads:
                read = locus_dict[locus].get_read(read)
                read.add_cpg(cpg_item)

    with open(args.output, "w") as out:
        for locus in locus_dict:
            locus_dict[locus].aggregate_methylation(args.len, args.aggregation)
            # chr, chr.start, chr,end, chr.pos, read, read.pos, methylation, strand
            locus_object = locus_dict[locus]
            if hasattr(locus_object, "i_methylation"):
                out.write("{:s}\t{:d}\t{:d}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\n".format(
                            locus_object.chr, locus_object.start, locus_object.end, locus_object.empty_b_methylation, locus_object.empty_a_methylation, locus_object.empty_i_methylation, locus_object.insertion_b_methylation, locus_object.insertion_a_methylation, locus_object.i_methylation))
    end_time = time.time()
    print("--- %s hours ---" % ((end_time - start_time)/3600))

if __name__ == '__main__':
    main()