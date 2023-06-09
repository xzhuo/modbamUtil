import os
import sys
import argparse
from mpire import WorkerPool

class Read:
    def __init__(self, line):
        self.chr, self.start, self.end, self.read, self.CGins, self.readtype = line.split()
        self.cpg = []

    def add_cpg(self, pos, rel_pos, methylation, type):
        cpg = {"pos": pos,"rel_pos": rel_pos, "methylation": methylation, "type": type}
        self.cpg.append(cpg)

    def mean_methylation(self,type,length):
        methylation_list = []
        if type == "b":
            methylation_list = [i["methylation"] for i in self.cpg if i["rel_pos"]<0 and i["rel_pos"]>length * -1 and pos > 0]
        elif type == "a":
            methylation_list = [i["methylation"] for i in self.cpg if i["rel_pos"]>0 and i["rel_pos"]<length and pos > 0]
        elif type == "i":
            methylation_list = [i["methylation"] for i in self.cpg if i["type"] == "i"]
        else:
            sys.exit("unknow region type, should be b, a, or i.")

        mean_methylation = sum(methylation_list)/len(methylation_list)
        return mean_methylation


def main():
    parser = argparse.ArgumentParser(description='Calculate mean methylation of both flanking and insertion regions of all reads')
    parser.add_argument('-p', '--position_file', type=str, required=True, help='input file with each CpG in all reads')
    parser.add_argument('-g', '--groupby_file', type=str, required=True, help='input file of all reads with read type annotation')
    parser.add_argument('-o', '--output', type=str, required=True, help='output file')
    parser.add_argument('-l', '--len', type=int, default=50, help='length of flanking regions to average methylation')
    parser.add_argument('-t', '--threads', type=int, default=1, help='multi-threading')
    args = parser.parse_args()
    position_file = os.path.abspath(args.position_file)
    groupby_file = os.path.abspath(args.groupby_file)
    if not os.path.exists(position_file):
        raise ValueError("--position file does not exist!")
    if not os.path.exists(groupby_file):
        raise ValueError("--groupby file does not exist!")
    output = os.path.abspath(args.output)

    read_dict = {}
    with open(groupby_file, 'r') as f:  # read the region file
        for line in f.readlines():
            line_item = Read(line)
            read_dict[line_item.read] = line_item

if __name__ == '__main__':
    main()