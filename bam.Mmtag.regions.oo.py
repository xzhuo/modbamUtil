import sys
import os
import argparse
import pysam
from itertools import repeat
from mpire import WorkerPool
import time

class Interval:
    def __init__(self, line, type):
        self.line = line
        self.type = type
        self.dissect_line()

    def dissect_line(self):
        line_list = self.line.split('\t')
        if self.type == "vcf":
            # example vcf line:
            # chr1    136603  svim_asm.INS.2  A       AGGCTGACCTCTGTCCGCGTGGGAGGGGCCGGTGTGAGGCAAGGGCTCG       .       PASS    SVTYPE=INS;END=136603;SVLEN=48;READS=HG00741#2#JAHALX010000154.1        GT      0/1     .       -1      -1      .       .
            vcf_dict = dict(item.split("=") for item in line_list[7].split(";"))
            self.chr = line_list[0]
            self.pos = int(line_list[1])
            self.start = self.pos
            self.end = self.pos
            self.indel_id = line_list[2]
            self.sv_len = int(vcf_dict['SVLEN'])
            self.mei = line_list[13] if len(line_list) == 15 else '.'
            self.mei_strand = line_list[14] if len(line_list) == 15 else '.'
        elif self.type == "bed":
            # example bed file line:
            # chr1    90257   90258   4       HG00621,HG01952,HG01978,HG03516 1       0       1       1       1
            self.chr = line_list[0]
            self.start = int(line_list[1])
            self.end = int(line_list[2])
        else:
            sys.exit("unknow interval type!")
    
    def attach_modbase_list(self, window, bam):
        flanking_window = (self.start - window if self.start > window else 0, self.end + window)
        self.modbase_list = []
        for read in bam.fetch(self.chr, flanking_window[0], flanking_window[1]):
            if read.is_supplementary or read.is_secondary or read.is_unmapped:
                continue
            all_ref_pos = read.get_reference_positions()
            read_ref_start = min(all_ref_pos, key=lambda x:abs(x-flanking_window[0]))
            read_ref_end = min(all_ref_pos, key=lambda x:abs(x-flanking_window[1]))
            ref_pos = min(all_ref_pos, key=lambda x:abs(x-self.start))
            get_pos = convert_pos(read)
            query_flanking_start = get_pos['find_query'][read_ref_start]
            query_flanking_end = get_pos['find_query'][read_ref_end]
            query_start = get_pos['find_query'][ref_pos]
            query_name = read.query_name
            modbase_key = ('C', 1, 'm') if read.is_reverse else ('C', 0, 'm')
            strand = '-' if read.is_reverse else '+'
            try:
                modbase_list = read.modified_bases[modbase_key]
                modbase_query_list = [j[0] for j in list(filter(lambda i: i[0] >= query_flanking_start and i[0] < query_flanking_end, modbase_list))]
                modbase_ref_list = [get_pos['find_ref'][i] if i in get_pos['find_ref'] else -1 for i in modbase_query_list]
                modbase_rel_pos_list = [i - query_start for i in modbase_query_list]
                modbase_perc_list = [j[1]/255 for j in list(filter(lambda i: i[0] >= query_flanking_start and i[0] < query_flanking_end, modbase_list))]
                modbase_pos_perc = list(zip(modbase_rel_pos_list, modbase_ref_list, modbase_perc_list))
            except:
                modbase_pos_perc = []
            self.modbase_list.append({'read': query_name, 'strand': strand, 'modbase_pos_perc': modbase_pos_perc})


def convert_pos(read):
    result_dict = {'find_ref': {}, 'find_query': {}}
    for i in read.get_aligned_pairs(matches_only=True):
        result_dict['find_ref'][i[0]] = i[1]
        result_dict['find_query'][i[1]] = i[0]
    return result_dict


def intersect_methylation(bam_file, interval, window):
    """
    summarize the methylation of each CpG per read in the defined regions
    """

    out_list = []
    bam = pysam.AlignmentFile(bam_file, threads = 4, check_sq=False)
    interval.attach_modbase_list(window, bam)
    for i in interval.modbase_list:
        for j in i['modbase_pos_perc']:
            # chr, chr.start, chr,end, chr.pos, read, read.pos, methylation, strand
            out_list.append([interval.chr, interval.start, interval.end, j[1], i['read'], j[0], j[2], i['strand']])
    return out_list


def main():
    parser = argparse.ArgumentParser(description='calculate CpG methylation average of inserted regions in the bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file with Mm and Ml tags')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='multi-threading')
    parser.add_argument('-r', '--region', type=str, required=True,
                        help='a vcf/bed file of genomeic regions that will be used to summarize the methylation')
    parser.add_argument('-f', '--form', type=str, required=True,
                        help='The file format of the genomic region file, only vcf or bed like files are supported.')
    parser.add_argument('-l', '--len', type=int, default=50,
                        help='length and coordinate offset from vcf to the bam file. Only insertions with coordinate offset < len and SVlength difference < len will be considered')
    parser.add_argument('-w', '--window', type=int, default=2000,
                        help='flanking window on both ends of identified insrtions of the vcf file. CpG methylation within this window on both ends will be summarized')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output bed like txt file storing the methylation data in the defined regions')

    args = parser.parse_args()
    bam_file = os.path.abspath(args.bam)
    region_file = os.path.abspath(args.region)
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    if not os.path.exists(region_file):
        raise ValueError("--region file does not exist!")
    start_time = time.time()
    interval_array = []
    outputs = []
    with open(region_file, 'r') as f:  # read the region file
        for line in f.readlines():
            line_item = Interval(line.strip(), args.form)
            interval_array.append(line_item)

    if args.threads == 1:
        for i in interval_array:
            outputs.append(intersect_methylation(bam_file, i, args.window))
    else:
        with WorkerPool(n_jobs=args.threads) as pool:
            outputs = pool.imap(intersect_methylation, zip(repeat(bam_file), interval_array, repeat(args.window)), iterable_len=len(interval_array), progress_bar=True)

    with open(args.out, "w") as out:
        for out_list in outputs:
            last_status = ''
            for i in out_list:
                # perl -lape '$status = $F[5]<0?"b":"a";$status = "i" if $status eq "a" && $last_status eq "b" && $F[3] == -1;$status = "i" if $last_status eq "i" && $F[3] == -1;$last_status=$status;$_.="\t$status"'
                status
                if i[5] < 0:
                    status = 'b'
                elif i[3] == -1 and (last_status == 'b' or last_status == 'i'):
                    status = 'i'
                else:
                    status = 'a'
                last_status = status
                # chr, chr.start, chr,end, chr.pos, read, read.pos, methylation, strand
                out.write("{:s}\t{:d}\t{:d}\t{:d}\t{:s}\t{:d}\t{:.2f}\t{:s}\t{:s}\n".format(
                        i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],status))
    end_time = time.time()
    print("--- %s hours ---" % ((end_time - start_time)/3600))


if __name__ == '__main__':
    main()
