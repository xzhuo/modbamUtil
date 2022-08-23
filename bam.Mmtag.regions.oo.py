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
            self.start = int(line_list[1])
            self.end = int(line_list[1])
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

def intersect_methylation(bam_file, interval, window, len_offset):
    """
    summarize the methylation of each CpG per read in the defined regions
    """

    bam = pysam.AlignmentFile(bam_file, threads = 4, check_sq=False)
    out_list = []
    flanking_window = (interval.start - window, interval.end + window)
    for read in bam.fetch(interval.chr, flanking_window[0], flanking_window[1]):
        if read.is_supplementary or read.is_secondary or read.is_unmapped:
            continue
        aligned_pairs = read.get_aligned_pairs()
        query_flanking_window = 
        query_name = read.query_name
        modbase_key = ('C', 1, 'm') if read.is_reverse else ('C', 0, 'm')
        strand = '-' if read.is_reverse else '+'
        try:
            modbase_list = read.modified_bases[modbase_key] # a list of tuples
            modbase_pos_list = [j[0] for j in list(filter(lambda i: i[0] >= query[0] and i[0] < query[1], modbase_list))]
            modbase_perc_list = [j[1]/255 for j in list(filter(lambda i: i[0] >= query[0] and i[0] < query[1], modbase_list))]
            last_match = None
            for i in read.get_aligned_pairs():
                if i[1] < interval.start - window or i[1] > interval.end + window:
                    continue
                if i[0] is None:
                    # deletion
                    try:
                        deletion_length +=1
                    except UnboundLocalError:
                        deletion_length = 1
                elif i[1] is None:
                    # insertion
                    try:
                        insertion_length +=1
                    except UnboundLocalError:
                        insertion_length = 1
                else:
                    # match
                    if last_match is not None:
                        if insertion_length > len_filter:
                            ref = (i[1],i[1])
                            query = (i[0] - insertion_length, i[0])
                            query_seq = read.query_sequence[query[0]:query[1]]
                            modbase_pos_list = [j[0] - query[0] for j in list(filter(lambda i: i[0] >= query[0] and i[0] < query[1], modbase_list))]
                            modbase_perc_list = [j[1]/255 for j in list(filter(lambda i: i[0] >= query[0] and i[0] < query[1], modbase_list))]
                            modbase_pos_string = ','.join(["%d" % i for i in modbase_pos_list])
                            modbase_string = ','.join(["%.2f" % i for i in modbase_perc_list])
                            modbase_count = len(modbase_perc_list)
                            if len(modbase_perc_list) > 0:
                                modbase_perc = sum(modbase_perc_list)/len(modbase_perc_list)
                            else:
                                modbase_perc = -1
                            out_list.append([read.reference_name, ref[0], ref[1], query_name, query[0], query[1], query[1]-query[0], strand, modbase_perc, modbase_count, modbase_string, modbase_pos_string, query_seq])

                        elif deletion_length > len_filter:
                            ref = (i[1] - deletion_length,i[1])
                            query = (i[0], i[0])

                    insertion_length = 0
                    deletion_length = 0
                    last_match = i
        except:
            pass

    for pileupcolumn in bam.pileup(interval.chr, interval.sv_pos - window, interval.sv_pos + window, truncate=True):
        ref_pos = pileupcolumn.reference_pos
        for pileupread in pileupcolumn.pileups:
            pairs = pileupread.alignment.get_aligned_pairs()
            try:
                query_sv_pos = [i[0] for i in pairs if i[1] == interval.sv_pos][0]
            except IndexError:
                query_sv_pos = 0.5
            query_name = pileupread.alignment.query_name
            modbase_key = ('C', 1, 'm') if pileupread.alignment.is_reverse else ('C', 0, 'm')
            strand = '-' if pileupread.alignment.is_reverse else '+'
            if pileupread.is_del or pileupread.is_refskip:
                continue
            query_pos = pileupread.query_position
            ins_len = pileupread.indel
            if ins_len == 0:
                try:
                    modbase_perc = [j[1]/255 for j in list(filter(lambda i: i[0] == query_pos, pileupread.alignment.modified_bases[modbase_key]))][0]
                    methyl_dict = {'chr': interval.chr, 'ref_pos': ref_pos, 'query_name': query_name, 'query_pos': query_pos, 'rel_query_pos': query_pos - query_sv_pos, 'modbase_perc': modbase_perc, 'strand': strand, 'id': interval.indel_id, 'sv_len': interval.sv_len, 'ins_len': ins_len, 'type': 'flanking', 'mei': interval.mei, 'mei_strand': interval.mei_strand}
                    out_list.append(methyl_dict)
                except:
                    pass

            elif ins_len >= 50 and ref_pos - interval.sv_pos < len_offset:
                query = (pileupread.query_position, pileupread.query_position + ins_len)
                # query_seq = pileupread.alignment.query_sequence[query[0]:query[1]]
                try:
                    for j in list(filter(lambda i: i[0] >= query[0] and i[0] < query[1], pileupread.alignment.modified_bases[modbase_key])):
                        modbase_perc = j[1]/255
                        methyl_dict = {'chr': interval.chr, 'ref_pos': ref_pos, 'query_name': query_name, 'query_pos': j[0], 'rel_query_pos': j[0] - query_sv_pos, 'modbase_perc': modbase_perc, 'strand': strand, 'id': interval.indel_id, 'sv_len': interval.sv_len, 'ins_len': ins_len, 'type': 'insertion', 'mei': interval.mei, 'mei_strand': interval.mei_strand}
                        out_list.append(methyl_dict)
                except:
                    pass

    bam.close()

    return out_list



def main():
    parser = argparse.ArgumentParser(description='calculate CpG methylation average of inserted regions in the bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file with Mm and Ml tags')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='multi-threading')
    parser.add_argument('-r', '--region', type=str, required=True,
                        help='a vcf/bed file of genomeic regions that will be used to summarize the methylation')
    parser.add_argument('-f', '--file', type=str, required=True,
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
            line_item = Interval(line.strip(), "vcf")
            interval_array.append(line_item)

    with WorkerPool(n_jobs=args.threads) as pool:
        outputs = pool.imap(intersect_methylation, zip(repeat(bam_file), interval_array, repeat(args.window), repeat(args.len)), iterable_len=len(interval_array), progress_bar=True)

    with open(args.out, "w") as out:
        for out_list in outputs:
            for me in out_list:
                out.write("{:s}\t{:d}\t{:s}\t{:d}\t{:d}\t{:d}\t{:d}\t{:s}\t{:.2f}\t{:s}\t{:s}\t{:s}\t{:s}\n".format(
                    me['chr'], me['ref_pos'], me['query_name'], me['query_pos'], me['rel_query_pos'], me['sv_len'], me['ins_len'], me['strand'], me['modbase_perc'], me['id'], me['type'], me['mei'], me['mei_strand']))
    end_time = time.time()
    print("--- %s hours ---" % ((end_time - start_time)/3600))


if __name__ == '__main__':
    main()
