# Create a unaligned bam file from a fasta file, and attach data from an associated bed file in auxiliary tags.
import os
import argparse
import re
import pysam
from Bio import SeqIO

def process_bed(bed):
    '''
    example bed file from intersectBed:
    chr1	82713606	82716730	chr1_82256025_82256027	0	+	chr1	82714191	82714192	64.3	Total	8	5	3	62.51
    chr1	82713606	82716730	chr1_82256025_82256027	0	+	chr1	82714438	82714439	73.1	Total	8	6	2	75.01
    '''
    ml_dict = {}
    with open(bed, "r") as f:
        for line in f:
            line = line.strip().split()
            name = line[3]
            pos = int(line[7]) - int(line[1])
            score = float(line[9])
            ml = round(score*255/100)
            if name not in ml_dict:
                ml_dict[name] = {"pos": [], "ml": []}
            ml_dict[name]["pos"].append(pos)
            ml_dict[name]["ml"].append(ml)
    return ml_dict

def attach (fasta, ml_dict, out, sample):
    header = { 'HD': {'VN':'1.6','SO':'unsorted','GO':'query'}}
    with pysam.AlignmentFile(out, "w", header=header) as outf:
        for record in SeqIO.parse(fasta, "fasta"):
            a = pysam.AlignedSegment()
            a.query_name = sample + ":" + record.id if sample else record.id
            a.query_sequence = str(record.seq)
            c_list = []
            for match in re.finditer(r'C', a.query_sequence, flags=re.IGNORECASE):
                c_list.append(match.start())
            try: 
                pos_dict = ml_dict[record.id]
                mm_list = [c_list.index(pos_dict["pos"][i]) - c_list.index(pos_dict["pos"][i-1]) - 1 if i > 0 else c_list.index(pos_dict["pos"][i]) for i in range(len(pos_dict["pos"]))]
                mm_tag = 'C+m?,' + ','.join([str(i) for i in mm_list]) + ';'
                a.flag = 4
                a.set_tag('MM', mm_tag, value_type='Z')
                a.set_tag('ML', value = pos_dict["ml"])
                print("imported to the sam file:", record.id)
                outf.write(a)
            except KeyError:
                print("No CG site for ", record.id)

def main():
    parser = argparse.ArgumentParser(description='create a bam file from a fasta file and a bed file')
    parser.add_argument('-f', '--fasta', type=str, required=True,
                        help='input fasta file')
    parser.add_argument('-b', '--bed', type=str, required=True,
                        help='input bed file')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output unaligned bam with with values from the bed file in auxiliary tags')
    parser.add_argument('-s', '--sample', type=str, required=False,
                        help='sample name used as prefix for the read names')

    args = parser.parse_args()
    fasta_file = os.path.abspath(args.fasta)
    bed_file = os.path.abspath(args.bed)

    if not os.path.exists(fasta_file):
        raise ValueError("--fasta file does not exist!")
    if not os.path.exists(bed_file):
        raise ValueError("--bed file does not exist!")
    ml_dict = process_bed(bed_file)
    attach(fasta_file, ml_dict, args.out, args.sample)



if __name__ == '__main__':
    main()