# Create a unaligned bam file from a fasta file, and attach data from an associated bed file in auxiliary tags.
import os
import argparse
import re
import tempfile
import pybedtools
import pysam
# from Bio import SeqIO


def standardize_tmp_bed(bed, bed4):
    '''
    standardize bed file to 4 columns
    '''
    records = []
    with open(bed, "r") as f:
        for line in f:
            record = line.strip().split()
            if len(record) < 3:
                raise ValueError("bed file must have at least 3 columns")
            elif len(record) == 3:
                record.append(".")
            records.append(record[:4])

    with open(bed4, "w") as f:
        for record in records:
            f.write("\t".join(record) + "\n")

def process_beds(methyl, bed, column, name):
    '''
    intersect methyl file with bed file first.
    example bed file from intersectBed:
    chr1	82713606	82716730	chr1_82256025_82256027	chr1	82714191	82714192	64.3	Total	8	5	3	62.51
    chr1	82713606	82716730	chr1_82256025_82256027	chr1	82714438	82714439	73.1	Total	8	6	2	75.01
    '''

    # intersect bed file with methyl file with -wa -wb option. 
    methyl_file = pybedtools.BedTool(methyl)
    bed_file = pybedtools.BedTool(bed)
    intersect_bed = bed_file.intersect(methyl_file, wa=True, wb=True)

    ml_dict = {}
    with open(intersect_bed.fn, "r") as f:
        for line in f:
            line = line.strip().split()
            coordinates = line[0] + ":" + line[1] + "-" + line[2]
            id = line[3] if name == 'name' else coordinates if name == 'coordinates' else coordinates + "_" + line[3]
            pos = int(line[5]) - int(line[1])
            score = float(line[column + 3])
            ml = round(score*255/100)
            if id not in ml_dict:
                ml_dict[id] = {"pos": [], "ml": [], "chr": line[0], "start": int(line[1]), "end": int(line[2])}
            ml_dict[id]["pos"].append(pos)
            ml_dict[id]["ml"].append(ml)
    return ml_dict

def attach (fasta, ml_dict, out, prefix):
    header = { 'HD': {'VN':'1.6','SO':'unsorted','GO':'query'}}
    with pysam.AlignmentFile(out, "w", header=header) as outf:
        fasta = pysam.FastaFile(fasta)
        for id in ml_dict:
            a = pysam.AlignedSegment()
            pos_dict = ml_dict[id]
            a.query_sequence = fasta.fetch(reference=pos_dict["chr"], start=pos_dict["start"], end=pos_dict["end"])
            a.query_name = prefix + ":" + id if prefix else id
            c_list = []
            for match in re.finditer(r'C', a.query_sequence, flags=re.IGNORECASE):
                c_list.append(match.start())
            try:
                ref_CpG_list = [x for x in pos_dict["pos"] if x in c_list]
                mm_list = [c_list.index(ref_CpG_list[i]) - c_list.index(ref_CpG_list[i-1]) - 1 if i > 0 else c_list.index(ref_CpG_list[i]) for i in range(len(ref_CpG_list))]
                mm_tag = 'C+m?,' + ','.join([str(i) for i in mm_list]) + ';'
                a.flag = 4
                a.set_tag('MM', mm_tag, value_type='Z')
                a.set_tag('ML', value = pos_dict["ml"])
                print("imported to the sam file:", id)
                outf.write(a)
            except KeyError:
                print("No CG site for ", id)

def main():
    parser = argparse.ArgumentParser(description='create a bam file from a bed file. Similar to bedtools getfasta, but return a unmapped bam file with methylation percentage in auxiliary tags.')
    parser.add_argument('-f', '--fasta', type=str, required=True,
                        help='input fasta file')
    parser.add_argument('-m', '--methyl', type=str, required=True,
                        help='input methylation percentage file. For example, bed file from pb-CpG-tools or bedmethyl file from modkit pileup with traditional preset.')
    parser.add_argument('-b', '--bed', type=str, required=True,
                        help='input bed file used to intersect and extract regions')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output unmapped bam with with values from the bed file in auxiliary tags')
    parser.add_argument('-c', '--column', type=int, default=11,
                        help='column from the methyl file used for methylation percentage that will be attached to the ML tag. Default is 4 (4 for bed from pb-CpG-tools, 11 for bedmethyl from modkit).')
    parser.add_argument('-n', '--name', type=str, default='both',choices=['coordinates', 'name', 'both'],
                        help='How to name reads from bed file. Coordinates (chrom:start-end) or name (column 4 from the bed file) or both. Must be one of: coordiates, name, or both. Default is both.')
    parser.add_argument('-p', '--prefix', type=str, required=False,
                        help='prefix for the all read names')
    
    args = parser.parse_args()
    fasta_file = os.path.abspath(args.fasta)
    methyl_file = os.path.abspath(args.methyl)
    bed_file = os.path.abspath(args.bed)

    if not os.path.exists(fasta_file):
        raise ValueError("--fasta file does not exist!")
    if not os.path.exists(methyl_file):
        raise ValueError("--methyl file does not exist!")
    if not os.path.exists(bed_file):
        raise ValueError("--bed file does not exist!")

    # tmp_bed4_file = bed_file + ".tmp"
    tmp_bed4_file = tempfile.NamedTemporaryFile(mode='w', delete=False)

    standardize_tmp_bed(bed_file, tmp_bed4_file.name)
    ml_dict = process_beds(methyl_file, tmp_bed4_file.name, args.column, args.name)
    attach(fasta_file, ml_dict, args.out, args.prefix)
    tmp_bed4_file.close()
    os.unlink(tmp_bed4_file.name)


if __name__ == '__main__':
    main()