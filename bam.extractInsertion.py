import os
import argparse
import pysam
from Bio.Seq import Seq

def find_insertion(read, soft_clip=True, len_threshold=50):
    """
    find the insertion query positions in the read
    """
    curr_pos = 0
    insertions = []
    for operation, length in read.cigartuples:
        if operation == 1 and length > len_threshold:
            insertions.append((curr_pos, curr_pos + length))
        if soft_clip and operation == 4 and length > len_threshold:
            insertions.append((curr_pos, curr_pos + length))
        if operation == 0 or operation == 1 or operation == 4 or operation == 7 or operation == 8:
            curr_pos += length
    return insertions

def extract_insertion(bam_file, region_file, sample, out, extend):
    with open(region_file, "r") as f:
        regions = [line.strip().split() for line in f]
    bam = pysam.AlignmentFile(bam_file, "rb")
    header = bam.header
    header["HD"]["VN"] = "1.6"
    header["HD"]["SO"] = "unsorted"
    header["HD"]["GO"] = "query"
    with pysam.AlignmentFile(out, "wb", header=header) as outf:
        for region in regions:
            region_id = str(sample) + ":" + region[0] + ":" + region[1] + "-" + region[2]
            print("Extracting reads from region: ", region_id)
            for read in bam.fetch(region[0], int(region[1]), int(region[2])):
                if read.is_supplementary or read.is_secondary or read.is_unmapped:
                    continue
                else:
                    if read.cigartuples:
                        aligned_pairs = read.get_aligned_pairs(matches_only=True)
                        query_range = [x[0] for x in aligned_pairs if x[1]>int(region[1])-extend and x[1]<int(region[2])+extend]
                        range_start = min(query_range)-1
                        range_end = max(query_range)+1
                        insertions = find_insertion(read)
                        if insertions:
                            for start, end in insertions:
                                if (start >= range_start and start <= range_end) or (end >= range_start and end <= range_end):
                                    a = pysam.AlignedSegment()
                                    a.query_name = read.query_name
                                    a.query_sequence = Seq(read.query_sequence[start:end]).reverse_complement() if read.is_reverse else read.query_sequence[start:end]
                                    a.query_qualities = read.get_forward_qualities()[start:end]
                                    a.flag = 4
                                    a.set_tag("HP", region_id)
                                    outf.write(a)
                                    print(f"Extracted subseq from {read.query_name} inserted in {region_id}")
                        else:
                            print("No insertion for read: ", read.query_name)
                            continue
                    else:
                        print("No cigar string for read: ", read.query_name)
                        continue

    bam.close()

def main():
    parser = argparse.ArgumentParser(description='Extract insertion and soft clipped sections intersecting coordinates from a bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file')
    parser.add_argument('-r', '--region', type=str, required=True,
                        help='a bed file of genomeic regions that will be used to intersect with the bam file')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output extracted seq in unmapped bam format')
    parser.add_argument('-s', '--sample', type=str, required=False,
                        help='sample name used as prefix for the read names')
    parser.add_argument('-e', '--extend', type=int, default=20,
                        help='extend the region by this many bp before extracting the sequence')

    args = parser.parse_args()
    bam_file = os.path.abspath(args.bam)
    region_file = os.path.abspath(args.region)
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    if not os.path.exists(region_file):
        raise ValueError("--region file does not exist!")
    extract_insertion(bam_file, region_file, args.sample, args.out, args.extend)

if __name__ == '__main__':
    main()