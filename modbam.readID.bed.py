from collections import defaultdict
from modbampy import ModBam
import argparse

parser = argparse.ArgumentParser(description='calculate CpG methylation average of inserted regi    ons in the bam file')
parser.add_argument('-b', '--bam', type=str, required=True,
                    help='input bam file with Mm and Ml tags')
parser.add_argument('-c', '--chrom', type=str, required=True,
                    help='input chromSize file')
parser.add_argument('-w', '--window', type=int, default=100000000,
                    help='processing window size')

args = parser.parse_args()

args.size

chrom = "chr1"
start = 120680
end = 120685

with ModBam(args.bam) as bam:
    for read in bam.reads(chrom, start, end):
        # print(read.mod_sites)
        output = defaultdict(dict)
        for pos_mod in read.mod_sites:
            """read_id,
            reference position,
            query (read) position,
            reference strand (+ or -),
            modification strand (0 or 1, as defined in the HTSlib tag specification. This is invariable 0),
            canonical base associated with modification,
            modified base,
            modified-base score (scaled to 0-255)."""
            pos = pos_mod[1] if pos_mod[3] == "+" else pos_mod[1] - 1

            if pos_mod[1] > start and pos_mod[1] <= end:
                if pos_mod[7]/255 > 0.5:
                    output[chrom][pos]["methylated"].append(pos_mod[0])
                else:
                    output[chrom][pos]["unmethylated"].append(pos_mod[0])
            