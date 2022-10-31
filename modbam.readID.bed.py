from collections import defaultdict
from modbampy import ModBam
import argparse
import time

parser = argparse.ArgumentParser(description='calculate CpG methylation average of inserted regi    ons in the bam file')
parser.add_argument('-b', '--bam', type=str, required=True,
                    help='input bam file with Mm and Ml tags')
parser.add_argument('-c', '--chrom', type=str, required=True,
                    help='input chromSize file')
parser.add_argument('-w', '--window', type=int, default=100000000,
                    help='processing window size')
parser.add_argument('-o', '--out', type=str, required=True,
                    help='output bed like txt file storing the methylation data in the defined regions')

args = parser.parse_args()

size_list = []
start_time = time.time()
with open(args.size, 'r') as f:  # read the region file
    for line in f.readlines():
        chrom_list = line.strip().split()
        chrom_dict = {"chrom": chrom_list[0], "start": 0, "end": int(chrom_list[1])}
        size_list.append(chrom_dict)

with ModBam(args.bam) as bam:
    output = defaultdict(dict)
    for chrom_size in size_list:
        chrom = chrom_size["chrom"]
        start = chrom_size["start"]
        end = chrom_size["end"]
        for read in bam.reads(chrom, start, end):
            # print(read.mod_sites)
            for pos_mod in read.mod_sites:
                """read_id,
                reference position,
                query (read) position,
                reference strand (+ or -),
                modification strand (0 or 1, as defined in the HTSlib tag specification. This is invariable 0),
                canonical base associated with modification,
                modified base,
                modified-base score (scaled to 0-255)."""
                pos = pos_mod[1]
                strand = pos_mod[3]
                if not output[chrom][pos].has_key("strand"):
                    output[chrom][pos]["strand"] = strand
                if pos_mod[1] > start and pos_mod[1] <= end:
                    if pos_mod[7]/255 > 0.5:
                        output[chrom][pos]["methylated"].append(pos_mod[0])
                    else:
                        output[chrom][pos]["unmethylated"].append(pos_mod[0])
    
    with open(args.out, "w") as out:
        for chrom in output:
            for pos in output[chrom]:
                out.write("{:s}\t{:i}\t{:s}\t{:s}\t{:s}\n".format(
                    chrom, pos, output[chrom][pos]["strand"], ",".join(output[chrom][pos]["methylated"]), ",".join(output[chrom][pos]["unmethylated"])))

end_time = time.time()
print("--- %s hours ---" % ((end_time - start_time)/3600))
