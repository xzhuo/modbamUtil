from curses import window
import os
import argparse
import pysam
import time

def intersect_methylation(bam_file, vcf_file, window, len_offset, out_file):
    """
    summarize the methylation of each CpG per read in the defined regions
    """
    bam = pysam.AlignmentFile(bam_file, threads = 8, check_sq=False)
    out_list = []
    with open(vcf_file, 'r') as f:  # read the bed file
        vcf_list = f.readlines()
        vcf_dict = dict(item.split("=") for item in vcf_list[7].split(";"))
        sv_len = int(vcf_dict['SVLEN'])
        for pileupcolumn in bam.pileup(vcf_list[0], vcf_list[1] - window, vcf_list[1] + window):
            ref_pos = pileupcolumn.reference_pos
            for pileupread in pileupcolumn.pileups:
                query_name = pileupread.alignment.query_name
                modbase_key = ('C', 1, 'm') if pileupread.alignment.is_reverse else ('C', 0, 'm')
                strand = '-' if pileupread.alignment.is_reverse else '+'
                if pileupread.is_del or pileupread.is_refskip:
                    continue
                if pileupread.indel == 0:
                    query_pos = pileupread.query_position
                    try:
                        modbase_perc = [j[1]/255 for j in list(filter(lambda i: i[0] == query_pos, pileupread.alignment.modified_bases[modbase_key]))][0]
                        dict = {'chr': vcf_list[0], 'ref_pos': ref_pos, 'query_name': query_name, 'query_pos': query_pos, 'modbase_perc': modbase_perc, 'strand': strand}
                    except:
                        pass

                elif pileupread.indel - sv_len < len_offset and ref_pos - vcf_list[1] < len_offset:
                    pass




                    last_match = None
                    for i in pileupread.alignment.get_aligned_pairs():
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
                                    query_seq = pileupread.alignment.query_sequence[query[0]:query[1]]
                                    modbase_pos_list = [j[0] - query[0] for j in list(filter(lambda i: i[0] >= query[0] and i[0] < query[1], modbase_list))]
                                    modbase_perc_list = [j[1]/255 for j in list(filter(lambda i: i[0] >= query[0] and i[0] < query[1], modbase_list))]
                                    modbase_pos_string = ','.join(["%d" % i for i in modbase_pos_list])
                                    modbase_string = ','.join(["%.2f" % i for i in mod

    bam.close()

    with open(out_file, "w") as out:
        for line in out_list:
            out.write("{:s}\t{:d}\t{:d}\t{:s}\t{:d}\t{:d}\t{:d}\t{:s}\t{:.4f}\t{:d}\t{:s}\t{:s}\t{:s}\n".format(
                line[0],
                line[1],
                line[2],
                line[3],
                line[4],
                line[5],
                line[6],
                line[7],
                line[8],
                line[9],
                line[10],
                line[11],
                line[12]))

def main():
    parser = argparse.ArgumentParser(description='calculate CpG methylation average of inserted regions in the bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file with Mm and Ml tags')
    parser.add_argument('-r', '--regions', type=str, required=True,
                        help='a vcf file of genomeic regions that will be used to summarize the methylation')
    parser.add_argument('-l', '--len', type=int, default=50,
                        help='length and coordinate offset from vcf to the bam file. Only insertions with coordinate offset < len and SVlength difference < len will be considered')
    parser.add_argument('-w', '--window', type=int, default=2000,
                        help='flanking window on both ends of identified insrtions of the vcf file. CpG methylation within this window on both ends will be summarized')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output bed like txt file storing the methylation data in the defined regions')

    args = parser.parse_args()
    bam_file = os.path.abspath(args.bam)
    vcf_file = os.path.abspath(args.regions)
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    if not os.path.exists(vcf_file):
        raise ValueError("--vcf file does not exist!")
    start_time = time.time()
    intersect_methylation(bam_file, vcf_file, args.window, args.len, args.out)
    end_time = time.time()
    print("--- %s seconds ---" % (end_time - start_time))


if __name__ == '__main__':
    main()
