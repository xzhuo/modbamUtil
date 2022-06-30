import os
import argparse
import pysam
import time

def methylation_calculation(bam_file, out_file, len_filter):
    """
    calculate the methylation average of the inserted regions in the bam file
    """

    bam = pysam.AlignmentFile(bam_file, threads = 8, check_sq=False)
    out_list = []
    for read in bam.fetch(until_eof=True):
        if read.is_supplementary or read.is_secondary or read.is_unmapped:
            continue
        query_name = read.query_name
        modbase_key = ('C', 1, 'm') if read.is_reverse else ('C', 0, 'm')
        strand = '-' if read.is_reverse else '+'
        try:
            modbase_list = read.modified_bases[modbase_key] # a list of tuples
            last_match = None
            for i in read.get_aligned_pairs():
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
    parser.add_argument('-l', '--len', type=int, default=50,
                        help='lenght filter of the bam file. Only reads with length >= len will be considered')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output bed like txt file storing the methylation average of the inserted regions')

    args = parser.parse_args()
    bam_file = os.path.abspath(args.bam)
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    start_time = time.time()
    methylation_calculation(bam_file, args.out, args.len)
    end_time = time.time()
    print("--- %s seconds ---" % (end_time - start_time))


if __name__ == '__main__':
    main()
