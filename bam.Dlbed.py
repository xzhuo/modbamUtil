import argparse
import pysam
import time
# from mpire import WorkerPool # multiprocessing

# convert bam file with Ml/Mm tags to bed file with methylation information in format: chr, start, end, methylated position array, unmethylated position array.
# One line per read.

def process_read(read):
    line = []
    if not (read.is_supplementary or read.is_secondary or read.is_unmapped):
        chrom = read.reference_name
        query_name = read.query_name
        start = read.reference_start
        end = read.reference_end
        strand = '-' if read.is_reverse else '+'
        modbase_key = ('C', 1, 'm') if read.is_reverse else ('C', 0, 'm')
        try:
            modbase_list = read.modified_bases[modbase_key]
            modbase_methy_list = [j[0] for j in list(filter(lambda i: i[1]/255 >= 0.5, modbase_list))]
            modbase_unmet_list = [j[0] for j in list(filter(lambda i: i[1]/255 < 0.5, modbase_list))]
            modbase_methy_string = ','.join(map(str, modbase_methy_list))
            modbase_unmet_string = ','.join(map(str, modbase_unmet_list))
        except:
            modbase_methy_string = ''
            modbase_unmet_string = ''
        line = [chrom, start, end, modbase_methy_string, modbase_unmet_string, strand]
    return line

def main():
    parser = argparse.ArgumentParser(description='calculate CpG methylation average of inserted regions in the bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file with Mm and Ml tags')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='output file')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='number of threads')
    args = parser.parse_args()
    bam_file = args.bam
    output_file = args.output
    threads = args.threads
    bam = pysam.AlignmentFile(bam_file, threads = 4, check_sq=False)
    num_reads = bam.count()
    print('total reads: {}'.format(num_reads))
    start_time = time.time()
    output = []
    threads = 1
    if threads == 1:
        for read in bam.fetch():
            output.append(process_read(read))
    # bam file multiprocessing is not working yet...
    # https://github.com/pysam-developers/pysam/issues/950
    # else:
    #     with WorkerPool(n_jobs=threads) as pool:
    #         output = pool.imap(process_read, (read for read in bam.fetch()), iterable_len = num_reads, progress_bar=True)

    with open(output_file, "w") as out:
        for line in output:
            if line:
                out.write('\t'.join(map(str, line)) + '\n')
    end_time = time.time()
    print("It took: --- %s hours ---" % ((end_time - start_time)/3600))

if __name__ == '__main__':
    main()
