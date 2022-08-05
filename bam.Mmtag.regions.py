import os
import argparse
import pysam
from itertools import repeat
from mpire import WorkerPool
import time

def intersect_methylation(bam_file, vcf_line, window, len_offset):
    """
    summarize the methylation of each CpG per read in the defined regions
    """

    bam = pysam.AlignmentFile(bam_file, threads = 4, check_sq=False)
    vcf_list = vcf_line.split('\t')
    vcf_dict = dict(item.split("=") for item in vcf_list[7].split(";"))
    sv_len = int(vcf_dict['SVLEN'])
    indel_id = vcf_list[2]
    sv_pos = int(vcf_list[1])
    out_list = []
    mei = vcf_list[13] if len(vcf_list) == 14 else '.'
    for pileupcolumn in bam.pileup(vcf_list[0], sv_pos - window, sv_pos + window, truncate=True):
        ref_pos = pileupcolumn.reference_pos
        for pileupread in pileupcolumn.pileups:
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
                    methyl_dict = {'chr': vcf_list[0], 'ref_pos': ref_pos, 'query_name': query_name, 'query_pos': query_pos, 'modbase_perc': modbase_perc, 'strand': strand, 'id': indel_id, 'sv_len': sv_len, 'ins_len': ins_len, 'type': 'flanking', 'mei': mei}
                    out_list.append(methyl_dict)
                except:
                    pass

            elif ins_len >= 50 and ref_pos - sv_pos < len_offset:
                query = (pileupread.query_position, pileupread.query_position + ins_len)
                # query_seq = pileupread.alignment.query_sequence[query[0]:query[1]]
                try:
                    for j in list(filter(lambda i: i[0] >= query[0] and i[0] < query[1], pileupread.alignment.modified_bases[modbase_key])):
                        modbase_perc = j[1]/255
                        methyl_dict = {'chr': vcf_list[0], 'ref_pos': ref_pos, 'query_name': query_name, 'query_pos': j[0], 'modbase_perc': modbase_perc, 'strand': strand, 'id': indel_id, 'sv_len': sv_len, 'ins_len': ins_len, 'type': 'insertion', 'mei': mei}
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
    parser.add_argument('-v', '--vcf', type=str, required=True,
                        help='a vcf file of genomeic regions that will be used to summarize the methylation')
    parser.add_argument('-l', '--len', type=int, default=50,
                        help='length and coordinate offset from vcf to the bam file. Only insertions with coordinate offset < len and SVlength difference < len will be considered')
    parser.add_argument('-w', '--window', type=int, default=2000,
                        help='flanking window on both ends of identified insrtions of the vcf file. CpG methylation within this window on both ends will be summarized')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output bed like txt file storing the methylation data in the defined regions')

    args = parser.parse_args()
    bam_file = os.path.abspath(args.bam)
    vcf_file = os.path.abspath(args.vcf)
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    if not os.path.exists(vcf_file):
        raise ValueError("--vcf file does not exist!")
    start_time = time.time()
    vcf_array = []
    outputs = []
    with open(vcf_file, 'r') as f:  # read the bed file
        for line in f.readlines():
            vcf_array.append(line.strip())

    with WorkerPool(n_jobs=args.threads) as pool:
        outputs = pool.imap(intersect_methylation, zip(repeat(bam_file), vcf_array, repeat(args.window), repeat(args.len)), iterable_len=len(vcf_array), progress_bar=True)

    with open(args.out, "w") as out:
        for out_list in outputs:
            for me in out_list:
                out.write("{:s}\t{:d}\t{:s}\t{:d}\t{:d}\t{:d}\t{:s}\t{:.2f}\t{:s}\t{:s}\t{:s}\n".format(
                    me['chr'], me['ref_pos'], me['query_name'], me['query_pos'], me['sv_len'], me['ins_len'], me['strand'], me['modbase_perc'], me['id'], me['type'], me['mei']))
    end_time = time.time()
    print("--- %s seconds ---" % (end_time - start_time))


if __name__ == '__main__':
    main()
