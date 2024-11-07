# Scripts to process BAM files with methylation modification tags (ML and MM tags).

## bam.Mmtag.regions.oo.py
usage: bam.Mmtag.regions.oo.py [-h] -b BAM [-t THREADS] -r REGION -f FORM [-l LEN_OFFSET] [-d DEPTH_FILTER]
                               [-w WINDOW] -o OUT

calculate CpG methylation average of inserted regions in the bam file

options:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     input bam file with Mm and Ml tags
  -t THREADS, --threads THREADS
                        multi-threading
  -r REGION, --region REGION
                        a vcf/bed file of genomic regions that will be used to summarize the methylation
  -f FORM, --form FORM  The format of the genomic region file, only vcf or bed-like files are supported.
  -l LEN_OFFSET, --len_offset LEN_OFFSET
                        length and coordinate offset from vcf to the bam file. Only insertions with coordinate
                        offset < len and SVlength difference < len will be considered
  -d DEPTH_FILTER, --depth_filter DEPTH_FILTER
                        regions with read depth higher than the threshold will be ignored
  -w WINDOW, --window WINDOW
                        flanking window on both ends of identified insertions of the vcf file. CpG methylation
                        within this window on both ends will be summarized
  -o OUT, --out OUT     output bed like txt file storing the methylation data in the defined regions



## read.methylation.mean.py 
usage: read.methylation.mean.py [-h] [-p POSITION_FILE] [-g GROUPBY_FILE] [-i INTEGRATED_FILE] -o OUTPUT
                                [-s START] [-l LEN] [-w WINDOW_NUMBER] [-c | --locus | --no-locus] [-t THREADS]

Calculate the mean methylation of both flanking and insertion regions of all reads

options:
  -h, --help            show this help message and exit
  -p POSITION_FILE, --position_file POSITION_FILE
                        input file with each CpG in all reads
  -g GROUPBY_FILE, --groupby_file GROUPBY_FILE
                        input file of all reads with read-type annotation
  -i INTEGRATED_FILE, --integrated_file INTEGRATED_FILE
                        The single input file of all reads if it contains all the info
  -o OUTPUT, --output OUTPUT
                        output file
  -s START, --start START
                        position outside of insertion to start averaging methylations
  -l LEN, --len LEN     length of flanking regions to average methylation
  -w WINDOW_NUMBER, --window_number WINDOW_NUMBER
                        number of sliding windows within the flanking region length. Only compatible while using
                        integrated file option
  -c, --locus, --no-locus
                        if true, process locus by locus instead of chr by chr (default: False)
  -t THREADS, --threads THREADS
                        multi-threading


## bam.extractInsertion.py
usage: bam.extractInsertion.py [-h] -b BAM -r REGION -o OUT [-s SAMPLE] [-e EXTEND]

Extract insertion and soft clipped sections intersecting coordinates from a bam file

options:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     input bam file
  -r REGION, --region REGION
                        a bed file of genomic regions that will be used to intersect with the bam file
  -o OUT, --out OUT     output extracted seq in unmapped bam format
  -s SAMPLE, --sample SAMPLE
                        sample name used as a prefix for the read names
  -e EXTEND, --extend EXTEND
                        extend the region by this many bp before extracting the sequence


## bam.creation.py
usage: bam.creation.py [-h] -f FASTA -m METHYL -b BED -o OUT [-c COLUMN]
                       [-n {coordinates,name,both}] [-p PREFIX]

create a bam file from a bed file. Similar to bedtools getfasta, but return a
unmapped bam file with methylation percentage in auxiliary tags.

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        input fasta file
  -m METHYL, --methyl METHYL
                        input methylation percentage file. For example, bed
                        file from pb-CpG-tools or bedmethyl file from modkit
                        pileup with traditional preset.
  -b BED, --bed BED     input bed file used to intersect and extract regions
  -o OUT, --out OUT     output unmapped bam with values from the bed file
                        in auxiliary tags
  -c COLUMN, --column COLUMN
                        column from the methyl file used for methylation
                        percentage that will be attached to the ML tag.
                        Default is 4 (4 for bed from pb-CpG-tools, 11 for
                        bedmethyl from modkit).
  -n {coordinates,name,both}, --name {coordinates,name,both}
                        How to name reads from bed file. Coordinates
                        (chrom:start-end) or name (column 4 from the bed file)
                        or both. Must be one of: coordinates, name, or both.
                        The default is both.
  -p PREFIX, --prefix PREFIX
                        prefix for all read names
