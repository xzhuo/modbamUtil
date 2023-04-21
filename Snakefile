SAMPLES = ["bam"]

rule all:
    input:
        expand("{sample}.pdf", sample=SAMPLES),

rule pbmm2:
    input:
        bam = "{sample}.bam",
        ref = "/storage1/fs1/hprc/Active/xzhuo/ref/hg38.fa"
    output:
        bam = "{sample}.{ref}.bam",
        bai = "{sample}.{ref}.bam.bai"
    threads:
        8
    params:
        sort_threads = 2,
        alignment_threads = {threads} - {params.sort_threads}
    shell:
        "pbmm2 align {input.ref} {input.bam} {output.bam} --preset CCS --sort -j {params.alignment_threads} -J {params.sort_threads}"
    container:
        "quay.io/biocontainers/pbmm2:1.10.0--h9ee0642_0"

rule pb-CpG-tools:
    input:
        bam = "{sample}.{ref}.bam"
    output:
        model = "{sample}.{ref}.model.combined.bed",
        count = "{sample}.{ref}.count.combined.bed"
    threads:
        8
    params:
        model = "/pb-CpG-tools-v2.2.0-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite"
    run:
        shell("aligned_bam_to_cpg_scores --bam {input} --output-prefix {sample}.{ref}.model --model {params.model} --threads {threads}")
        shell("aligned_bam_to_cpg_scores --bam {input} --output-prefix {sample}.{ref}.count --pileup-mode count --threads {threads}")
    container:
        "docker.io/xiaoyuz/modbamutil:latest"

rule ggpairs_wgbs:
    input:
        model = "{sample}.{ref}.model.combined.bed",
        count = "{sample}.{ref}.count.combined.bed",
        wgbs = "HG002.wgbs.methylC.CpG.txt"
    output:
        bed = "HG002.wgbs.{sample}.methylC.CpG.txt",
        pdf = "{sample}.pdf"
    run:
        bedtools intersect -a {input.wgbs} -b {input.count} -loj |cut -f1-5,9,10 |bedtools intersect -a stdin -b {input.model} -loj |cut -f1-7,11,12> {output.bed}
        perl -lane 'print join("\t",$F[0],$F[1],$F[2],$F[3],$F[9],$F[11],$F[13],$F[15]) if $F[4]>=5 && $F[10] >=5 && $F[12]>=5 && $F[14]>=5 && $F[16]>=5' HG002.wgbs.hifi.megalodon.revio_count.revio_model.sequel_count.sequel_model.methylC.CpG.txt > HG002.wgbs.revio_count.revio_model.sequel_count.sequel_model.ggpairs.txt

        Rscript ggpairs.revio.r HG002.wgbs.revio_count.revio_model.sequel_count.sequel_model.ggpairs.txt