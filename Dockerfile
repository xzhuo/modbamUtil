# FROM conda/miniconda3
FROM condaforge/miniforge3

RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && mamba install -y wget git samtools bedtools minimap2 pysam mpire scipy pbjasmine pbmm2 pbtk 
    # pybedtools 
RUN pip install modbedtools biopython

# https://stackoverflow.com/questions/36996046/how-to-prevent-dockerfile-caching-git-clone
# --build-arg CACHEBUST=$(date +%s)
ARG CACHEBUST=1
RUN git clone https://github.com/xzhuo/modbamUtil.git
# RUN mkdir /modbamutil/
# ADD ./ /modbamutil

RUN wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v3.0.0/pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu.tar.gz
RUN tar -xzf pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu.tar.gz
RUN rm pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu.tar.gz

RUN wget https://github.com/nanoporetech/modkit/releases/download/v0.6.1/modkit_v0.6.1_u16_x86_64.tar.gz
RUN tar -xzf modkit_v0.6.1_u16_x86_64.tar.gz
RUN rm modkit_v0.6.1_u16_x86_64.tar.gz


ENV PATH="$PATH:/dist_modkit_v0.6.1_8fa79e3:/pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu/bin"
