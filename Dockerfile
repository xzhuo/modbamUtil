FROM conda/miniconda3

RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -y wget git pysam mpire scipy
RUN pip install modbampy

# https://stackoverflow.com/questions/36996046/how-to-prevent-dockerfile-caching-git-clone
# --build-arg CACHEBUST=$(date +%s)
ARG CACHEBUST=1
RUN git clone https://github.com/xzhuo/modbamUtil.git
# RUN mkdir /modbamutil/
# ADD ./ /modbamutil

RUN wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.2.0/pb-CpG-tools-v2.2.0-x86_64-unknown-linux-gnu.tar.gz
RUN tar -xzf pb-CpG-tools-v2.2.0-x86_64-unknown-linux-gnu.tar.gz
ENV PATH="$PATH:/pb-CpG-tools-v2.2.0-x86_64-unknown-linux-gnu/bin"
