FROM conda/miniconda3

RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -y git pysam mpire scipy
RUN pip install modbampy

# https://stackoverflow.com/questions/36996046/how-to-prevent-dockerfile-caching-git-clone
ARG CACHEBUST=1
RUN git clone https://github.com/xzhuo/modbamUtil.git
# RUN mkdir /modbamutil/
# ADD ./ /modbamutil