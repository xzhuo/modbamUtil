FROM conda/miniconda3

RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -y git pysam mpire
RUN pip install modbampy scipy

RUN git clone https://github.com/xzhuo/modbamUtil.git
# RUN mkdir /modbamutil/
# ADD ./ /modbamutil