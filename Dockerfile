FROM conda/miniconda3

RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -y pysam mpire

RUN mkdir /modbamutil/
ADD ./ /modbamutil