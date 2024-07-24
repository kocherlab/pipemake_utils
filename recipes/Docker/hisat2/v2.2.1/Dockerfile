FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/hisat2/bin:$PATH
RUN conda create -n hisat2 -c conda-forge -c bioconda hisat2=2.2.1
RUN echo "source activate hisat2" > ~/.bashrc
