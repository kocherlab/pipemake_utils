FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/samtools/bin:$PATH
COPY samtools.yml .
RUN conda env create -f samtools.yml
RUN echo "source activate samtools" > ~/.bashrc
