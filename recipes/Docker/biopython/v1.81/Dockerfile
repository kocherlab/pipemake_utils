FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/biopython/bin:$PATH
COPY biopython.yml .
RUN conda env create -f biopython.yml
RUN echo "source activate biopython" > ~/.bashrc
