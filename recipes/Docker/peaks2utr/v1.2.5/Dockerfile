FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/peaks2utr/bin:$PATH
COPY peaks2utr.yml .
RUN conda env create -f peaks2utr.yml
RUN echo "source activate peaks2utr" > ~/.bashrc
