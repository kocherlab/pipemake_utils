FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/subread/bin:$PATH
COPY subread.yml .
RUN conda env create -f subread.yml
RUN echo "source activate subread" > ~/.bashrc
