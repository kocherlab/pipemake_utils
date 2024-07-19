FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/macs3/bin:$PATH
COPY macs3.yml .
RUN conda env create -f macs3.yml
RUN echo "source activate macs3" > ~/.bashrc
