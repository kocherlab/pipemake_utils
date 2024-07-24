FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/bwa/bin:$PATH
RUN conda create -n bwa -c conda-forge -c bioconda bwa=0.7.18
RUN echo "source activate bwa" > ~/.bashrc
