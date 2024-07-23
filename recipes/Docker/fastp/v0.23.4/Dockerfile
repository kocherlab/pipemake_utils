FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/fastp/bin:$PATH
RUN conda create -n fastp -c conda-forge -c bioconda fastp=0.23.4
RUN echo "source activate fastp" > ~/.bashrc
