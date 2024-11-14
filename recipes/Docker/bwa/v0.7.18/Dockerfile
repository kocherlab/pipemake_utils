FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/bwa/bin:$PATH
RUN conda create -n bwa -c conda-forge bioconda::bwa=0.7.18 bioconda::samtools=1.21
RUN echo "source activate bwa" > ~/.bashrc
