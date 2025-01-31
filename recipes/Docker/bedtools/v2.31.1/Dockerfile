FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/bedtools/bin:$PATH
RUN conda create -n bedtools -c conda-forge bioconda::bedtools=2.31.1
RUN echo "source activate bedtools" > ~/.bashrc
