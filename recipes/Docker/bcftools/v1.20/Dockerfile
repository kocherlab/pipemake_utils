FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/bcftools/bin:$PATH
RUN conda create -n bcftools -c conda-forge bioconda::bcftools=1.20
RUN echo "source activate bcftools" > ~/.bashrc
