FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/bcftools/bin:$PATH
COPY bcftools.yml .
RUN conda env create -f bcftools.yml
RUN echo "source activate bcftools" > ~/.bashrc
