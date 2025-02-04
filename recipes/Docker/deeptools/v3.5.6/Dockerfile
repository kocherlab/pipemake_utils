FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/deeptools/bin:$PATH
RUN conda create -n deeptools -c conda-forge -c bioconda deeptools=3.5.6
RUN echo "source activate deeptools" > ~/.bashrc
