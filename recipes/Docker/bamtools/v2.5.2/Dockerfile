FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/bamtools/bin:$PATH
RUN conda create -n bamtools -c conda-forge bioconda::bamtools=2.5.2
RUN echo "source activate bamtools" > ~/.bashrc
