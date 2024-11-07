FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/eggnog-mapper/bin:$PATH
RUN conda create -n eggnog-mapper -c conda-forge bioconda::eggnog-mapper=2.1.12
RUN echo "source activate eggnog-mapper" > ~/.bashrc
