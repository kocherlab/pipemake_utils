FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/mafft/bin:$PATH
RUN conda create -n mafft -c conda-forge -c bioconda mafft=7.525
RUN echo "source activate mafft" > ~/.bashrc
