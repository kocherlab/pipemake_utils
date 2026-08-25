FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/bioconda/bin:$PATH
RUN mamba create -n bioconda bioconda::orthofinder=3.1.5 python=3.11
RUN echo "source activate bioconda" > ~/.bashrc
