FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/tabix/bin:$PATH
RUN mamba create -n tabix -c conda-forge bioconda::tabix=1.11
RUN echo "source activate tabix" > ~/.bashrc
