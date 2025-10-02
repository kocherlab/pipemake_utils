FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/miniprot/bin:$PATH
RUN mamba create -n miniprot -c conda-forge bioconda::miniprot=0.18
RUN echo "source activate miniprot" > ~/.bashrc