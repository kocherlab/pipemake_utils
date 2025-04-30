FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/ncbi-blast/bin:$PATH
RUN mamba create -n ncbi-blast -c conda-forge bioconda::blast=2.16.0
RUN echo "source activate ncbi-blast" > ~/.bashrc