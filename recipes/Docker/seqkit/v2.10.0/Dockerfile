FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/seqkit/bin:$PATH
RUN mamba create -n seqkit -c conda-forge bioconda::seqkit=2.10.0
RUN echo "source activate seqkit" > ~/.bashrc