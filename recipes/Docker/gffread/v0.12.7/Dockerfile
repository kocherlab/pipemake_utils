FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/gffread/bin:$PATH
RUN mamba create -n gffread -c conda-forge bioconda::gffread=0.12.7
RUN echo "source activate gffread" > ~/.bashrc