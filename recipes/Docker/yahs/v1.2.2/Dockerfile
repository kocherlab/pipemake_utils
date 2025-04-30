FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/yahs/bin:$PATH
RUN mamba create -n yahs -c conda-forge bioconda::yahs=1.2.2
RUN echo "source activate yahs" > ~/.bashrc