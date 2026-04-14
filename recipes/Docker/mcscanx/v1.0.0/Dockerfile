FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/mcscanx/bin:$PATH
RUN mamba create -n mcscanx -c conda-forge bioconda::mcscanx
RUN echo "source activate mcscanx" > ~/.bashrc