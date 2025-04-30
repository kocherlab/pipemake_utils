FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/diamond/bin:$PATH
RUN mamba create -n diamond -c conda-forge bioconda::diamond=2.1.11
RUN echo "source activate diamond" > ~/.bashrc