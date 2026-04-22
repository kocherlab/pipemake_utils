FROM condaforge/mambaforge
ENV PATH=/opt/conda/envs/macse/bin:$PATH
RUN mamba create -n macse bioconda::macse=2.07
RUN echo "source activate macse" > ~/.bashrc