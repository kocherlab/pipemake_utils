FROM condaforge/mambaforge
ENV PATH="/opt/conda/envs/clipkit/bin:$PATH"
RUN mamba create -n clipkit bioconda::clipkit
RUN echo "source activate clipkit" > ~/.bashrc
