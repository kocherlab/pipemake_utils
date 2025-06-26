FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/pigz/bin:$PATH
RUN mamba create -n pigz conda-forge::pigz=2.8
RUN echo "source activate pigz" > ~/.bashrc