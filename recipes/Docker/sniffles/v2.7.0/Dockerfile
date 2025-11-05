FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/sniffles/bin:$PATH
RUN mamba create -n sniffles -c conda-forge bioconda::sniffles=2.7.0
RUN echo "source activate sniffles" > ~/.bashrc