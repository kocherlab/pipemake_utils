FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/sniffles2/bin:$PATH
RUN mamba create -n sniffles2 -c conda-forge bioconda::sniffles=2.6.3
RUN echo "source activate sniffles2" > ~/.bashrc
