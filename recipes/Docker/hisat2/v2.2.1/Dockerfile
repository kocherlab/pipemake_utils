FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/hisat2/bin:$PATH
RUN mamba create -n hisat2 -c conda-forge -c bioconda hisat2=2.2.1 bioconda::samtools=1.21
RUN echo "source activate hisat2" > ~/.bashrc