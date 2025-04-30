FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/bwa-mem2/bin:$PATH
RUN mamba create -n bwa-mem2 -c conda-forge bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.21
RUN echo "source activate bwa-mem2" > ~/.bashrc
