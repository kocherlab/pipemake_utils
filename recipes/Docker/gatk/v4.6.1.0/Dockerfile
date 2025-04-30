FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/gatk4/bin:$PATH
RUN mamba create -n gatk4 -c conda-forge bioconda::gatk4=4.6.1.0 samtools=1.21
RUN echo "source activate gatk4" > ~/.bashrc