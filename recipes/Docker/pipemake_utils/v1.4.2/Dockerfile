FROM condaforge/mambaforge
ENV PATH=/opt/conda/envs/pipemake_utils/bin:$PATH
RUN mamba create -n pipemake_utils -c bioconda conda-forge::unzip kocherlab::pipemake_utils=1.4.2
RUN echo "source activate pipemake_utils" > ~/.bashrc
 