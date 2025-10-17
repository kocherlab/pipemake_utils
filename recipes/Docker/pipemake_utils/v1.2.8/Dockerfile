FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/pipemake_utils/bin:$PATH
RUN mamba create -n pipemake_utils -c conda-forge -c bioconda -c kocherlab conda-forge::unzip kocherlab::pipemake_utils=1.2.8
RUN echo "source activate pipemake_utils" > ~/.bashrc
 