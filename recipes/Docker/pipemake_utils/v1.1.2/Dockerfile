FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/pipemake_utils/bin:$PATH
RUN conda create -n pipemake_utils -c conda-forge -c bioconda -c kocherlab kocherlab::pipemake_utils=1.1.2
RUN echo "source activate pipemake_utils" > ~/.bashrc
