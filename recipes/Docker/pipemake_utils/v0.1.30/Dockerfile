FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/pipemake_utils/bin:$PATH
COPY pipemake_utils.yml .
RUN conda env create -f pipemake_utils.yml
RUN echo "source activate pipemake_utils" > ~/.bashrc
