FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/shapeit4/bin:$PATH
COPY shapeit4.yml .
RUN conda env create -f shapeit4.yml
RUN echo "source activate shapeit4" > ~/.bashrc
