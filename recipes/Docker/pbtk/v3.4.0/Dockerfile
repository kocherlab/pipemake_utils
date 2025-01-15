FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/pbtk/bin:$PATH
RUN conda create -n pbtk -c conda-forge -c bioconda pbtk=3.4.0
RUN echo "source activate pbtk" > ~/.bashrc
