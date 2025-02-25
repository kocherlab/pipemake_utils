FROM continuumio/miniconda3
ENV PATH=/opt/conda/envs/flair/bin:$PATH
RUN conda create -n flair -c conda-forge -c bioconda bioconda::flair=2.0.0
RUN echo "source activate flair" > ~/.bashrc
