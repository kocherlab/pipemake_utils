FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/sortmerna/bin:$PATH
WORKDIR /opt/
RUN wget https://github.com/sortmerna/sortmerna/releases/download/v4.3.7/sortmerna-4.3.7-conda-linux-64.tar.bz2
RUN conda create -n sortmerna /opt/sortmerna-4.3.7-conda-linux-64.tar.bz2
RUN echo "source activate sortmerna" > ~/.bashrc
RUN mkdir /opt/DBs
RUN wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
RUN tar -xvf database.tar.gz --directory /opt/DBs
RUN chmod -R 755 /opt/DBs
RUN rm database.tar.gz