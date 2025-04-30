FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/sortmerna/bin:$PATH
RUN conda create -n sortmerna -c conda-forge bioconda::sortmerna=4.3.6
RUN echo "source activate sortmerna" > ~/.bashrc
RUN mkdir /opt/DBs
RUN wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
RUN tar -xvf database.tar.gz --directory /opt/DBs
RUN chmod -R 755 /opt/DBs
RUN rm database.tar.gz