FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/eggnog-mapper/bin:$PATH
RUN conda create -n eggnog-mapper -c conda-forge bioconda::eggnog-mapper=2.1.13
RUN sed -i 's/eggnogdb\.embl\.de/eggnog5.embl.de/g' /opt/conda/envs/eggnog-mapper/bin/download_eggnog_data.py
RUN echo "source activate eggnog-mapper" > ~/.bashrc
