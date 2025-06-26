FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/ncbi-genome-tools/bin:$PATH
RUN mamba create -n ncbi-genome-tools python=3.9
RUN echo "source activate ncbi-genome-tools" > ~/.bashrc
RUN wget https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/add_utrs_to_gff/add_utrs_to_gff.py
RUN mv add_utrs_to_gff.py /opt/conda/envs/ncbi-genome-tools/bin/add_utrs_to_gff.py
RUN chmod +x /opt/conda/envs/ncbi-genome-tools/bin/add_utrs_to_gff.py