FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/genespace/bin:$PATH
RUN mamba create -n genespace conda-forge::r-base=4.4.3 bioconda::bioconductor-biocinstaller orthofinder=2.5.5 bioconda::bioconductor-rtracklayer conda-forge::r-devtools conda-forge::tar conda-forge::gzip conda-forge::liblzma-devel conda-forge::openjdk conda-forge::gxx=14.3.0 conda-forge::r-data.table
RUN echo "source activate genespace" > ~/.bashrc
RUN wget https://github.com/wyp1125/MCScanX/archive/refs/tags/v1.0.0.tar.gz
RUN tar -xzvf v1.0.0.tar.gz
RUN cd MCScanX-1.0.0 && make && mv MCScanX /opt/conda/envs/genespace/bin/ && mv MCScanX_h /opt/conda/envs/genespace/bin/ && mv duplicate_gene_classifier /opt/conda/envs/genespace/bin/ && mv downstream_analyses/* /opt/conda/envs/genespace/bin/
RUN R -e 'devtools::install_github("jtlovell/GENESPACE", upgrade="always")'

