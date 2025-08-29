FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/blobtools/bin:$PATH
RUN mamba create -n blobtools -c conda-forge bioconda::blobtools=1.1.1
RUN echo "source activate blobtools" > ~/.bashrc
RUN wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
RUN mkdir -p /opt/taxdump/
RUN tar -xzf new_taxdump.tar.gz -C /opt/taxdump/
RUN rm new_taxdump.tar.gz
RUN blobtools nodesdb --nodes /opt/taxdump/nodes.dmp --names /opt/taxdump/names.dmp
RUN rm -r /opt/taxdump/