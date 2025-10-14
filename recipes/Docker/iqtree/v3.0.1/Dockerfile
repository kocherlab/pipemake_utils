FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/iqtree/bin:$PATH
RUN mamba create -n iqtree -c conda-forge bioconda::iqtree=3.0.1
RUN echo "source activate iqtree" > ~/.bashrc
