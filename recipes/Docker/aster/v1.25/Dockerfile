FROM condaforge/mambaforge
ENV PATH=/opt/conda/envs/aster/bin:$PATH
RUN mamba create -n aster -c conda-forge bioconda::aster=1.25
RUN echo "source activate aster" > ~/.bashrc
