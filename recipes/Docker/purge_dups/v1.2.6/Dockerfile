FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/purge_dups/bin:$PATH
RUN conda create -n purge_dups -c conda-forge -c bioconda purge_dups=1.2.6 minimap2=2.28 busco=5.8.2 unzip
ENV LD_LIBRARY_PATH /opt/conda/envs/purge_dups/lib
RUN sed -i '/^export/d' /opt/conda/envs/purge_dups/bin/run_busco2
RUN sed -i '/^#export/d' /opt/conda/envs/purge_dups/bin/run_busco2
RUN sed -i '/^run_BUSCO/d' /opt/conda/envs/purge_dups/bin/run_busco2
RUN echo "busco -i \$fasta -o \$output_prefix -l \$lineage -m genome -c \$cpus" >> /opt/conda/envs/purge_dups/bin/run_busco2
RUN wget https://github.com/dfguan/runner/archive/refs/heads/master.zip
RUN unzip master.zip
RUN rm master.zip
RUN cd runner-master && pip install . && cd ..
RUN rm -r runner-master
RUN echo "source activate purge_dups" > ~/.bashrc