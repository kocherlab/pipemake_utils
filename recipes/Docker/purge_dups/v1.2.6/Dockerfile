FROM continuumio/miniconda3
RUN apt-get update -y && apt-get install wget zip bzip2 libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev build-essential libz-dev -y
ENV PATH /opt/conda/envs/purge_dups/bin:$PATH
RUN conda install -c conda-forge mamba
RUN mamba create -n purge_dups -c conda-forge -c bioconda minimap2=2.28 busco python=3.7
ENV LD_LIBRARY_PATH /opt/conda/envs/purge_dups/lib
RUN wget https://github.com/dfguan/purge_dups/archive/refs/heads/master.zip
RUN unzip master.zip
RUN rm master.zip
RUN cd purge_dups-master/src && make && cd ../..
RUN cp -rf purge_dups-master/bin/* /opt/conda/envs/purge_dups/bin
RUN cp -rf purge_dups-master/scripts/* /opt/conda/envs/purge_dups/bin
RUN sed -i '/^export/d' /opt/conda/envs/purge_dups/bin/run_busco2
RUN sed -i '/^run_BUSCO/d' /opt/conda/envs/purge_dups/bin/run_busco2
RUN echo "busco -i \$fasta -o \$output_prefix -l \$lineage -m genome -c \$cpus" >> /opt/conda/envs/purge_dups/bin/run_busco2
RUN rm -rf purge_dups-master
RUN wget https://github.com/dfguan/runner/archive/refs/heads/master.zip
RUN unzip master.zip
RUN cd runner-master && pip install . && cd ..
RUN echo "source activate purge_dups" > ~/.bashrc