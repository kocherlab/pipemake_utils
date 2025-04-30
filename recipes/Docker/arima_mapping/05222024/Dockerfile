FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/arima_mapping/bin:$PATH
RUN conda create -n arima_mapping -c conda-forge conda-forge::perl=5.32.1 bioconda::samtools=1.21 conda-forge::unzip
RUN echo "source activate arima_mapping" > ~/.bashrc
RUN wget https://github.com/ArimaGenomics/mapping_pipeline/archive/refs/heads/master.zip
RUN unzip master.zip
RUN rm master.zip
RUN sed -i '1s|.*|#!/usr/bin/env perl|' mapping_pipeline-master/filter_five_end.pl
RUN sed -i '1s|.*|#!/usr/bin/env perl|' mapping_pipeline-master/two_read_bam_combiner.pl
RUN sed -i '1s|.*|#!/usr/bin/env perl|' mapping_pipeline-master/get_stats.pl
RUN chmod -R 755 mapping_pipeline-master
RUN mv mapping_pipeline-master/* /opt/conda/envs/arima_mapping/bin
