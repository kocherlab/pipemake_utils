FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/juicer_tools/bin:$PATH
RUN mamba create -n juicer_tools -c conda-forge conda-forge:openjdk=23.0.2
RUN echo "source activate juicer_tools" > ~/.bashrc
RUN wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.19.02.jar
RUN mv juicer_tools_1.19.02.jar /opt/juicer_tools.jar