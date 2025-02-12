FROM continuumio/miniconda3
ENV PATH /opt/conda/envs/sleap/bin:$PATH
RUN conda create -n sleap -c sleap -c nvidia -c conda-forge sleap::sleap=1.2.8
RUN conda clean -afy
RUN rm -rf /opt/conda/pkgs
RUN find /opt/conda/ -type f -name '*.a' -delete
RUN find /opt/conda/ -type f -name '*.pyc' -delete
RUN find /opt/conda/ -type f -name '*.js.map' -delete
RUN rm -rf /tmp/* /var/tmp/*
RUN echo "source activate sleap" > ~/.bashrc
ENV LD_LIBRARY_PATH=/opt/conda/envs/sleap/lib
