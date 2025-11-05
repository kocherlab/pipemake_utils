FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/sniffles2_plot/bin:$PATH
RUN mamba create -n sniffles2_plot unzip python=3.11
RUN echo "source activate sniffles2_plot" > ~/.bashrc
RUN wget https://github.com/aewebb80/Sniffles2_plot/archive/refs/heads/master.zip
RUN unzip master.zip
RUN cd Sniffles2_plot-master/ && pip install .
CMD ["sniffles2_plot"]
