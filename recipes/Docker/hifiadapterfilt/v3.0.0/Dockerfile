FROM continuumio/miniconda3
ENV PATH=/opt/conda/envs/hifiadapterfilt/bin:$PATH
RUN conda create -n hifiadapterfilt -c conda-forge bioconda::hifiadapterfilt=3.0.0
RUN echo "source activate hifiadapterfilt" > ~/.bashrc