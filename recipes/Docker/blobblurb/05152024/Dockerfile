FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/blobblurb/bin:$PATH
RUN mamba create -n blobblurb unzip python=3.11
RUN echo "source activate blobblurb" > ~/.bashrc
RUN pip install pandas numpy==1.26.4 scipy
RUN wget https://github.com/aewebb80/blobblurb/archive/refs/heads/master.zip
RUN unzip master.zip
RUN rm master.zip
RUN sed -i '1s|.*|#!/usr/bin/env python|' blobblurb-master/blobblurb.py
RUN mv blobblurb-master/blobblurb.py /opt/conda/envs/blobblurb/bin/blobblurb
RUN rm -rf blobblurb-master/
RUN chmod -R 755 /opt/conda/envs/blobblurb/bin/blobblurb