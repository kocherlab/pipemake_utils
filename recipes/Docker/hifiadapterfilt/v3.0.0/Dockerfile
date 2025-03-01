FROM continuumio/miniconda3
ENV PATH=/opt/conda/envs/hifiadapterfilt/bin:$PATH
RUN conda create -n hifiadapterfilt -c conda-forge bioconda::hifiadapterfilt=3.0.0
RUN echo "source activate hifiadapterfilt" > ~/.bashrc
RUN sed -i "s/DBpath=.*/DBpath=\/opt\/conda\/envs\/hifiadapterfilt\/bin\/DB/g" /opt/conda/envs/hifiadapterfilt/bin/hifiadapterfilt.sh
RUN sed -i "s/sed 's\/\\\.fq\/\/'; done)/sed 's\/\\\.fq\/\/' | xargs -- basename; done)/g" /opt/conda/envs/hifiadapterfilt/bin/hifiadapterfilt.sh
RUN sed -i "s/read_path=\\$.*/read_path=\\$\(dirname \\$\{prefix\}\)/g" /opt/conda/envs/hifiadapterfilt/bin/hifiadapterfilt.sh
RUN sed -i "s/readssprefix/readsprefix/g" /opt/conda/envs/hifiadapterfilt/bin/hifiadapterfiltFCS.sh
RUN sed -i '$d' /opt/conda/envs/hifiadapterfilt/bin/hifiadapterfiltFCS.sh