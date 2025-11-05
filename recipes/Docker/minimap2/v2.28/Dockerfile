FROM ubuntu:latest
MAINTAINER aewebb80
ENV DEBIAN_FRONTEND noninteractive
ENV PATH="${PATH}:/opt/minimap2/bin"
WORKDIR /opt/minimap2/bin
RUN apt-get update -y && apt-get install wget bzip2 build-essential libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev -y
RUN wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2
RUN tar -xjf minimap2-2.28_x64-linux.tar.bz2
RUN mv minimap2-2.28_x64-linux/minimap2 .
RUN chmod +x minimap2
RUN mv minimap2-2.28_x64-linux/paftools.js .
RUN chmod +x paftools.js
RUN mv minimap2-2.28_x64-linux/k8 .
RUN chmod +x k8
RUN rm -rf minimap2-2.28_x64-linux/
RUN rm minimap2-2.28_x64-linux.tar.bz2
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
RUN tar -xjf samtools-1.21.tar.bz2
RUN cd samtools-1.21/ && ./configure --prefix=/opt/minimap2/ && make && make install
RUN rm -rf samtools-1.21/
RUN rm samtools-1.21.tar.bz2
CMD ["minimap2"]