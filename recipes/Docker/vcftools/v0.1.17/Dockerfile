FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="${PATH}:/opt/vcftools/"
WORKDIR /opt/vcftools/

RUN apt-get update && apt-get install -y wget build-essential pkg-config zlib1g-dev perl

RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.17/vcftools-0.1.17.tar.gz
RUN tar -xzf vcftools-0.1.17.tar.gz
RUN rm -rf vcftools-0.1.17.tar.gz
RUN cd vcftools-0.1.17/ && ./configure && make && make install


CMD ["vcftools"]

