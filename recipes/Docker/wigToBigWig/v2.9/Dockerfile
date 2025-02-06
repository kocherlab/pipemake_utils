FROM ubuntu:latest
MAINTAINER aewebb80
ENV DEBIAN_FRONTEND noninteractive
ENV PATH="${PATH}:/opt/wigtobigwig/"
WORKDIR /opt/wigtobigwig/
RUN apt-get update -y && apt-get install wget curl -y 
RUN wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
RUN chmod +x wigToBigWig
CMD ["wigToBigWig"]