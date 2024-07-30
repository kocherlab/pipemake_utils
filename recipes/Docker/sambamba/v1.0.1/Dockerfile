FROM ubuntu:latest
MAINTAINER aewebb80
ENV DEBIAN_FRONTEND noninteractive
ENV PATH="${PATH}:/opt/sambamba/"
WORKDIR /opt/sambamba/
RUN apt-get update -y && apt-get install wget gzip -y
RUN wget https://github.com/biod/sambamba/releases/download/v1.0.1/sambamba-1.0.1-linux-amd64-static.gz
RUN gunzip sambamba-1.0.1-linux-amd64-static.gz
RUN mv sambamba-1.0.1-linux-amd64-static sambamba
RUN chmod +x sambamba
CMD ["sambamba"]
