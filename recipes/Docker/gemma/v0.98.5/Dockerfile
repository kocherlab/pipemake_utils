FROM ubuntu:latest
MAINTAINER aewebb80
ENV DEBIAN_FRONTEND noninteractive
ENV PATH="${PATH}:/opt/gemma/"
WORKDIR /opt/gemma/
RUN apt-get update -y && apt-get install wget gzip -y && wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz && gunzip gemma-0.98.5-linux-static-AMD64.gz && mv gemma-0.98.5-linux-static-AMD64 gemma && chmod +x gemma
CMD ["gemma"]
