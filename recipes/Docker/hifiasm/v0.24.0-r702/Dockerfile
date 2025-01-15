FROM ubuntu:latest
MAINTAINER aewebb80
ENV DEBIAN_FRONTEND noninteractive
ENV PATH="${PATH}:/opt/hifiasm/"
WORKDIR /opt/hifiasm/
RUN apt-get update -y && apt-get install wget zip build-essential libz-dev -y
RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/0.24.0.zip
RUN unzip 0.24.0.zip 
RUN cd hifiasm-0.24.0 && make && cd ..
RUN mv hifiasm-0.24.0/hifiasm .
RUN chmod +x hifiasm
CMD ["hifiasm"]
