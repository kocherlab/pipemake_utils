FROM ubuntu:latest
MAINTAINER aewebb80
ENV DEBIAN_FRONTEND noninteractive
ENV PATH="${PATH}:/opt/jellyfish/bin"
WORKDIR /opt/jellyfish/bin
RUN apt-get update -y && apt-get install wget build-essential -y 
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.1/jellyfish-2.3.1.tar.gz
RUN tar -xzf jellyfish-2.3.1.tar.gz
RUN cd jellyfish-2.3.1/ && ./configure --prefix=/opt/jellyfish/ && make -j 4 && make install
RUN chmod +x jellyfish
CMD ["jellyfish"]
