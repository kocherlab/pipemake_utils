FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/aster/bin:$PATH"
WORKDIR /opt/aster/
RUN apt-get update -y && apt-get install wget unzip build-essential -y
RUN wget https://github.com/chaoszhang/ASTER/archive/refs/heads/master.zip
RUN unzip master.zip
RUN cd ASTER-master && make && mv bin/ /opt/aster/

