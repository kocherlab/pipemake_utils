FROM ubuntu:latest
MAINTAINER aewebb80
ENV DEBIAN_FRONTEND noninteractive
ENV PATH="${PATH}:/opt/admixture/"
WORKDIR /opt/admixture/
RUN apt-get update -y && apt-get install wget zip -y 
RUN wget https://github.com/NovembreLab/admixture/archive/refs/heads/master.zip
RUN unzip master.zip 
RUN mv admixture-master/releases/admixture_linux-1.3.0/admixture .
RUN chmod +x admixture
CMD ["admixture"]
