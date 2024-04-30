FROM ubuntu:latest
MAINTAINER aewebb80
ENV DEBIAN_FRONTEND noninteractive
ENV PATH="${PATH}:/opt/plink/"
WORKDIR /opt/plink/
RUN apt-get update -y && apt-get install wget zip -y && wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20240418.zip && unzip plink2_linux_x86_64_20240418.zip && chmod +x plink2
CMD ["plink2"]
