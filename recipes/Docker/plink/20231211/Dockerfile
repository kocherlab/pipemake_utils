FROM ubuntu:latest
MAINTAINER aewebb80
ENV DEBIAN_FRONTEND noninteractive
ENV PATH="${PATH}:/opt/plink/"
WORKDIR /opt/plink/
RUN apt-get update -y && apt-get install wget zip -y && wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip && unzip plink_linux_x86_64_20231211.zip && chmod +x plink && chmod +x prettify
CMD ["plink"]
