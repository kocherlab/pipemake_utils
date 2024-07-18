FROM ubuntu:latest
MAINTAINER aewebb80
ENV DEBIAN_FRONTEND noninteractive
ENV PATH="${PATH}:/opt/selscan/"
WORKDIR /opt/selscan/
RUN apt-get update -y && apt-get install wget tar -y && wget https://github.com/szpiech/selscan/archive/refs/tags/v2.0.0.tar.gz && tar -xzf v2.0.0.tar.gz && mv selscan-2.0.0/bin/linux/selscan selscan && mv selscan-2.0.0/bin/linux/norm norm && chmod +x selscan && chmod +x norm
CMD ["selscan"]