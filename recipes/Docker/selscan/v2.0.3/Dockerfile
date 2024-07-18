FROM ubuntu:latest
MAINTAINER aewebb80
ENV DEBIAN_FRONTEND noninteractive
ENV PATH="${PATH}:/opt/selscan/"
WORKDIR /opt/selscan/
RUN apt-get update -y && apt-get install wget tar -y && wget https://github.com/szpiech/selscan/archive/refs/tags/v2.0.3.tar.gz && tar -xzf v2.0.3.tar.gz && mv selscan-2.0.3/bin/linux/selscan-2.0.3 selscan && mv selscan-2.0.3/bin/linux/norm norm && chmod +x selscan && chmod +x norm
CMD ["selscan"]