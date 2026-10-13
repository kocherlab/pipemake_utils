FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="${PATH}:/opt/selscan/"
WORKDIR /opt/selscan/
RUN apt-get update -y && apt-get install wget -y && wget https://github.com/szpiech/selscan/releases/download/v3.0.1/selscan-v3.0.1-linux
RUN mv selscan-v3.0.1-linux selscan && chmod +x selscan
CMD ["selscan"]