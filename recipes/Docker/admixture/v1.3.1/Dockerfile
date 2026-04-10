FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/admixture:$PATH"
WORKDIR /opt/admixture/
RUN apt-get update -y && apt-get install wget -y 
RUN wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.1.tar.gz
RUN tar -xzf admixture_linux-1.3.1.tar.gz
RUN mv admixture_linux-1.3.1/admixture .
RUN chmod +x admixture
CMD ["admixture"]

