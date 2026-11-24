FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/kmc/bin:$PATH"
WORKDIR /opt/kmc

RUN apt-get update -y \
 && apt-get install -y --no-install-recommends wget ca-certificates \
 && wget -q https://github.com/refresh-bio/KMC/releases/download/v3.2.4/KMC3.2.4.linux.x64.tar.gz \
 && tar -xzf KMC3.2.4.linux.x64.tar.gz \
 && rm KMC3.2.4.linux.x64.tar.gz \
 && apt-get purge -y wget \
 && apt-get autoremove -y \
 && rm -rf /var/lib/apt/lists/*

CMD ["kmc"]
