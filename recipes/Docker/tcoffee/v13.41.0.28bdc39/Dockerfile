FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV HOME=/root

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        t-coffee \
        ca-certificates \
        perl \
        wget \
        build-essential \
        cmake \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/TimoLassmann/kalign/archive/refs/tags/v3.5.1.tar.gz && \
    tar -xzf v3.5.1.tar.gz && \
    cd kalign-3.5.1 && \
    mkdir build && cd build && \
    cmake .. && \
    make && \
    make test && \
    make install

RUN wget https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-linux-x86.v5.3 && \
    mv muscle-linux-x86.v5.3 /usr/local/bin/muscle && \
    chmod +x /usr/local/bin/muscle
    
RUN wget https://gitlab.com/sysimm/mafft/-/archive/main/mafft-main.tar.gz && \
    tar -xzf mafft-main.tar.gz && \
    cd mafft-main/core && \
    make && \
    make install

CMD ["t_coffee"]
