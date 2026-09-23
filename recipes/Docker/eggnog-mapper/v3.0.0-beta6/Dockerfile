# syntax=docker/dockerfile:1

# Recipe adapted from https://github.com/eggnogdb/eggnog-mapper/blob/main/apptainer/eggnog-mapper.def
FROM python:3.12-slim-bookworm

ARG EMAPPER_VERSION=3.0.0-beta6
ARG DIAMOND_VERSION=2.1.24
ARG MMSEQS_VERSION=latest
ARG HMMER_VERSION=3.4
ARG PRODIGAL_VERSION=2.6.3

LABEL org.opencontainers.image.title="eggNOG-mapper" \
      org.opencontainers.image.version="${EMAPPER_VERSION}" \
      org.opencontainers.image.description="Fast functional annotation using eggNOG orthology assignments" \
      org.opencontainers.image.source="https://github.com/eggnogdb/eggnog-mapper" \
      maintainer="j.huerta@upm.es" \
      diamond_version="${DIAMOND_VERSION}" \
      mmseqs2_version="${MMSEQS_VERSION}" \
      hmmer_version="${HMMER_VERSION}" \
      prodigal_version="${PRODIGAL_VERSION}"

ENV PATH=/opt/diamond:/opt/mmseqs/bin:/opt/hmmer/bin:/opt/prodigal:$PATH \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    HOME=/root \
    TMPDIR=/tmp \
    DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq

RUN apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        gcc g++ make \
        zlib1g-dev \
        libgomp1

RUN rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/diamond && \
    curl -sL \
        "https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VERSION}/diamond-linux64.tar.gz" \
        | tar -xz -C /opt/diamond

RUN chmod +x /opt/diamond/diamond

RUN mkdir -p /opt/mmseqs

RUN if [ "${MMSEQS_VERSION}" = "latest" ]; then \
        curl -sL https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz | tar -xz -C /opt; \
    else \
        curl -sL \
            "https://github.com/soedinglab/MMseqs2/releases/download/${MMSEQS_VERSION}/mmseqs-linux-avx2.tar.gz" \
            | tar -xz -C /opt; \
    fi

RUN mkdir -p /tmp/hmmer-build && \
    curl -sL "http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz" \
        | tar -xz -C /tmp/hmmer-build --strip-components=1

RUN cd /tmp/hmmer-build && ./configure --prefix=/opt/hmmer > /dev/null 2>&1

RUN cd /tmp/hmmer-build && make -j"$(nproc)" > /dev/null 2>&1

RUN cd /tmp/hmmer-build && make install > /dev/null 2>&1

RUN rm -rf /tmp/hmmer-build

RUN mkdir -p /opt/prodigal

RUN curl -sL \
        "https://github.com/hyattpd/Prodigal/releases/download/v${PRODIGAL_VERSION}/prodigal.linux" \
        -o /opt/prodigal/prodigal

RUN chmod +x /opt/prodigal/prodigal

RUN python3.12 -m pip install --break-system-packages --no-cache-dir "cython>=3.0"

RUN mkdir -p /opt/emapper-src

RUN curl -sL \
        "https://github.com/eggnogdb/eggnog-mapper/archive/refs/tags/v${EMAPPER_VERSION}.tar.gz" \
        | tar -xz -C /opt/emapper-src --strip-components=1

RUN python3.12 -m pip install --break-system-packages --no-cache-dir /opt/emapper-src

RUN rm -rf /opt/emapper-src

RUN emapper.py --version

RUN diamond --version

RUN mmseqs version

RUN hmmsearch -h 2>&1 | head -2

RUN prodigal -v 2>&1 | head -2

ENTRYPOINT ["emapper.py"]
