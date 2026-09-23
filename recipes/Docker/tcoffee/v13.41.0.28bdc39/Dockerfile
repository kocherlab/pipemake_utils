FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV HOME=/root

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        t-coffee \
        ca-certificates \
        perl \
    && rm -rf /var/lib/apt/lists/*

CMD ["t_coffee"]
