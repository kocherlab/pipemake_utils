FROM ubuntu:24.04
RUN apt-get update && apt-get install --no-install-recommends -y build-essential cmake wget unzip
RUN groupadd -r hyphyuser && useradd -r -g hyphyuser hyphyuser
WORKDIR /hyphy
RUN wget --no-check-certificate https://github.com/veg/hyphy/archive/refs/tags/2.5.97.tar.gz
RUN wget --no-check-certificate https://github.com/veg/hyphy-analyses/archive/refs/heads/master.zip
RUN tar -xzf 2.5.97.tar.gz
RUN rm 2.5.97.tar.gz
RUN unzip master.zip
RUN mv hyphy-analyses-master hyphy-analyses
RUN mv hyphy-2.5.97/cmake /hyphy/cmake
RUN mv hyphy-2.5.97/src /hyphy/src
RUN mv hyphy-2.5.97/contrib /hyphy/contrib
RUN mv hyphy-2.5.97/res /hyphy/res
RUN mv hyphy-2.5.97/CMakeLists.txt /hyphy
RUN mv hyphy-2.5.97/tests /hyphy/tests
RUN chown -R hyphyuser:hyphyuser .
RUN cmake . && make -j install
USER hyphyuser