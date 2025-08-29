FROM ubuntu:latest
ENV PATH="${PATH}:/opt/prank/bin"
WORKDIR /opt/prank/
RUN apt-get update -y && apt-get install wget -y
RUN wget http://wasabiapp.org/download/prank/prank.linux64.170427.tgz
RUN tar -xvzf prank.linux64.170427.tgz
RUN mv prank/bin/ /opt/prank/bin/
RUN rm -rf prank
RUN rm prank.linux64.170427.tgz
RUN chmod -R +x /opt/prank/bin
CMD ["prank"]
