FROM ubuntu:latest
ENV PATH="${PATH}:/opt/prank/bin"
WORKDIR /opt/prank/
RUN apt-get update -y && apt-get install wget -y
RUN wget http://wasabiapp.org/download/prank/development/prank.linux64.170703.tgz
RUN tar -xvzf prank.linux64.170703.tgz
RUN mv  development/prank/bin/ /opt/prank/bin/
RUN rm -rf development/
RUN rm prank.linux64.170703.tgz
RUN chmod -R +x /opt/prank/bin
CMD ["prank"]
