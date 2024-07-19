FROM ubuntu:latest
MAINTAINER aewebb80
ENV DEBIAN_FRONTEND noninteractive
ENV PATH="${PATH}:/opt/star/"
WORKDIR /opt/star/
RUN apt-get update -y && apt-get install wget zip libgomp1 -y
RUN wget https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip
RUN unzip STAR_2.7.11b.zip
RUN mv STAR_2.7.11b/Linux_x86_64/STAR* . 
RUN chmod +x STAR
RUN chmod +x STARlong
CMD ["STAR"]
