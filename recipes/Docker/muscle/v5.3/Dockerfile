FROM ubuntu:latest
ENV PATH="${PATH}:/opt/muscle/"
WORKDIR /opt/muscle/
RUN apt-get update -y && apt-get install wget -y
RUN wget https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-linux-x86.v5.3
RUN mv muscle-linux-x86.v5.3 muscle
RUN chmod +x muscle
CMD ["muscle"]
