FROM ubuntu:latest

# Update package list and install necessary tools
RUN apt-get update && \
    apt-get install -y build-essential cmake git libboost-all-dev

RUN cd /home/ && git clone https://github.com/Qiaseq-research/UMI-nea.git

RUN cd /home/UMI-nea/UMI-nea && make

RUN mkdir /home/UMI-nea/data
