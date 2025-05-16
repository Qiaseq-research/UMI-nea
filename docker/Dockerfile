FROM python:3.8

RUN apt-get -y update && \
    apt-get -y install apt-utils sudo wget curl git vim bc less \
    gcc g++ cmake \
    libtool build-essential clang \
    zlib1g-dev libbz2-dev \
    libboost-all-dev libboost-iostreams-dev

RUN apt-get -y install samtools
RUN apt-get -y install bwa

COPY requirements.txt /tmp/

RUN pip install umi_tools
RUN pip install --requirement /tmp/requirements.txt
COPY network.py /usr/local/lib/python3.8/site-packages/umi_tools/

RUN mkdir -p /Download/
WORKDIR /Download/
RUN git clone -b develop --single-branch https://github.com/Qiaseq-research/UMI-nea.git
RUN cd UMI-nea/UMI-nea && make

RUN git clone https://github.com/fhlab/UMIC-seq.git
COPY UMIC-seq.py UMIC-seq/

CMD ["/bin/bash"]