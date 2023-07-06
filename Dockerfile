FROM ubuntu:20.04

MAINTAINER Pedro Barbosa <pedro.barbosa@medicina.ulisboa.pt>

WORKDIR /tools
RUN apt-get update && apt-get install -y git bzip2 wget libssl-dev libcurl4-openssl-dev build-essential zlib1g-dev libbz2-dev liblzma-dev
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py310_22.11.1-1-Linux-x86_64.sh && \
    bash Miniconda3-py310_22.11.1-1-Linux-x86_64.sh -b -p miniconda3

ENV PATH="$PATH:/tools/miniconda3/bin/"
ADD conda_environment.yml /tmp/
RUN conda env update -n base -f /tmp/conda_environment.yml
RUN git clone https://github.com/PedroBarbosa/VETA.git
RUN cd VETA && pip install .

RUN wget -c https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && \
    tar -xf bcftools-1.17.tar.bz2 && \
    tar -xf htslib-1.17.tar.bz2

WORKDIR /tools/htslib-1.17
RUN ./configure --prefix=/tools/miniconda3/envs/veta/lib && make && make install

WORKDIR /tools/bcftools-1.17
RUN ./configure --prefix=/tools/miniconda3/envs/veta/lib && make && make install
ENV PATH="/tools/miniconda3/envs/veta/lib/bin:${PATH}"

WORKDIR /tools
RUN rm -rf *sh htslib* bcftools*
