FROM ubuntu:20.04

MAINTAINER Pedro Barbosa <pedro.barbosa@medicina.ulisboa.pt>

WORKDIR /tools
RUN apt-get update && apt-get install -y git bzip2 wget libssl-dev libcurl4-openssl-dev build-essential zlib1g-dev libbz2-dev
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p miniconda3

ENV PATH="$PATH:/tools/miniconda3/bin/"
ADD conda_environment.yml /tmp/
RUN conda env update -n base -f /tmp/conda_environment.yml
RUN git clone https://github.com/PedroBarbosa/VETA.git
RUN cd VETA && pip install .
