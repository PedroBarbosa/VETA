FROM ubuntu:18.04

MAINTAINER Pedro Barbosa <pedro.barbosa@medicina.ulisboa.pt>

WORKDIR /tools

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p miniconda3

ENV PATH="$PATH:/tools/miniconda3/bin/"
ADD requirements.* /tmp/

#RUN apt-get update && apt-get install -y build-essential zlib1g-dev libbz2-dev
RUN conda install --yes --file /tmp/requirements.conda
RUN pip install -r /tmp/requirements.pip
RUN git clone https://github.com/PedroBarbosa/VETA.git
#RUN echo 'alias veta="python /tools/VETA/src/veta.py"' >> ~/.bashrc
