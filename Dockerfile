FROM mcfonsecalab/python37_bio:latest

MAINTAINER Pedro Barbosa <pedro.barbosa@medicina.ulisboa.pt>

WORKDIR /tools
ADD requirements.* /tmp/
RUN apt-get update && apt-get install -y build-essential zlib1g-dev libbz2-dev
RUN conda install --yes --file /tmp/requirements.conda
RUN pip install -r /tmp/requirements.pip
RUN git clone https://github.com/PedroBarbosa/VETA.git
RUN echo 'alias veta="python /tools/VETA/src/veta.py"' >> ~/.bashrc
