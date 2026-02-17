FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Installer Python, outils de bioinfo et BLAST+
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    git \
    wget \
    curl \
    prodigal \
    hmmer \
    ncbi-blast+ \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Créer un alias python -> python3
RUN ln -s /usr/bin/python3 /usr/bin/python

# Installer DGR-package et Biopython
RUN pip3 install --no-cache-dir DGR-package biopython pandas

WORKDIR /data
CMD ["bash"]