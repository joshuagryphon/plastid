FROM ubuntu:18.04
LABEL maintainer "Joshua Griffin Dunn"

ARG DEBIAN_FRONTEND=noninteractive

ENV HOME=/root
ENV MPLBACKEND=agg

RUN apt-get update && apt-get install \
    --assume-yes \
    --verbose-versions \
    --allow-change-held-packages \
    -o Dpkg::Options::="--force-confdef" \
    build-essential \
    git \
    libssl-dev \
    sudo \
    wget \
    python-dev \
    python3.6-dev \
    python3.7-dev \
    python3.8-dev \
    python \
    python-pip \
    python3.6 \
    python3.7 \
    python3.8 \
    zlib1g-dev 


ENV PROJECT_HOME=/usr/src/plastid
WORKDIR $PROJECT_HOME
COPY requirements.txt ./
RUN pip install --upgrade pip && pip install -r requirements.txt

COPY . .
RUN python setup.py clean && pip install -e .
RUN python -c "from plastid import *"


CMD ["/bin/bash"]
