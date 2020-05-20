FROM ubuntu:18.04
LABEL maintainer "Joshua Griffin Dunn"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install \
        --assume-yes \
        --verbose-versions \
        --allow-change-held-packages \
        -o Dpkg::Options::="--force-confdef" \
        build-essential \
        curl \
        git \
        gfortran \
        libatlas3-base \
        libatlas-base-dev \
        libfreetype6 \
        libfreetype6-dev \
        liblapack3 \
        liblapack-dev \
        libpng16-16 \
        libpng-dev \
        libssl-dev \
        pkg-config \
        python3 \
        python3-pip \
        python3.6-dev \
        python3.7-dev \
        python3.8-dev \
        python3.6 \
        python3.7 \
        python3.8 \
        software-properties-common \
        sudo \
        vim \
        zlib1g-dev 

# Install Python 3.5 as a legacy version
RUN add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install \
        --assume-yes \
        --verbose-versions \
        --allow-change-held-packages \
        -o Dpkg::Options::="--force-confdef" \
        python3.5 \
        python3.5-dev

# Upgrade pip3 and install tox, which will handle installation of all
# dependencies inside virtual environments running various version of Python
RUN pip3 install --upgrade pip \
    && pip3 install tox


# Copy source code into project
ENV PROJECT_HOME=/usr/src/plastid
WORKDIR $PROJECT_HOME
COPY . .

# Install default in Python 3
RUN python3 setup.py clean && pip3 install -e .

# Verify successful build by importing
RUN python3 -c "from plastid import *"

# Download data required to run full test suite
RUN curl -L -o plastid/test/plastid_test_data.tar.bz2 \
        https://www.dropbox.com/s/h17go7tnas4hpby/plastid_test_data.tar.bz2?dl=0 \
    && cd plastid/test \
    && tar -jxvf plastid_test_data.tar.bz2 \
    && rm plastid_test_data.tar.bz2

# Force build of C-extensions for all test versions
# RUN tox -r --notest


# Set some useful variables
ENV HOME=/root
ENV MPLBACKEND=agg
ENV HOSTNAME=plastid

# Boot into bash terminal rather than run tests, because tests are slow
# and sometimes we only want to run a subset of them
CMD ["/bin/bash"]
