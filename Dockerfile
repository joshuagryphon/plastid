FROM phusion/baseimage:focal-1.1.0
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
        libbz2-dev \
        libcurl4 \
        libcurl4-openssl-dev \
        libfreetype6 \
        libfreetype6-dev \
        liblapack3 \
        liblapack-dev \
        liblzma-dev \
        libpng16-16 \
        libpng-dev \
        libssl-dev \
        pkg-config \
        software-properties-common \
        sudo \
        vim \
        zlib1g-dev 

# Install Python 3.6-3.9
RUN add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install \
        --assume-yes \
        --verbose-versions \
        --allow-change-held-packages \
        -o Dpkg::Options::="--force-confdef" \
        python \
        python-dev \
        python3.6 \
        python3.6-dev \
        python3.7 \
        python3.7-dev \
        python3.7-venv \
        python3.7-distutils \
        python3.9 \
        python3.9-dev \
        python3.9-venv \
        python3.9-distutils

# Copy source code into project
ENV PROJECT_HOME=/usr/src/plastid
WORKDIR $PROJECT_HOME
COPY . .

# Upgrade pip and install tox, which will handle installation of all
# dependencies inside virtual environments running various version of Python
#
# Note, we're installing pip and everything below under Python 3.6,
# while the system default version is 3.8.
#
# This may seem weird, but plastid is tested against package versions
# that no longer build under 3.8.
RUN curl -o get-pip.py -sSL https://bootstrap.pypa.io/get-pip.py \
    && python3.6 get-pip.py "pip==21.3.1" \
    && pip install tox

# We do this do enable documentation to build inside the container, without
# having to activate a tox env
RUN pip install -r requirements.txt \
    && python3.6 setup.py clean \
    && pip install -e .

# Verify successful build by importing
RUN python3.6 -c "from plastid import *"

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
