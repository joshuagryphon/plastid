FROM phusion/baseimage:focal-1.2.0
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
        bowtie \
        python \
        python-dev \
        python3.6 \
        python3.6-dev \
        python3.6-distutils \
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
RUN curl -o get-pip.py -sSL https://bootstrap.pypa.io/get-pip.py \
    && python3 get-pip.py "pip==22.0.4" \
    && pip install -r requirements-test.txt

# Download data required to run full test suite
RUN curl -L -o plastid/test/plastid_test_data.tar.bz2 \
        https://www.dropbox.com/s/np3wlfvp6gx8tb8/2022-05-04.plastid-test-data.tar.bz2?dl=0 \
    && cd plastid/test \
    && tar -jxvf plastid_test_data.tar.bz2 \
    && rm plastid_test_data.tar.bz2

# Configure test environments & verify build
RUN cat requirements.txt | sed -e "s/=.*//" >requirements-latest.txt \
    && tox -r --notest

# Set some useful variables
ENV HOME=/root
ENV MPLBACKEND=agg
ENV HOSTNAME=plastid

# Boot into bash terminal rather than run tests, because tests are slow
# and sometimes we only want to run a subset of them
CMD ["/bin/bash"]
