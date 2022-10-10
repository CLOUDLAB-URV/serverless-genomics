############################################################
# Dockerfile for variant caller - cloudbutton genomics use case
# Based on ubuntu 20.04
############################################################
ARG FUNCTION_DIR="/function"
ARG BINS_DIR="/function/bin"

# Set the base image to Ubuntu
FROM ubuntu:20.04

# set user
USER root

# add credentials for boto3/aws
RUN mkdir -p /root/.aws && touch /root/.aws/credentials
#COPY aws_credentials /root/.aws/credentials
ENV AWS_SHARED_CREDENTIALS_FILE=/root/.aws/credentials

# add lithops config
RUN mkdir -p /root/.lithops && touch /root/.lithops/config
#COPY lithops_config /root/.lithops/config
ENV LITHOPS_CONFIG_FILE=/root/.lithops/config

# Create function work directories
ARG FUNCTION_DIR
ARG BINS_DIR
RUN mkdir -p ${FUNCTION_DIR}
RUN mkdir -p ${BINS_DIR}

# Add timezone so that tzdata does not hang to ask for timezone
# https://grigorkh.medium.com/fix-tzdata-hangs-docker-image-build-cdb52cc3360d
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Add basic utilities and compilers
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    apt-utils \
    zip \
    unzip \
    git \
    curl \
    gawk \
    wget \
    cmake \
    gcc \
    gfortran \
    apt-transport-https \
    ca-certificates \
    software-properties-common

# install libraries and dependencies needed for python/ R environments
RUN apt-get install -y \
    zlib1g-dev \
    libbz2-dev \
    libgsl0-dev \
    libncurses5-dev \
    libmysqlclient-dev \
    default-libmysqlclient-dev \
    autotools-dev \
    automake \
    redis-tools \
    libblas-dev \
    liblapack-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libjpeg-dev \
    zlib1g-dev \
    libfontconfig1-dev \
    libpq-dev \
    libudunits2-dev \
    dirmngr \
    libglu1-mesa-dev \
    freeglut3-dev \
    mesa-common-dev \
    libgdal-dev

############################################################
# Install and configure python environment
# version=3.9

# Install py38 from deadsnakes repository
# install pip3
# install python packages + lithops
RUN add-apt-repository ppa:deadsnakes/ppa \
    && apt-get install -y python3.8 \
    && apt-get install -y python3-pip \
    && pip3 install --no-cache --upgrade pip wheel setuptools \
    && pip install --no-cache --upgrade \
        numpy \
        pandas \
        lithops \
        cloudpickle==2.1.0 \
        redis \
        boto3 \
        "lithops[aws]"


# link python to python3 binary
RUN ln -s /usr/bin/python3 /usr/bin/python

############################################################
# Install and configure R environment
# Version=latest available

# install r-base
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
    apt-get install --yes --no-install-recommends r-base  \
	  build-essential && \
	  apt-get update --fix-missing

# install R packages
RUN Rscript -e "install.packages('plyr', dependencies=TRUE)" \
    && Rscript -e "install.packages('tidyverse', dependencies=TRUE)" \
    && Rscript -e "install.packages('stringi', dependencies=TRUE)" \
    && Rscript -e "install.packages('plotly', dependencies=TRUE)" \
    && Rscript -e "install.packages('rgl', dependencies=TRUE)" \
    && Rscript -e "install.packages('viridis', dependencies=TRUE)"


############################################################
# Install packages requiring compilation
# List: sratoolkit.v.3.0.0,

# Download, compile and install SRA toolkit
ARG FUNCTION_DIR
ARG BINS_DIR
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz \
	&& tar zxvf sratoolkit.3.0.0-ubuntu64.tar.gz  \
    && cd sratoolkit.3.0.0-ubuntu64/bin \
    && mv * ${BINS_DIR} \
    && cd -

############################################################
# Finalize

# clean up
RUN apt-get autoremove && apt-get clean

# Update permissions for binaries and append location to PATH
RUN chmod -R ugoa+rwx ${BINS_DIR}
ENV PATH="${BINS_DIR}:${PATH}"
ENV HOME=/root
