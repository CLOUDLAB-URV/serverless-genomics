# Dockerfile for genomics with Kubernetes backend
ARG FUNCTION_DIR="/function"
ARG BINS_DIR="/function/bin"
FROM python:3.8-slim-buster
ARG FUNCTION_DIR
ARG BINS_DIR

# Install build dependencies
RUN apt-get update && apt-get install -y \
        zip \
        g++ \
	gcc \
	make \
	cmake \
	autoconf \
	automake \
	unzip \
	perl \
    git \
	wget \
	libssl-dev \
	libncurses5-dev \
	zlib1g-dev \
	libxslt-dev \
	libxml2-dev \
	zlib1g-dev \
    liblzma-dev \
	libbz2-dev \
        && rm -rf /var/lib/apt/lists/*

# Create function work directories
RUN mkdir -p ${FUNCTION_DIR}
RUN mkdir -p ${BINS_DIR}

# Compile gem-mapper
RUN git clone --recursive https://github.com/smarco/gem3-mapper.git \
    && cd gem3-mapper \
    && sh configure \
    && make \
    && cd .. \
    && mv gem3-mapper/bin/* ${BINS_DIR}
     
# Compile samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2 -O samtools.tar.bz2 \
    && tar -xjvf samtools.tar.bz2 \
    && cd samtools-1.15 \
    && ./configure prefix=${FUNCTION_DIR} \
    && make \
    && make install \
    && cd ..

# Download, compile and install gztool
RUN wget https://github.com/circulosmeos/gztool/archive/refs/tags/v1.4.2.zip \
    && unzip v1.4.2.zip \
    && cd gztool-1.4.2 \
    && gcc -O3 -o gztool gztool.c -lz -lm \
    && cp gztool ${BINS_DIR}

# Download, compile and install fastq-dump
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz \
	&& tar zxvf sratoolkit.3.0.0-ubuntu64.tar.gz  \
    && cd /sratoolkit.3.0.0-ubuntu64/bin \
    && mv * ${BINS_DIR}

#Install dependencies for fastq-dump
RUN apt-get update  \
    && apt-get --quiet install --yes libxml-libxml-perl

# Add sinple
COPY SiNPle-0.5 ${BINS_DIR}/SiNPle-0.5
COPY com ${BINS_DIR}/com
        
# Install function dependencies
RUN pip install --upgrade setuptools six pip \
    && pip install --no-cache-dir \
        flask \
        pika \
        ibm-cos-sdk \
        redis \
        gevent \
        requests \
        PyYAML \
        kubernetes \
        numpy \
        cloudpickle \
        ps-mem \
        boto3 \
        httplib2 \
        scipy \
        pandas \
        kafka-python \
        pyarrow \
        fastparquet \
        tblib

# Add mpileup merge script
COPY mpileup_merge_reducev3.sh ${BINS_DIR}/mpileup_merge_reducev3.sh
COPY mpileup_merge_reducev3_nosinple.sh ${BINS_DIR}/mpileup_merge_reducev3_nosinple.sh
COPY parse_gem_maxindex_minimapfile_stdin_v2.sh ${BINS_DIR}/parse_gem_maxindex_minimapfile_stdin_v2.sh
COPY map_index_and_filter_map_file_cmd_awsruntime.sh ${BINS_DIR}/map_index_and_filter_map_file_cmd_awsruntime.sh
COPY map_file_index_correction.sh ${BINS_DIR}/map_file_index_correction.sh
COPY map_file_index_correction_v2.sh ${BINS_DIR}/map_file_index_correction_v2.sh
COPY map_file_index_correction_v3.sh ${BINS_DIR}/map_file_index_correction_v3.sh
COPY gempileup_v7.sh ${BINS_DIR}/gempileup_v7.sh
COPY gempileup_merge.sh ${BINS_DIR}/gempileup_merge.sh
COPY gempileup_run.sh ${BINS_DIR}/gempileup_run.sh

#INDEX CORRECTION
COPY binary_reducer.sh ${BINS_DIR}/binary_reducer.sh
COPY filter_merged_index.sh ${BINS_DIR}/filter_merged_index.sh
COPY merge_gem_alignment_metrics.sh ${BINS_DIR}/merge_gem_alignment_metrics.sh

ENV PYTHONUNBUFFERED TRUE

# Include global arg in this stage of the build
ARG FUNCTION_DIR
ARG BINS_DIR

RUN apt update && apt install \
    gawk \
    && rm -rf /var/lib/apt/lists/* \
    && apt-cache search linux-headers-generic

# Copy Lithops proxy and lib to the container image.
ENV APP_HOME /lithops

# Set working directory to function root directory
WORKDIR $APP_HOME

# Update permissions for binaries and append location to PATH
RUN chmod -R ugoa+rwx ${BINS_DIR}
ENV PATH="${BINS_DIR}:${PATH}"

COPY lithops_k8s.zip .
RUN unzip lithops_k8s.zip && rm lithops_k8s.zip
