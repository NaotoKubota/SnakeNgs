#####################################################################
# Dockerfile to build container images for ngsFetch
# Based on python 3.11-buster
#####################################################################

FROM python:3.11-buster

# File Author / Maintainer
LABEL maintainer="Naoto Kubota <naotok@ucr.edu>"

# Install dependencies first
RUN apt-get -qq update && \
	apt-get -qq -y install \
	make wget less vim gcc git jq aria2

# mkdir /src
RUN mkdir /src

# git clone ZhengLab RNA-seq pipeline
RUN git clone https://github.com/NaotoKubota/SnakeNgs.git -b v0.1.1 && \
	mv SnakeNgs /src/SnakeNgs && \
	chmod a+x /src/SnakeNgs/*.bash

# Set PATH
ENV PATH $PATH:/src:/src/SnakeNgs

# Upgrade pip
RUN pip install --upgrade pip

# Install ffq (version 0.3.1)
RUN pip install ffq==0.3.1

# Set working directory
WORKDIR /home

# bash
CMD ["bash"]
