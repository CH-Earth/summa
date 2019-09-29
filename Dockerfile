FROM ubuntu:xenial

# install only the packages that are needed
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    software-properties-common python-software-properties \
    ca-certificates \
    git \
    make \
    libnetcdff-dev \
    liblapack-dev \
    vim

# install gfortran-6
RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y \
    && apt-get update \
    && apt-get install -y --no-install-recommends gfortran-6 \
    && apt-get clean

# set environment variables for docker build
ENV F_MASTER /code
ENV FC gfortran
ENV FC_EXE gfortran
ENV INCLUDES -I/usr/include
ENV LIBRARIES '-L/usr/lib -lnetcdff -llapack -lblas'

# add code directory
WORKDIR /code
ADD . /code

# fetch tags and build summa
RUN git fetch --tags && make -C build/ -f Makefile

# run summa when running the docker image
WORKDIR bin
ENTRYPOINT ["./summa.exe"]
