# use the zesty distribution, which has gcc-6
FROM ubuntu:zesty

# install only the packages that are needed
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    git \
    make \
    gfortran-6 \
    libnetcdff-dev \
    liblapack-dev \
    && apt-get clean

# set environment variables for docker build
ENV F_MASTER /code
ENV FC gfortran
ENV FC_EXE gfortran
ENV FC_ENV gfortran-6-docker

# add code directory
WORKDIR /code
ADD . /code

# fetch tags and build summa
RUN git fetch --tags && make -C build/ -f Makefile

# run summa when running the docker image
WORKDIR bin
ENTRYPOINT ["./summa.exe"]
