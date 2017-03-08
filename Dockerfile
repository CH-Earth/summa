FROM phusion/baseimage:0.9.19

USER root

RUN apt-get update && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    libnetcdff-dev \
    netcdf-bin \
    liblapack-dev \
    gfortran-6 \
    gcc && apt-get clean

ENV FC gfortran-6
ENV F_MASTER /code
ENV NCDF_PATH /usr
ENV LD_LIBRARY_PATH ${NCDF_PATH}/lib

WORKDIR /code
ADD . /code

RUN make -C build/ -f Makefile
