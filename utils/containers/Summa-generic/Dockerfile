FROM ubuntu:latest

WORKDIR /code
RUN apt-get update -y && \
    apt-get upgrade -y && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y software-properties-common \
    libnetcdf-dev \
    libnetcdff-dev \
    liblapack-dev \
    libopenblas-dev \
    cmake \
    g++ \
    git \
    libssl-dev \
    make \
    gfortran \
    wget \
    python3-pip \
    valgrind \
    gdb &&\
    apt-get autoclean

# Install Sundials
WORKDIR /opt
RUN wget https://github.com/LLNL/sundials/archive/refs/tags/v7.1.1.tar.gz
RUN tar -xzf v7.1.1.tar.gz
WORKDIR /opt/sundials-7.1.1
RUN mkdir build/
WORKDIR /opt/sundials-7.1.1/build
RUN cmake ../ -DBUILD_FORTRAN_MODULE_INTERFACE=ON \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_INSTALL_PREFIX=/usr/local/sundials
RUN make -j 4
RUN make install


WORKDIR /code
