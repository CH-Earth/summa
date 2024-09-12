FROM ubuntu:22.04

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
    libboost-all-dev \
    hdf5-tools \
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

# Install Armadillo in conventional linux directory (/opt)
WORKDIR /opt
RUN wget https://sourceforge.net/projects/arma/files/armadillo-10.5.0.tar.xz
RUN tar -xf armadillo-10.5.0.tar.xz
WORKDIR /opt/armadillo-10.5.0
RUN cmake . -D DETECT_HDF5=true -DCMAKE_C_FLAGS="-DH5_USE_110_API"
RUN make install
# Enable HDF5 Support
RUN sed -i '121s/^\/\/ #define ARMA_USE_HDF5/#define ARMA_USE_HDF5/' /usr/include/armadillo_bits/config.hpp

# Set Entry Point
WORKDIR /code

