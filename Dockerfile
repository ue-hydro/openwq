
FROM ubuntu:24.04
WORKDIR /code

# Update package manager
RUN apt-get -y update && apt-get upgrade -y && \
    DEBIAN_FRONTEND="noninteractive" \
    apt-get -y install g++ \
    gdb \
    gdbserver \
    cmake \
    wget \
    xz-utils \
    liblapack-dev \
    libblas-dev \
    libopenblas-dev \
    libboost-all-dev \
    hdf5-tools \
    software-properties-common \
    libnetcdf-dev \
    libnetcdff-dev \
    liblapack-dev \
    gfortran \
    git \
    nano

# Install PHREEQC-RM
WORKDIR /opt
RUN git clone https://github.com/usgs-coupled/phreeqcrm_use_cmake.git phreeqcrm -b modern
WORKDIR /opt/phreeqcrm
RUN cmake . -P build_phreeqcrm.cmake -DCMAKE_INSTALL_PREFIX=/usr/local/phreeqcrm


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
RUN wget http://sourceforge.net/projects/arma/files/armadillo-10.6.0.tar.xz
RUN tar -xf armadillo-10.6.0.tar.xz
WORKDIR /opt/armadillo-10.6.0
RUN cmake . -D DETECT_HDF5=true -DCMAKE_C_FLAGS="-DH5_USE_110_API"
RUN make install
# Enable HDF5 Support
RUN sed -i '121s/^\/\/ #define ARMA_USE_HDF5/#define ARMA_USE_HDF5/' /usr/include/armadillo_bits/config.hpp
# Set Entry Point
WORKDIR /code