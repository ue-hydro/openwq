#! /bin/bash

# Specify the path to where Armadillo will be installed
ARMADILLO_VERSION="10.6.0"
ROOT_DIR="$(pwd)/armadillo-${ARMADILLO_VERSION}"
INSTALL_DIR="${ROOT_DIR}/inst"
BUILD_DIR="${ROOT_DIR}/build"

# Check if the install directory is specified
if [ -z "$INSTALL_DIR" ]; then
  echo "Please specify the path to where Armadillo will be installed"
  exit 1
fi

# Create the directories for installation
mkdir -p $ROOT_DIR
mkdir -p $INSTALL_DIR
mkdir -p $BUILD_DIR

# Check if the tar file is already downloaded
if [ -f "${BUILD_DIR}/armadillo-${ARMADILLO_VERSION}.tar.xz" ]; then
  echo "Armadillo tar file already downloaded. Skipping Installation."
  echo "To re-install, please remove the armadillo-${ARMADILLO_VERSION} " \
       "directory and re-run the script"
  exit 0
fi

# Download the Armadillo source code
cd $BUILD_DIR
wget http://sourceforge.net/projects/arma/files/armadillo-${ARMADILLO_VERSION}.tar.xz
tar -xf armadillo-${ARMADILLO_VERSION}.tar.xz
cd armadillo-${ARMADILLO_VERSION}
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DDETECT_HDF5=ON \
    -DCMAKE_C_FLAGS="-DH5_USE_110_API" -DCMAKE_CXX_FLAGS="-" .
make install









