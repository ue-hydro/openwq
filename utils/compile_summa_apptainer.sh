#! /bin/bash
# File: compile_summa_apptainer.sh
# Authour: Kyle Klenk
# Date: 2024-10-04

###############################################
# Prerequisites
#  - Install apptainer
#  - Buld the container from the Apptainerfile.def
#     - sudo apptainer build openwq.sif Apptainerfile.def 
#  - Once you have an openwq.sif file, you are ready to run this script.
###############################################

SUMMA_DIR=$(realpath $(pwd)/../../../../../)
OPENWQ_DIR="/code/build/source/openwq/openwq"
apptainer exec --bind $SUMMA_DIR:/code --pwd $OPENWQ_DIR ../openwq.sif \
    /bin/bash -c "mkdir build && cd build && cmake .. && make"
