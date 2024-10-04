#! /bin/bash
# File: compile_summa_apptainer.sh
# Authour: Kyle Klenk
# Date: 2024-10-04

###############################################
# Prerequisites
#  - Install apptainer
#  - Buld the container from the Apptainerfile.def
#     - sudo apptainer build openwq.sif Apptainerfile.def 
#  - Use the container to compile SUMMA-OPENWQ
#  - Once you have SUMMA-OPENWQ compiled, you are ready to run this script.
###############################################

SUMMA_DIR=$(realpath $(pwd)/../../../../../)
OPENWQ_DIR="/code/build/source/openwq/openwq"

FILE_MANAGER="/code/synthetic_tests/2_nrTrans_instS_PorMedia/summa/summa/SUMMA/summa_fileManager_OpenWQ_systheticTests_BGQ.txt"
MASTER_JSON="/code/synthetic_tests/2_nrTrans_instS_PorMedia/summa/openWQ_master.json"

apptainer exec --bind $SUMMA_DIR:/code --pwd $OPENWQ_DIR/bin \
    --env master_json=$MASTER_JSON ../openwq.sif \
    ./summa_openwq_debug -g 1 1 -m $FILE_MANAGER
