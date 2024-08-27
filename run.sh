#!/bin/bash

BUILD_DIR=$(pwd)/build

# generate a 6-digit code based on the Unix timestamp and name the run directory after that
RUN_DIR=runs/${1}/run_$(date +"%s" | tail -c 7)

# create the run directory if it doesn't exist
mkdir -p $RUN_DIR

# copy the parameter file to it
cp ${2} ${RUN_DIR}/${2}

cd ${RUN_DIR}

if [ ${1}=="substructure" ]
then
    ${BUILD_DIR}/compute_substructure ${2}
#    addqueue -c "4 hours" ${BUILD_DIR}/compute_substructure ${2}
elif [ ${1}=="EEC" ]
then
    ${BUILD_DIR}/compute_EEC ${2}
#    addqueue -c "4 hours" ${BUILD_DIR}/compute_EEC ${2}
fi