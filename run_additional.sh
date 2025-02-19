#!/bin/bash

BUILD_DIR=$(pwd)/build

# provide existing directory to run more files for
RUN_DIR=runs/${1}/${2}

cd ${RUN_DIR}

# addqueue -c "1 day" ${BUILD_DIR}/compute_substructure ${2} 1
# addqueue -c "1 day" ${BUILD_DIR}/compute_EEC ${2} 1

for VARIABLE in {31..90}
do
    addqueue -c "1 day" ${BUILD_DIR}/compute_EEC params_EEC.ini $VARIABLE
    # addqueue -c "1 day" ${BUILD_DIR}/compute_substructure params_sub.ini $VARIABLE
done

#addqueue -c "1 day" ${BUILD_DIR}/compute_substructure ${2}

# if [ "${1}"=="substructure" ]; then
#     addqueue -c "4 hours" ${BUILD_DIR}/compute_substructure ${2}
# elif [ "${1}"=="EEC" ]; then
#     addqueue -c "4 hours" ${BUILD_DIR}/compute_EEC ${2}
# fi
