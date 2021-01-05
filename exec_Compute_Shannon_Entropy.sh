#! /bin/sh
################################################################################
#
# Read inputs
#
RATMANL_DATA_DIR=${1}
ATMF006_DATA_DIR=${2}
RATMANL_PREFIX=${3}
ATMF006_PREFIX=${4}
VARNAME=${5}
TARGET_DX=${6}
NC_OUTFILE_NAME=${7}
PYTHON_PROG=${8}
#
################################################################################
#
# Load conda envrionment
#
source activate bhoover-GDAS_information_content_fastkde
#
################################################################################
#
# Run PYTHON_PROG
#
python3 ${PYTHON_PROG} << EOF
${RATMANL_DATA_DIR}
${ATMF006_DATA_DIR}
${RATMANL_PREFIX}
${ATMF006_PREFIX}
${VARNAME}
${TARGET_DX}
${NC_OUTFILE_NAME}
EOF
#
################################################################################
