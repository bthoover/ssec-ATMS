#! /bin/sh
################################################################################
#
# Define inputs
#
WORKING_DIR=`pwd`
#
RATMANL_YYYY=2020
RATMANL_MM=04
RATMANL_DD=17
RATMANL_CYC=t12
#
ATMF006_YYYY=2020
ATMF006_MM=04
ATMF006_DD=17
ATMF006_CYC=t06
#
PYTHON_PROG=Compute_Shannon_Entropy.py
#
################################################################################
#
# Generate inputs
#
RATMANL_DATE=${RATMANL_YYYY}${RATMANL_MM}${RATMANL_DD}
#
RATMANL_DATA_DIR=/data/users/bhoover/ATMS/GDAS_ensemble/native_to_prs_levs/
ATMF006_DATA_DIR=/data/users/bhoover/ATMS/GDAS_ensemble/native_to_prs_levs/
RATMANL_PREFIX=gdas.${RATMANL_CYC}.ratmanl_plev.
ATMF006_PREFIX=gdas.${ATMF006_CYC}.atmf006_plev.
VARNAME=U
TARGET_DX=1.0
NC_OUTFILE_NAME=gdas.shannon_entropy.${VARNAME}.ratmanl.${RATMANL_CYC}.atmf006.${ATMF006_CYC}.${RATMANL_DATE}.nc
#
# Run exec script
#
bash ${WORKING_DIR}/exec_Compute_Shannon_Entropy.sh ${RATMANL_DATA_DIR} ${ATMF006_DATA_DIR} ${RATMANL_PREFIX} ${ATMF006_PREFIX} ${VARNAME} ${TARGET_DX} ${NC_OUTFILE_NAME} ${PYTHON_PROG}
#
################################################################################
