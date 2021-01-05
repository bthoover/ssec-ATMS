#! /bin/bash -f

########################################################################################################################
#
# SET RUNTIME NCL ENVIRONMENT
#
module load ncl
#setenv NCARG_ROOT /usr/local/ncl-6.3.0
#setenv PATH $NCARG_ROOT/bin:$PATH
#
########################################################################################################################
NATIVE_DATA_NAME=${1}
LAT_VAR_NAME=${2}
LON_VAR_NAME=${3}
VC1_VAR_NAME=${4}
VC2_VAR_NAME=${5}
U_VAR_NAME=${6}
V_VAR_NAME=${7}
T_VAR_NAME=${8}
Q_VAR_NAME=${9}
PS_VAR_NAME=${10}
ZS_VAR_NAME=${11}
OUTF_NAME=${12}
NCL_SCRIPT=${13}
########################################################################################################################
#
# Define environment variables for NCL_SCRIPT
#
export native_data_name="${NATIVE_DATA_NAME}"
export LAT_var_name="${LAT_VAR_NAME}"
export LON_var_name="${LON_VAR_NAME}"
export VC1_var_name="${VC1_VAR_NAME}"
export VC2_var_name="${VC2_VAR_NAME}"
export U_var_name="${U_VAR_NAME}"
export V_var_name="${V_VAR_NAME}"
export T_var_name="${T_VAR_NAME}"
export Q_var_name="${Q_VAR_NAME}"
export PS_var_name="${PS_VAR_NAME}"
export ZS_var_name="${ZS_VAR_NAME}"
export outf_name="${OUTF_NAME}"
#
# Run NCL_SCRIPT
#
ncl ${NCL_SCRIPT}
#
########################################################################################################################
