#! /bin/csh -f

set YYYY = 2020
set MM   = 04
foreach DD ( 17 )
#foreach HH ( 00 06 12 18 )
foreach HH ( 06 12 )
foreach XXX ( 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 )

#foreach DD ( 17 )
#foreach HH ( 00 )
#foreach XXX ( 001 )

set WORKDIR = /scratch/users/bhoover/OUTNEMSIO/enkfgdas.${YYYY}${MM}${DD}/${HH}/mem${XXX}
set NATIVE_DATA_NAME = ${WORKDIR}/gdas.t${HH}.atmf006.mem${XXX}.nc
set LAT_VAR_NAME = lat
set LON_VAR_NAME = lon
set VC1_VAR_NAME = VCOORD1
set VC2_VAR_NAME = VCOORD2
set U_VAR_NAME   = U
set V_VAR_NAME   = V
set T_VAR_NAME   = T
set Q_VAR_NAME   = Q
set PS_VAR_NAME  = sp
set ZS_VAR_NAME  = sz
set OUTF_NAME = ${WORKDIR}/gdas.t${HH}.atmf006_plev.mem${XXX}.nc
set NCL_SCRIPT = compute_gdas_pressure_level_data.ncl
./exec_compute_gdas_pressure_level_data.bash ${NATIVE_DATA_NAME} ${LAT_VAR_NAME} ${LON_VAR_NAME} ${VC1_VAR_NAME} ${VC2_VAR_NAME} ${U_VAR_NAME} ${V_VAR_NAME} ${T_VAR_NAME} ${Q_VAR_NAME} ${PS_VAR_NAME} ${ZS_VAR_NAME} ${OUTF_NAME} ${NCL_SCRIPT}

set WORKDIR = /scratch/users/bhoover/OUTNEMSIO/enkfgdas.${YYYY}${MM}${DD}/${HH}/mem${XXX}
set NATIVE_DATA_NAME = ${WORKDIR}/gdas.t${HH}.ratmanl.mem${XXX}.nc
set LAT_VAR_NAME = lat
set LON_VAR_NAME = lon
set VC1_VAR_NAME = VCOORD1
set VC2_VAR_NAME = VCOORD2
set U_VAR_NAME   = U
set V_VAR_NAME   = V
set T_VAR_NAME   = T
set Q_VAR_NAME   = Q
set PS_VAR_NAME  = sp
set ZS_VAR_NAME  = sz
set OUTF_NAME = ${WORKDIR}/gdas.t${HH}.ratmanl_plev.mem${XXX}.nc
set NCL_SCRIPT = compute_gdas_pressure_level_data.ncl
./exec_compute_gdas_pressure_level_data.bash ${NATIVE_DATA_NAME} ${LAT_VAR_NAME} ${LON_VAR_NAME} ${VC1_VAR_NAME} ${VC2_VAR_NAME} ${U_VAR_NAME} ${V_VAR_NAME} ${T_VAR_NAME} ${Q_VAR_NAME} ${PS_VAR_NAME} ${ZS_VAR_NAME} ${OUTF_NAME} ${NCL_SCRIPT}

end
end
end
