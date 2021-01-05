load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; DESCRIPTION
;
; 1. Computes post-processed variables (geopotential height, vorticity, divergence)
;
; 2. Converts variables to standard pressure surfaces
;
; 3. Outputs pressure-level data to file
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; INPUT PARAMETERS
;
; native_data_name................................. input file-name for data on native grid (hybrid sigma-pressure coordinate)
; LAT_var_name..................................... Name of latitude variable
; LON_var_name..................................... Name of longitude variable
; VC1_var_name..................................... Name of variable coordinate-1 variable (vertical coordinate ak(p) value)
; VC2_var_name..................................... Name of variable coordinate-2 variable (vertical coordinate bk(p) value)
; U_var_name....................................... Name of zonal wind variable
; V_var_name....................................... Name of merid wind variable
; T_var_name....................................... Name of temperature variable
; Q_var_name....................................... Name of mixing ratio variable
; PS_var_name...................................... Name of surface pressure variable
; ZS_var_name...................................... Name of surface elevation variable
; outf_name........................................ output file-name on new grid (standard pressure surfaces)
native_data_name = getenv("native_data_name");
LAT_var_name = getenv("LAT_var_name");
LON_var_name = getenv("LON_var_name");
VC1_var_name = getenv("VC1_var_name");
VC2_var_name = getenv("VC2_var_name");
U_var_name = getenv("U_var_name");
V_var_name = getenv("V_var_name");
T_var_name = getenv("T_var_name");
Q_var_name = getenv("Q_var_name");
PS_var_name = getenv("PS_var_name");
ZS_var_name = getenv("ZS_var_name");
outf_name = getenv("outf_name");
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("EXTRACTING DATA FROM NATIVE GRID");
;
; Open input file
;
native_file = addfile(native_data_name,"r");
;
; Load data: "n" indicates native data (includes single dimensions, single format)
;
nLAT = native_file->$LAT_var_name$(:);..................................................... Latitude [lat]
nLON = native_file->$LON_var_name$(:);..................................................... Longitude [lon]
nVC1 = native_file->$VC1_var_name$(:);..................................................... Variable cooordinate-1 [lev+1]
nVC2 = native_file->$VC2_var_name$(:);..................................................... Variable cooordinate-2 [lev+1]
nU = native_file->$U_var_name$(:,:,:,:);................................................... Zonal wind [time,lev,lat,lon]
nV = native_file->$V_var_name$(:,:,:,:);................................................... Merid wind [time,lev,lat,lon]
nT = native_file->$T_var_name$(:,:,:,:);................................................... Temperature [time,lev,lat,lon]
nQ = native_file->$Q_var_name$(:,:,:,:);................................................... Mixing ratio [time,lev,lat,lon]
nPS = native_file->$PS_var_name$(:,:,:);................................................... Surface pressure [time,lev,lat,lon]
nZS = native_file->$ZS_var_name$(:,:,:);................................................... Surface elevation [time,lev,lat,lon]
;
; Remove single dimensions from variables and change to double format
;
LAT = todouble(rm_single_dims(nLAT))
LON = todouble(rm_single_dims(nLON))
VC1 = todouble(rm_single_dims(nVC1))
VC2 = todouble(rm_single_dims(nVC2))
U   = todouble(rm_single_dims(nU))
V   = todouble(rm_single_dims(nV))
T   = todouble(rm_single_dims(nT))
Q   = todouble(rm_single_dims(nQ))
PS  = todouble(rm_single_dims(nPS))
ZS  = todouble(rm_single_dims(nZS))
delete(nLAT)
delete(nLON)
delete(nVC1)
delete(nVC2)
delete(nU)
delete(nV)
delete(nT)
delete(nQ)
delete(nPS)
delete(nZS)
;
; Define dimension sizes (use T as an example 3D array)
;
dSizes = dimsizes(T) ;..................................................................... Dimension sizes
nlev = dSizes(0) ;......................................................................... Number of vertical levels (dashed levels)
nlev2 = dSizes(0)+1 ;...................................................................... Number of vertical levels (solid levels)
nlat = dSizes(1) ;......................................................................... Number of latitude points
nlon = dSizes(2) ;......................................................................... Number of longitude points
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("COMPUTE PRESSURE AND VIRTUAL TEMPERATURE")
;
; Compute pressure at each [lev,lat,lon] gridpoint from hybrid sigma-pressure coordinates.
; the vertical coordinate coefficients are generated on "solid" levels, while we need the
; pressure computed on "dashed" levels, which requires some careful interpolation. Using
; Allen Lenzen's code example to accomplish this.
;
P_sol = new((/nlev2,nlat,nlon/),typeof(T)) ;............................................... Pressure (solid levels) [lev,lat,lon]
P_sol_k = new((/nlev2,nlat,nlon/),typeof(T)) ;............................................. Pressure (solid levels, interpolated) [lev,lat,lon]
P_k = new((/nlev,nlat,nlon/),typeof(T)) ;.................................................. Pressure (solid levels, interpolated) [lev,lat,lon]
P_dash = new((/nlev,nlat,nlon/),typeof(T)) ;............................................... Pressure (dashed levels) [lev,lat,lon]
;
cp = 1004.6 ;.............................................................................. Heat capacity of air at constant pressure
r = 287.05 ;............................................................................... Dry gas constant
capa = r / cp ;............................................................................ Kappa
grav = 9.80655 ;........................................................................... Gravity
;
; ALL ELEMENT-WISE OPERATIONS IN HORIZONTAL (j,i)-SPACE DONE SIMULTANEOUSLY: FAST 
; REMOVING P_dash() OPERATIONS BECAUSE THEY SEEM TO BE REDUNDANT
;
do k = 0, nlev2-1
    P_sol(k,:,:) = (VC1(k)*0.01) + (VC2(k)*PS(:,:))
    P_sol_k(k,:,:) = (P_sol(k,:,:)/1000.0)^capa
end do
do k = 0, nlev-1
    P_k(k,:,:) = (((P_sol(k+1,:,:)*P_sol_k(k+1,:,:))-(P_sol(k,:,:)*P_sol_k(k,:,:))) / (1. + capa))/(P_sol(k+1,:,:)-P_sol(k,:,:))
    P_dash(k,:,:) = (1000.0)*P_k(k,:,:)^(1./capa)
end do
;
; Compute virtual temperature as Tv = T * (1.0 + 0.61*W)
;
TV = new(dimsizes(T),typeof(T)) ;.......................................................... Virtual temperate [lev,lat,lon]
TV = T * (1.0 + 0.61*Q)
;
; Computing surface (virtual) temperature by following dry adiabatic lapse rate (Poisson's
; Equation), assuming well-mixed moisture in lowest layer (same Q value as layer above 
; surface)
;
TS = new(dimsizes(PS),typeof(T)) ;......................................................... Surface temperature [lat,lon]: APPROXIMATION
TSV = new(dimsizes(PS),typeof(T)) ;........................................................ Surface virtual temperature [lat,lon]: APPROXIMATION
T_bottom = rm_single_dims(T(nlev-1,:,:)) ;................................................. Temperature at bottom level
P_bottom = rm_single_dims(P_dash(nlev-1,:,:)) ;............................................ Pressure at bottom level
Q_bottom = rm_single_dims(Q(nlev-1,:,:)) ;................................................. Mixing ratio at bottom level
TS = (T_bottom*(PS/P_bottom)^capa)
TSV = TS  * (1.0 + 0.61*Q_bottom)
; 
; Create P_ARRAY, TV_ARRAY that are reversed in k-dimension and have PS, TS as bottom-lev
;
P_ARRAY = new((/nlev+1,nlat,nlon/),typeof(T)) ;............................................ Pressure array for geopotential height calculation
TV_ARRAY = new((/nlev+1,nlat,nlon/),typeof(T)) ;........................................... Virtual temperature array for geopotential height calcuation
P_ARRAY(0,:,:) = PS
TV_ARRAY(0,:,:) = TSV
P_ARRAY(1:nlev,:,:) = P_dash(::-1,:,:)
TV_ARRAY(1:nlev,:,:) = TV(::-1,:,:)
ZS_ARRAY = ZS
;
; Name dimensions
;
P_ARRAY!0  = "level"
P_ARRAY!1  = "latitude"
P_ARRAY!2  = "longitude"
TV_ARRAY!0 = "level"
TV_ARRAY!1 = "latitude"
TV_ARRAY!2 = "longitude"
ZS_ARRAY!0 = "latitude"
ZS_ARRAY!1 = "longitude"
;
; Re-arrange into [lon,lat,lev] format, and [lon,lat] format for ZS_ARRAY
;
nP_ARRAY = P_ARRAY(longitude|:,latitude|:,level|:)
nTV_ARRAY = TV_ARRAY(longitude|:,latitude|:,level|:)
nZS_ARRAY = ZS_ARRAY(longitude|:,latitude|:)
delete(P_ARRAY)
delete(TV_ARRAY)
delete(ZS_ARRAY)
P_ARRAY = nP_ARRAY
TV_ARRAY = nTV_ARRAY
ZS_ARRAY = nZS_ARRAY
delete(nP_ARRAY)
delete(nTV_ARRAY)
delete(nZS_ARRAY)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("COMPUTE GEOPOTENTIAL HEIGHT")
;
; Compute geopotential height using hydro()
;
gh = hydro(P_ARRAY,TV_ARRAY,ZS_ARRAY) ;.................................................... Geopotential height [lon,lat,lev]
;
; Pull out interior levels (non-surface) to G
G = new((/nlon,nlat,nlev/),typeof(T)) ;.................................................... Geopotential height [lon,lat,lev]: Surface-level removed
do k = 1,nlev2-1
    G(:,:,k-1) = rm_single_dims(gh(:,:,k))
end do
G!0 = "longitude"
G!1 = "latitude"
G!2 = "level"
;
; Reverse k-dimension back to native top-to-bottom format
;
G = G(:,:,::-1)
;
; Return G to native [lev,lat,lon] format
;
nG = G(level|:,latitude|:,longitude|:)
delete(G)
G = nG
delete(nG)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("INTERPOLATE MODEL DATA TO STANDARD PRESSURE LEVELS")
;
; Use interp2p_n_Wrap() to move vertical coordinate to standard pressure levels
;
Pres_Levs = (/1000.0,975.0,950.0,925.0,900.0,850.0,800.0,750.0,700.0,650.0,600.0,550.0,\ ;. Pressure levels (hPa) to interpolate to
              500.0,450.0,400.0,350.0,300.0,250.0,200.0,150.0,100.0,50.0,25.0,10.0/)
UP = int2p_n_Wrap(P_dash,U,Pres_Levs,0,0) ;................................................ Zonal wind interpolated to Pres_Levs [lev,lat,lon]
VP = int2p_n_Wrap(P_dash,V,Pres_Levs,0,0) ;................................................ Merid wind interpolated to Pres_Levs [lev,lat,lon]
TP = int2p_n_Wrap(P_dash,T,Pres_Levs,0,0) ;................................................ Temperature interpolated to Pres_Levs [lev,lat,lon]
QP = int2p_n_Wrap(P_dash,Q,Pres_Levs,0,0) ;................................................ Mixing ratio interpolated to Pres_Levs [lev,lat,lon]
GP = int2p_n_Wrap(P_dash,G,Pres_Levs,0,0) ;................................................ Geopotential height interpolated to Pres_Levs [lev,lat,lon]
;
; Formalize coordinate variables
;
UP!0 = "level"
UP!1 = "latitude"
UP!2 = "longitude"
VP!0 = "level"
VP!1 = "latitude"
VP!2 = "longitude"
TP!0 = "level"
TP!1 = "latitude"
TP!2 = "longitude"
QP!0 = "level"
QP!1 = "latitude"
QP!2 = "longitude"
GP!0 = "level"
GP!1 = "latitude"
GP!2 = "longitude"
PS!0 = "latitude"
PS!1 = "longitude"
ZS!0 = "latitude"
ZS!1 = "longitude"
TS!0 = "latitude"
TS!1 = "longitude"
LAT!0 = "latitude"
LON!0 = "longitude"
Pres_Levs!0 = "level"
UP&level = Pres_Levs
VP&level = Pres_Levs
TP&level = Pres_Levs
QP&level = Pres_Levs
GP&level = Pres_Levs
Pres_Levs&level = Pres_Levs
UP&latitude = LAT
VP&latitude = LAT
TP&latitude = LAT
QP&latitude = LAT
GP&latitude = LAT
PS&latitude = LAT
ZS&latitude = LAT
TS&latitude = LAT
LAT&latitude = LAT
UP&longitude = LON
VP&longitude = LON
TP&longitude = LON
QP&longitude = LON
GP&longitude = LON
PS&longitude = LON
ZS&longitude = LON
TS&longitude = LON
LON&longitude = LON
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("COMPUTE VORTICITY, DIVERGENCE")
;
; Use uv2vr_cfd(), uv2dv_cfd() to compute vorticity and divergence from UP and VP. This
; uses finite center difference rather than spherical harmonics, but it can operate where
; there are missing values, while spherical harmonics cannot.
;
VOR = uv2vr_cfd(UP,VP,LAT,LON,3) ;......................................................... Vorticity [lev,lat,lon]
DIV = uv2dv_cfd(UP,VP,LAT,LON,3) ;......................................................... Divergence [lev,lat,lon]
;
VOR!0 = UP!0
VOR!1 = UP!1
VOR!2 = UP!2
DIV!0 = UP!0
DIV!1 = UP!1
DIV!2 = UP!2
;
copy_VarCoords(UP,VOR);
copy_VarCoords(UP,DIV);
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("CONVERTING TO 1x1 DEGREE GRID")

LAT_newgrd = ispan(-90,90,1)*1.
LAT_newgrd!0 = LAT!0
LAT_newgrd&latitude=LAT_newgrd
LON_newgrd = ispan(0,359,1)*1.
LON_newgrd!0 = LON!0
LON_newgrd&longitude=LON_newgrd

Opt = True
Opt@SrcFileName = "native_grid_file.nc"
Opt@DstFileName = "newgrd_grid_file.nc"
Opt@WgtFileName = "native_to_newgrd_weights_file.nc"
Opt@SrcForceOverwrite = True
Opt@DstForceOverwrite = True
Opt@SrcLargeFile = True
Opt@DstLargeFile = True
Opt@SrcGridLat = LAT
Opt@SrcGridLon = LON
Opt@DstGridLat = LAT_newgrd
Opt@DstGridLon = LON_newgrd
Opt@SrcRegional = False
Opt@DstRegional = False
Opt@ForceOverwrite = True
system("/bin/rm -f native_grid_file.nc")
system("/bin/rm -f newgrd_grid_file.nc")
system("/bin/rm -f native_to_newgrd_weights_file.nc")
print("... UP ...")
UP_newgrd = ESMF_regrid(UP,Opt)
UP_newgrd!0 = UP!0
UP_newgrd!1 = UP!1
UP_newgrd!2 = UP!2
print("... VP ...")
VP_newgrd = ESMF_regrid(VP,Opt)
VP_newgrd!0 = VP!0
VP_newgrd!1 = VP!1
VP_newgrd!2 = VP!2
print("... TP ...")
TP_newgrd = ESMF_regrid(TP,Opt)
TP_newgrd!0 = TP!0
TP_newgrd!1 = TP!1
TP_newgrd!2 = TP!2
print("... QP ...")
QP_newgrd = ESMF_regrid(QP,Opt)
QP_newgrd!0 = QP!0
QP_newgrd!1 = QP!1
QP_newgrd!2 = QP!2
print("... GP ...")
GP_newgrd = ESMF_regrid(GP,Opt)
GP_newgrd!0 = GP!0
GP_newgrd!1 = GP!1
GP_newgrd!2 = GP!2
print("... VOR ...")
VOR_newgrd = ESMF_regrid(VOR,Opt)
VOR_newgrd!0 = VOR!0
VOR_newgrd!1 = VOR!1
VOR_newgrd!2 = VOR!2
print("... DIV ...")
DIV_newgrd = ESMF_regrid(DIV,Opt)
DIV_newgrd!0 = DIV!0
DIV_newgrd!1 = DIV!1
DIV_newgrd!2 = DIV!2
print("... PS ...")
PS_newgrd = ESMF_regrid(PS,Opt)
PS_newgrd!0 = PS!0
PS_newgrd!1 = PS!1
print("... ZS ...")
ZS_newgrd = ESMF_regrid(ZS,Opt)
ZS_newgrd!0 = ZS!0
ZS_newgrd!1 = ZS!1
print("... TS ...")
TS_newgrd = ESMF_regrid(TS,Opt)
TS_newgrd!0 = TS!0
TS_newgrd!1 = TS!1
;
print(getvardims(UP_newgrd)+" "+dimsizes(UP_newgrd))
print(getvardims(LAT_newgrd)+" "+dimsizes(LAT_newgrd))
print(getvardims(LON_newgrd)+" "+dimsizes(LON_newgrd))
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("OUTPUT TO NETCDF FILE")
;
; Generate output file
;
setfileoption("nc","Format","LargeFile")
system("/bin/rm -f " + outf_name)
o_file = addfile(outf_name,"c") ;.......................................................... Output file handle
fAtt = True ;.............................................................................. File attributes
fAtt@title = "NCL generated netCDF file of FV3GDAS member data"
fAtt@Conventions = "None";
fileattdef(o_file,fAtt);
;
; Define variable dimensions, sizes
;
dNames = getvardims(UP_newgrd) ;.................................................................. Dimension names
dSizes = dimsizes(UP_newgrd) ;.................................................................... Dimension sizes
dUnlim = (/False,False,False/); ;.......................................................... Flags for unlimited-size dimensions
;
; Define new variable names
;
G_var_name = "GHT" ;....................................................................... Variable name for geopotential height
VOR_var_name = "VOR" ;..................................................................... Variable name for vorticity
DIV_var_name = "DIV" ;..................................................................... Variable name for divergence
TS_var_name = "ts" ;....................................................................... Variable name for (approximate) surface temperature
PL_var_name = "plev" ;..................................................................... Variable name for pressure levels
;
; Write to file
;
; 3D variables [lev,lat,lon]
filedimdef(o_file,dNames,dSizes,dUnlim);
filevardef(o_file,U_var_name,typeof(UP_newgrd),getvardims(UP_newgrd));
filevardef(o_file,V_var_name,typeof(VP_newgrd),getvardims(VP_newgrd))
filevardef(o_file,T_var_name,typeof(TP_newgrd),getvardims(TP_newgrd));
filevardef(o_file,Q_var_name,typeof(QP_newgrd),getvardims(QP_newgrd));
filevardef(o_file,G_var_name,typeof(GP_newgrd),getvardims(GP_newgrd));
filevardef(o_file,VOR_var_name,typeof(VOR_newgrd),getvardims(VOR_newgrd));
filevardef(o_file,DIV_var_name,typeof(DIV_newgrd),getvardims(DIV_newgrd));
; 2D variables [lat,lon]
filevardef(o_file,PS_var_name,typeof(PS_newgrd),getvardims(PS_newgrd));
filevardef(o_file,ZS_var_name,typeof(ZS_newgrd),getvardims(ZS_newgrd));
filevardef(o_file,TS_var_name,typeof(TS_newgrd),getvardims(TS_newgrd));
filevardef(o_file,LAT_var_name,typeof(LAT_newgrd),getvardims(LAT_newgrd));
filevardef(o_file,LON_var_name,typeof(LON_newgrd),getvardims(LON_newgrd));
; 1D variables
filevardef(o_file,PL_var_name,typeof(Pres_Levs),getvardims(Pres_Levs));
;
o_file->$U_var_name$ = UP_newgrd;
o_file->$V_var_name$ = VP_newgrd;
o_file->$T_var_name$ = TP_newgrd;
o_file->$Q_var_name$ = QP_newgrd;
o_file->$G_var_name$ = GP_newgrd;
o_file->$VOR_var_name$ = VOR_newgrd;
o_file->$DIV_var_name$ = DIV_newgrd;
o_file->$PS_var_name$ = PS_newgrd;
o_file->$ZS_var_name$ = ZS_newgrd;
o_file->$TS_var_name$ = TS_newgrd;
o_file->$LAT_var_name$ = LAT_newgrd;
o_file->$LON_var_name$ = LON_newgrd;
o_file->$PL_var_name$ = Pres_Levs;
;
delete(o_file)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("DONE!")




