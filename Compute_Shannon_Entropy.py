#####################################################################################
#
# Computes Shannon entropy at each [k,j,i] point for a given variable, for both the
# analysis (ratmanl) and corresponding background (atmf006, prior cycle). I am
# computing the Shannon entropy for both the analysis and background at the same
# time, because the PDF that the Shannon entropy is computed from should be produced
# using the same sampling range for both the analysis and background in order to
# do a fair comparison. So a single sampling range (x) is determined on each vertical
# level applied to all [j,i] points on that level for both the analysis and the
# background. Since the information content is computed as the difference in Shannon
# entropy between the two PDFs, the comparison is only fair if both PDFs use the same
# sampling range, hence why the analysis and background are handled together in this
# program.
#
# PDFs are estimated from the (20-member GDAS) ensemble using Kernel Density
# Estimation (KDE), which uses neighboring data in a histogram to produce an estimate
# of the true PDF underlying a sample histogram. This is achieved with fastKDE, which
# operates much more quickly than some other KDE approaches. This method is described
# in detail in the following papers:
#
# O’Brien, T. A., Kashinath, K., Cavanaugh, N. R., Collins, W. D. & O’Brien, J. P. A 
# fast and objective multidimensional kernel density estimation method: fastKDE. 
# Comput. Stat. Data Anal. 101, 148–160 (2016).
#
# O’Brien, T. A., Collins, W. D., Rauscher, S. A. & Ringler, T. D. Reducing the 
# computational cost of the ECF using a nuFFT: A fast and objective probability 
# density estimation method. Comput. Stat. Data Anal. 79, 222–234 (2014).
#
#####################################################################################
#
# Import modules
#
import numpy as np #................................................................. Array module
import matplotlib.pyplot as plt #.................................................... Array module
from netCDF4 import Dataset #........................................................ Array module
from fastkde import fastKDE #........................................................ Array module
import itertools #................................................................... Array module
#
#####################################################################################
#
# Read inputs
#
ratmanl_data_dir = input() #......................................................... Full path to directory containing ratmanl files
atmf006_data_dir = input() #......................................................... Full path to directory containing atmf006 files
ratmanl_prefix   = input() #......................................................... Prefix to ratmanl filenames
atmf006_prefix   = input() #......................................................... Prefix to atmf006 filenames (NOTE: Should be from a cycle 6 hrs prior to ratmanl files)
varname          = input() #......................................................... Variable name
target_dx        = input() #......................................................... Target granularity of range
nc_outfile_name  = input() #......................................................... Name of output netCDF file containing Shannon entropy
#ratmanl_data_dir='/data/users/bhoover/ATMS/GDAS_ensemble/native_to_prs_levs/'
#atmf006_data_dir='/data/users/bhoover/ATMS/GDAS_ensemble/native_to_prs_levs/'
#ratmanl_prefix='gdas.t12.ratmanl_plev.'
#atmf006_prefix='gdas.t06.atmf006_plev.'
#varname=U,V,T,Q,GHT,VOR,DIV,sp,sz,ts (lower-case are 2D)
#varname='U'
#target_dx='1.0'
#nc_outfile_name='gdas.shannon_entropy.ratmanl.t12.atmf006.t06.20200417.nc'
#
# Convert target_dx to a floating-point variable
#
target_dx = float(target_dx)
#
#####################################################################################
#
# Define a function to compute Shannon entropy:
# compute_H(v,H,it,k,x)
#
# INPUTS:
# v = Variable array, in (n_mem,n_lev,n_lat,n_lon) format
# H = Shannon entropy array, in (n_lev,n_lat,n_lon) format
# it = tuple containing (j,i) indices for (lon,lat)
# k = index for (lev)
# x = range to compute PDF over
#
# OUTPUTS:
# None, but the passed H array is modified so point (k,j,i) is
# assigned the value of Shannon entropy computed.
#
def compute_H(v,H,it,k,x):
    # Unpack it into j and i values
    j=it[0]
    i=it[1]
    # Create a vector of all members at one (k,j,i) point
    y=v[:,k,j,i].squeeze()
    if np.any(np.isnan(y)):
        H0=np.nan
    else:
        # Use fastKDE to create a pdf of y in the range defined by x
        y_pdf,x=fastKDE.pdf(y,axes=x)
        # Normalize y_pdf
        y_pdf = y_pdf/np.sum(y_pdf)
        # Initialize H
        H0=0.;
        # Loop through entire range on x, and sum contributions to
        # H0. Using base-2 log to keep entropy in units of bits
        for r in range(np.size(x)):
            if(y_pdf[r]>0.):
                H0 = H0 - y_pdf[r]*np.log2(y_pdf[r])
    # Assign H0 to H[k,j,i]
    H[k,j,i]=H0
    # Return
    return
#
#####################################################################################
#
# Read ensemble members into full (n_mem,n_lev,n_lat,n_lon) arrays
#
# First member: Establish (n_lev,n_lat,n_lon) arrays
mem='mem001'
# Background: atmf006
nc_file_bkg = atmf006_data_dir+atmf006_prefix+mem+'.nc' #............................ Background filename
nc_hdl_bkg  = Dataset(nc_file_bkg) #................................................. Background file-handle
bkg         = np.asarray(nc_hdl_bkg.variables[varname]).squeeze() #.................. Background data array
# Analysis: ratmanl
nc_file_ana = ratmanl_data_dir+ratmanl_prefix+mem+'.nc' #............................ Analysis filename
nc_hdl_ana  = Dataset(nc_file_ana) #................................................. Analysis file-handle
ana         = np.asarray(nc_hdl_ana.variables[varname]).squeeze() #.................. Analysis data array
# Grid data
lat = np.asarray(nc_hdl_bkg.variables['lat']).squeeze() #............................ Latitudes (1-D array, rectilinear grid)
lon = np.asarray(nc_hdl_bkg.variables['lon']).squeeze() #............................ Longitudes (1-D array, rectilinear grid)
lev = np.asarray(nc_hdl_bkg.variables['plev']).squeeze() #........................... Pressure levels (1-D array)
# Second member: Stack on first member along new axis=0
mem='mem002'
# Background: atmf006
nc_file_bkg = atmf006_data_dir+atmf006_prefix+mem+'.nc' #............................ Background filename
nc_hdl_bkg = Dataset(nc_file_bkg) #.................................................. Background file-handle
bkg = np.stack((bkg,np.asarray(nc_hdl_bkg.variables[varname]).squeeze()), #.......... Background data array
                axis=0)
# Analysis: ratmanl
nc_file_ana = ratmanl_data_dir+ratmanl_prefix+mem+'.nc' #............................ Analysis filename
nc_hdl_ana = Dataset(nc_file_ana) #.................................................. Analysis file-handle
ana = np.stack((ana,np.asarray(nc_hdl_ana.variables[varname]).squeeze()), #.......... Analysis data array
                axis=0)
# Additional members 3-20: concatenate along axis=0
mems=['mem003','mem004','mem005','mem006','mem007','mem008','mem009','mem010',
        'mem011','mem012','mem013','mem014','mem015','mem016','mem017','mem018',
        'mem019','mem020']
for mem in mems:
    # Background: atmf006
    nc_file_bkg = atmf006_data_dir+atmf006_prefix+mem+'.nc' #........................ Background filename
    nc_hdl_bkg = Dataset(nc_file_bkg) #.............................................. Background file-handle
    x=np.asarray(nc_hdl_bkg.variables[varname]).squeeze() #.......................... Background data array
    x=np.expand_dims(x,axis=0)
    bkg = np.concatenate((bkg,x),axis=0)
    # Analysis: ratmanl
    nc_file_ana = ratmanl_data_dir+ratmanl_prefix+mem+'.nc' #........................ Analysis filename
    nc_hdl_ana = Dataset(nc_file_ana) #.............................................. Analysis file-handle
    x=np.asarray(nc_hdl_ana.variables[varname]).squeeze() #.......................... Analysis data array
    x=np.expand_dims(x,axis=0)
    ana = np.concatenate((ana,x),axis=0)
# Values greater than 1.0e+10 are reassigned to np.nan
ana[np.where(ana>1.0e+10)]=np.nan
bkg[np.where(bkg>1.0e+10)]=np.nan
#
#####################################################################################
#
# Define dimension sizes: n_mem, n_lev, n_lat, n_lon
#
n_mem,n_lev,n_lat,n_lon=np.shape(bkg) #.............................................. Dimensions
#
# Create arrays to store Shannon entropy, H, of analysis and background
#
H_ana = np.zeros((n_lev,n_lat,n_lon)) #.............................................. Shannon entropy of analysis
H_bkg = np.zeros((n_lev,n_lat,n_lon)) #.............................................. Shannon entropy of background
#
# Loop over n_lev levels
#
for k in range(n_lev):
    print('Level:',k,'of',n_lev,'P=',lev[k])
    #
    # Define an appropriate range of values for x:
    #
    # 1) Define minimum and maximum values for variable
    # 2) Add 10% padding on either side to define sample range, xmin to xmax
    # 3) Choose a target granularity of range, target_dx
    # 4) Find appropriate power sample_pwr such that the range(xmin,xmax) across 
    #    2.**(sample_pwr)+1 points expresses a granularity closest to target_dx
    # 5) Define x based on xmin, xmax, and sample_pwr
    #
    xmin=min( #..................................................................... Minimum value on grid
             [
              np.floor(np.nanmin(bkg[:,k,:,:])) , 
              np.floor(np.nanmin(ana[:,k,:,:]))
             ]
            )
    xmax=max( #..................................................................... Maximum value on grid
             [
              np.ceil(np.nanmax(bkg[:,k,:,:])) , 
              np.ceil(np.nanmax(ana[:,k,:,:]))
             ]
            )
    # Add 10% padding to min/max range
    xmin = xmin - 0.1*abs(xmin)
    xmax = xmax + 0.1*abs(xmax)
    # Define a range of powers-of-2 to test
    pwr_range=range(20) #............................................................ Range of powers-of-2 to test
    # Initialize sample_pwr and sample_dx to absurd values
    sample_pwr = -9999. #............................................................ Current power-of-2
    sample_dx = 9999.*target_dx #.................................................... Current granularity of range
    # Loop over power-of-2 values
    for pwr in pwr_range:
        # Generate a sample range of size 2.**pwr + 1
        x=np.arange(xmin,xmax+0.01*target_dx,(xmax-xmin)/((2.**pwr))) #.............. Sample range
        # Define granularity of range 
        dx=x[1]-x[0] #............................................................... Granularity of range
        # If the granularity of range is closer to the target than the current
        # granularity of range, reassign sample_dx and sample_pwr to current
        # values
        if abs(dx-target_dx) < abs(sample_dx-target_dx):
            sample_dx = dx
            sample_pwr = pwr
    # Define range
    x=[list(np.arange(xmin,xmax+0.01*target_dx,(xmax-xmin)/((2.**sample_pwr))))] #... Range (as a list, required for fastkde)
    #
    # Loop through tuples of all (j,i) combinations and fill H_ana and H_bkg.
    # NOTE: Using itertools.product() is faster than nesting loops, which is
    #       why the loop is formulated this way
    #
    for it in itertools.product(range(n_lat),range(n_lon)):
        compute_H(bkg,H_bkg,it,k,x)
        compute_H(ana,H_ana,it,k,x)
#
#####################################################################################
#
# Output file
#
nc_out = Dataset( #.................................................................. Dataset object for output
                  nc_outfile_name  , # Dataset input: Output file name
                  "w"              , # Dataset input: Make file write-able
                  format="NETCDF4" , # Dataset input: Set output format to netCDF4
                )
# Dimensions
lat_dim  = nc_out.createDimension( #................................................... Output dimension
                               "latit" , # nc_out.createDimension input: Dimension name 
                               None    # nc_out.createDimension input: Dimension size limit ("None" == unlimited)
                             )
lon_dim  = nc_out.createDimension( #................................................... Output dimension
                               "longi" , # nc_out.createDimension input: Dimension name 
                               None    # nc_out.createDimension input: Dimension size limit ("None" == unlimited)
                             )
lev_dim = nc_out.createDimension( #.................................................... Output dimension
                               "plevl" , # nc_out.createDimension input: Dimension name 
                               None      # nc_out.createDimension input: Dimension size limit ("None" == unlimited)
                             )
# Variables
ratmanl_out = nc_out.createVariable( #............................................... Output variable
                                    "ratmanl_"+varname   , # nc_out.createVariable input: Variable name 
                                    "f8"      , # nc_out.createVariable input: Variable format 
                                    ( 
                                     "plevl" , # nc_out.createVariable input: Variable dimension
                                     "latit" , # nc_out.createVariable input: Variable dimension
                                     "longi"   # nc_out.createVariable input: Variable dimension
                                    )
                                   )
atmf006_out = nc_out.createVariable( #............................................... Output variable
                                    "atmf006_"+varname   , # nc_out.createVariable input: Variable name 
                                    "f8"      , # nc_out.createVariable input: Variable format 
                                    ( 
                                     "plevl" , # nc_out.createVariable input: Variable dimension
                                     "latit" , # nc_out.createVariable input: Variable dimension
                                     "longi"   # nc_out.createVariable input: Variable dimension
                                    )
                                   )
lat_out = nc_out.createVariable( #................................................................................ Output variable
                              "lat"   , # nc_out.createVariable input: Variable name 
                              "f8"    , # nc_out.createVariable input: Variable format 
                              ( 
                                "latit" # nc_out.createVariable input: Variable dimension
                              )
                            )
lon_out = nc_out.createVariable( #................................................................................ Output variable
                              "lon"   , # nc_out.createVariable input: Variable name 
                              "f8"    , # nc_out.createVariable input: Variable format 
                              ( 
                                "longi" # nc_out.createVariable input: Variable dimension
                              )
                            )
lev_out = nc_out.createVariable( #................................................................................ Output variable
                              "lev"   , # nc_out.createVariable input: Variable name 
                              "f8"    , # nc_out.createVariable input: Variable format 
                              ( 
                                "plevl" # nc_out.createVariable input: Variable dimension
                              )
                            )
ratmanl_out[:,:,:]=H_ana
atmf006_out[:,:,:]=H_bkg
lat_out[:]=lat
lon_out[:]=lon
lev_out[:]=lev
nc_out.close()
#
#####################################################################################
