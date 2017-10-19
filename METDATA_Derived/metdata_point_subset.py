#
# Extract data from netCDF files for use in SWAT.
#
# The netCDF files we're extracting from are the METDATA data with units
# converted for use with Envision. Since Envision and SWAT like the same units,
# we're not doing any unit conversions here like we would if extracting
# directly from the original METDATA netCDFs.
#
# We've also added netcdf files of daily mean relative humidity (which we're
# not using in Envision) into our Envision directories to make things easier 
# on our scripts.
#
# This script is intended to process a geographical subset of METDATA. For each
# climate variable, it reads all the matching files as a time series multi-file
# dataset, then loops over the set of points we want to use as SWAT weather 
# stations and writes an output file containing the data for that point.
#
# This script requires three command line arguments:
# 1. input directory of netCDF files
# 2. output directory where point subsets will be written (must already exist)
# 3. text file of latitude longitude points ("weather stations")
#
#
# Evan Linde, Oklahoma State University, 2017-10-11
#

from netCDF4 import MFDataset
import numpy as np
import sys
import os

# Get inputs and output locations from command line arguments
in_dir = sys.argv[1]
out_dir = sys.argv[2]    # This directory should already exist.
point_file = sys.argv[3] # Each line should have lon and lat separated 
                         # by a space, e.g. "-97.07 36.122".

# Dictionary of climate variables
# Each key is the variable name in the netcdf file name
# Each value is the variable name *inside* the netcdf file
clim_vars = {'pr':'precipitation',
             'rmean':'relative_humidity',
             'srad':'surface_downwelling_shortwave_flux_in_air',
             'tmmx':'air_temperature',
             'tmmn':'air_temperature',
             'vs':'wind_speed'}

# Read the point file into a numpy array
points = np.loadtxt(point_file)


for fvar in clim_vars.keys():
    ncvar = clim_vars[fvar]

    # file names will have the form of fvar_year.nc
    infilename_glob = fvar+'_*.nc' 

    # combine the input directory with filename glob for a full path glob
    infile_glob = os.path.join(in_dir, infilename_glob)

    # read multi-file dataset (all years for our current variable)
    mfd = MFDataset(infile_glob, aggdim='day')

    for p in points:
        # each point p is itself an array in the form of [lon, lat]
        # assign longitude and latitude to variables
        lon,lat = p[0],p[1]

        # base name for output file
        outfilename = '%.3f_%.3f_%s.txt' % (lat, lon, fvar)

        # full path for output file
        outfile = os.path.join(out_dir, outfilename)

        # find the lat and lon indices in the netcdf dataset
        # that are closest to our point
        latindex = np.where(abs(lat - mfd.variables["lat"][:]) < 0.020833)[0][0]
        lonindex = np.where(abs(lon - mfd.variables["lon"][:]) < 0.020833)[0][0]

        # the subset of the netcdf dataset that matches our point
        outdata = mfd.variables[ncvar][:,latindex,lonindex]

        # save the point data to a file
        np.savetxt(outfile, outdata, fmt='%.3f', newline='\r\n', header='19790101', comments='')
            

