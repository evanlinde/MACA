#
# Extract data from netCDF files for use in SWAT.
#
# The netCDF files we're extracting from are the MACAv2 data with units
# converted for use with Envision and longitudes in the Western hemisphere
# converted to negative values. Since Envision and SWAT like the same units,
# we're not doing any unit conversions here like we would if extracting
# directly from the original MACAv2 netCDFs.
#
# This script is intended to process an entire model directory. For each
# combination of climate variable and RCP, it reads all the matching files
# as a time series multi-file dataset, then loops over the set of points
# we want to use as SWAT weather stations and writes an output file containing
# the data for that point.
#
# This script requires three command line arguments:
# 1. input directory of netCDF files
# 2. output directory where point subsets will be written (must already exist)
# 3. text file of latitude longitude points ("weather stations")
#
#
# Evan Linde, Oklahoma State University, 2017-10-10
#

from netCDF4 import MFDataset
import numpy as np
import sys
import os

# Get inputs and output locations from command line arguments
in_dir = sys.argv[1]
out_dir = sys.argv[2]    # This directory should already exist.
point_file = sys.argv[3] # Each line should have lon and lat separated 
                         # by a space, e.g. "-97.12 36.004".


rcps = ["rcp45", "rcp85"]

# Dictionary of climate variables
# Each key is the variable name in the netcdf file name
# Each value is the variable name *inside* the netcdf file
clim_vars = {'pr':'precipitation', 
             'rhsmax':'relative_humidity',
             'rhsmin':'relative_humidity',
             'rsds':'surface_downwelling_shortwave_flux_in_air',
             'tasmax':'air_temperature',
             'tasmin':'air_temperature', 
             'uas':'eastward_wind',
             'vas':'northward_wind'}

# Read the point file into a numpy array
points = np.loadtxt(point_file)

# expression for base filename of input file containing both
# string formatters (e.g. '%s') and wildcards for globbing
infilename_expr = 'macav2metdata_%s_*_r?i1p1_%s_*.nc'


def get_point_data(file_glob, ncvar, lat, lon):
    # read multi-file dataset
    mfd = MFDataset(file_glob, aggdim='time')

    # find the lat and lon indices in the netcdf dataset
    # that are closest to our point
    latindex = np.where(abs(lat - mfd.variables["lat"][:]) < 0.020833)[0][0]
    lonindex = np.where(abs(lon - mfd.variables["lon"][:]) < 0.020833)[0][0]

    # the subset of the netcdf dataset that matches our point
    return mfd.variables[ncvar][:,latindex,lonindex]



for rcp in rcps:
    for fvar in clim_vars.keys():
        ncvar = clim_vars[fvar][0]
        varfmt = clim_vars[fvar][1]


        # replace the string formatters with values
        infilename_glob = infilename_expr % (fvar, rcp)

        # combine the input directory with filename glob for a full path glob
        infile_glob = os.path.join(in_dir, infilename_glob)

        for p in points:
            # each point p is itself an array in the form of [lon, lat]
            # assign longitude and latitude to variables
            lon,lat = p[0],p[1]

            # base name for output file
            outfilename = '%.3f_%.3f_%s_%s.txt' % (lat, lon, fvar, rcp)

            # full path for output file
            outfile = os.path.join(out_dir, outfilename)

            outdata = get_point_data(infile_glob, ncvar, lat, lon)

            # save the point data to a file
            np.savetxt(outfile, outdata, fmt=varfmt)
            

