#!/bin/bash
#
# Convert watershed-specific subsets from METDATA into a format usable
# in Envision. Envision-compatible netCDF files need to be in netCDF 3
# (classic) format and require dimension ordering like (time, lat, lon)
# instead of (time, lon, lat).
#
# Also convert relative humidity (for later conversion to SWAT) even though
# it's not needed for Envision.
#
# Additionally we will need to change some units and derive new variables.
#
# Convert temperature units from Kelvins to degrees Celsius.
#
# Calculate mean daily temperatures as:
#     tasmean = mean(tasmax, tasmin)
#
#
# This script requires the watershed name to be passed in as a command line
# argument.
#
# 
# Evan Linde, Oklahoma State University, 2017-10-09
#


if [ "$#" -ne 1 ]; then
    echo "This script requires the watershed directory name to be"
    echo "included as a command line argument."
    exit 1;
fi
watershed="$1" 

# directories
sourcedir="/data/public/datasets/MACA/METDATA_Derived/${watershed}"
destdir="/data/public/datasets/MACA/METDATA_Derived/${watershed}_Envision"
mkdir -p ${destdir}


# Create a temp directory
tmpdir=$(mktemp -d)

for year in {1979..2017}; do

    # Only correct dimensions for these variables
    for var in {pr,sph,srad,vs}; do
        infile="${sourcedir}/${var}_${year}.nc"
        outfile="${destdir}/${var}_${year}.nc"
        ncpdq -3 -a lat,lon ${infile} ${outfile}
    done

    # Rename precipitation variable so it will match the MACAv2 data
    # This isn't strictly necessary for Envision but makes it easier
    # to use both datasets.
    ncrename -v precipitation_amount,precipitation ${destdir}/pr_${year}.nc

    # Convert temperatures from Kelvins to degrees Celsius
    # and re-order dimensions
    for var in tmmn tmmx; do
        infile="${sourcedir}/${var}_${year}.nc"
        outfile="${destdir}/${var}_${year}.nc"
        tmpfile="${tmpdir}/${var}_${year}.nc"
        ncap2 -s 'air_temperature=(air_temperature - 273.15)' ${infile} ${tmpfile}
        ncatted -a units,air_temperature,o,char,'C' ${tmpfile}
        ncpdq -3 -a lat,lon ${tmpfile} ${outfile}
        rm ${tmpfile}
    done

    # Create the mean temperature file
    infiles=(${destdir}/tmm{x,n}_${year}.nc)  # array of two filenames
    outfile="${destdir}/tmmean_${year}.nc"
    nces ${infiles[@]} ${outfile}
    # Update metadata
    ncatted -a description,air_temperature,o,c,'Daily Mean Temperature' ${outfile}

    # Create mean relative humidity file
    infiles=(${sourcedir}/rm{in,ax}_${year}.nc)  # array of two filenames
    outfile="${destdir}/rmean_${year}.nc"
    tmpfile="${tmpdir}/rmean_${year}.nc"
    nces ${infiles[@]} ${tmpfile}
    ncpdq -3 -a lat,lon ${tmpfile} ${outfile}
    rm ${tmpfile}
    ncatted -a description,relative_humidity,o,c,'Daily Mean Relative Humidity' -a cell_methods,relative_humidity,o,c,'time: mean(interval: 24 hours)' ${outfile}

done  # end of years loop

# Cleanup: remove the temp directory we created at the beginning
rm -rf ${tmpdir}

