#!/bin/bash
#
# This script is a wrapper for the python script metdata_point_subset.py 
# which processes a directory containing a geographical subset of the METDATA
# dataset.
#
# This script performs some post-processing to create a SWAT temperature
# file with max and min temperatures.
#
# This script requires the watershed name to be passed in as a command
# line argument.
#

watershed="$1"
centroids_file="${watershed}_centroids.txt"

sourcedir="/data/public/datasets/MACA/METDATA_Derived/${watershed}_Envision"
destdir="/data/public/datasets/MACA/METDATA_Derived/${watershed}_SWAT"
mkdir -p ${destdir}

prefixes=($(awk '{printf("%.3f_%.3f\n",$2,$1)}' ${centroids_file}))

# Set PATH variable so that we use the Anaconda python distribution
PATH=/opt/anaconda3/bin:$PATH

python metdata_point_subset.py ${sourcedir} ${destdir} ${centroids_file}

# The point_subset.py script creates files with a single variable
# SWAT wants a two-column file with tasmax and tasmin separated by commas
for p in ${prefixes[@]}; do
    temperature_file="${destdir}/${p}_temperature.txt"
    paste -d, ${destdir}/${p}_tmm{x,n}.txt > ${temperature_file}
    # Get rid of the duplicated header (remove the comma and everything
    # after it on the first line) and remove mid-line carriage returns 
    # on all lines.
    sed -i '1s/,.*//;s/\r,/,/g' ${temperature_file}
done  # end point prefix loop
