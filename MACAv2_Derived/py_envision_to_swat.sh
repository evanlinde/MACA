#!/bin/bash
#
# This script is a wrapper for the python script point_subset.py which 
# processes a single model directory. The purpose of this script is to loop
# over the models and create output directories for the python script.
#
# This script requires the watershed name to be passed in as a command
# line argument.
#

watershed="$1"
centroids_file="${watershed}_centroids.txt"

# Skipping CCSM4 and NorESM1-M since these don't have relative humidity
MODELS=("bcc-csm1-1-m" "bcc-csm1-1" "BNU-ESM" "CanESM2" "CNRM-CM5" 
    "CSIRO-Mk3-6-0" "GFDL-ESM2G" "GFDL-ESM2M" "HadGEM2-CC365" "HadGEM2-ES365" 
    "inmcm4" "IPSL-CM5A-LR" "IPSL-CM5A-MR" "IPSL-CM5B-LR" "MIROC-ESM-CHEM" 
    "MIROC-ESM" "MIROC5" "MRI-CGCM3")

RCPS=("rcp45" "rcp85")

sourcepdir="/data/public/datasets/MACA/MACAv2_Derived/${watershed}_Envision"
destpdir="/data/public/datasets/MACA/MACAv2_Derived/${watershed}_SWAT"

prefixes=($(awk '{printf("%.3f_%.3f\n",$2,$1)}' ${centroids_file}))

# Set PATH variable so that we use the Anaconda python distribution
PATH=/opt/anaconda3/bin:$PATH

for model in ${MODELS[@]}; do
    sourcedir="${sourcepdir}/${model}"
    destdir="${destpdir}/${model}"
    mkdir -p ${destdir}
    python point_subset.py ${sourcedir} ${destdir} ${centroids_file}

    # The point_subset.py script creates files with a single variable
    # SWAT wants a two-column file with tasmax and tasmin separated by commas
    for p in ${prefixes[@]}; do
        for rcp in ${RCPS[@]}; do
            temperature_file="${destdir}/${p}_temperature_${rcp}.txt"
            paste -d, ${destdir}/${p}_tasm{ax,in}_${rcp}.txt > ${temperature_file}
            # Get rid of the duplicated header (remove the comma and everything
            # after it on the first line) and remove mid-line carriage returns 
            # on all lines.
            sed -i '1s/,.*//;s/\r,/,/g' ${temperature_file}
        done  # end rcp loop
    done  # end point prefix loop
done  # end model loop
