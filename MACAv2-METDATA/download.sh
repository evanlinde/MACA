#!/bin/bash
#
# Download netCDF files for MACAv2-METDATA
#

# Declare arrays for the models, year blocks, climate variables, and RCPs.
# These arrays will be used to build the URLs for the netCDF files we want
# to download.

MODELS=("bcc-csm1-1-m" "bcc-csm1-1" "BNU-ESM" "CanESM2" "CCSM4" "CNRM-CM5" 
    "CSIRO-Mk3-6-0" "GFDL-ESM2G" "GFDL-ESM2M" "HadGEM2-CC365" "HadGEM2-ES365" 
    "inmcm4" "IPSL-CM5A-LR" "IPSL-CM5A-MR" "IPSL-CM5B-LR" "MIROC-ESM-CHEM" 
    "MIROC-ESM" "MIROC5" "MRI-CGCM3" "NorESM1-M")

YEAR_BLOCKS=("2021_2025" "2026_2030" "2031_2035" "2036_2040" "2041_2045" 
    "2046_2050" "2051_2055" "2056_2060" "2061_2065" "2066_2070" "2071_2075" 
    "2076_2080" "2081_2085" "2086_2090" "2091_2095" "2096_2099")

VARS=("tasmax" "tasmin" "rhsmax" "rhsmin" "huss" "pr" "rsds" "uas" "vas")

RCPS=("rcp45" "rcp85")

# Beginning part of all the URLs we want
base_url="https://climate.northwestknowledge.net/MACAV2METDATA/MACAV2"

echo "Begin: $(date)"

for model in ${MODELS[@]}; do

    echo "Start downloading model ${model} $(date)"
    mkdir -p "${model}"

    # Go into each model directory to download files for that model. 
    # We can get back to our current directory later with the popd command.
    pushd "${model}"

    for wvar in ${VARS[@]}; do

        for block in ${YEAR_BLOCKS[@]}; do

            for rcp in ${RCPS[@]}; do

                # Build the URL we want to download using all the array 
                # variables that we're looping over.
                # First generate the filename
                filename="macav2metdata_${wvar}_${model}_r1i1p1_${rcp}_${block}_CONUS_daily.nc"
                # Then the full URL
                url="${base_url}/${model}/${filename}"

                # Start the download process backgrounded so we don't have to 
                # wait for it to finish before the next one starts.
                wget ${url} &

            done   # end rcp loop

            # Wait for backgrounded processes to finish so that we're not
            # trying to download ALL the files at once. We're doing this just
            # outside the rcp loop so that we only download two (length of 
            # RCPS array) files at a time.
            wait  

        done  # end year block loop

    done  # end wvar loop

    # Go back to the directory we were working in when we ran pushd earlier
    popd
    echo "Finish downloading model ${model} $(date)"

done  # end model loop

echo "End: $(date)"
