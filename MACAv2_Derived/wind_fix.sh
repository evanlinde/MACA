#!/bin/bash


if [ "$#" -ne 1 ]; then
    echo "This script requires the watershed directory name to be included as a command line argument."
    exit 1;
fi
watershed="$1" 

# Parent directories
sourcepdir="/data/public/datasets/MACAv2_Derived/${watershed}"
destpdir="/data/public/datasets/MACAv2_Derived/${watershed}_Envision"


# Declare arrays for the variables we'll need to loop over
MODELS=("bcc-csm1-1-m" "bcc-csm1-1" "BNU-ESM" "CanESM2" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "GFDL-ESM2G" "GFDL-ESM2M" "HadGEM2-CC365" "HadGEM2-ES365" "inmcm4" "IPSL-CM5A-LR" "IPSL-CM5A-MR" "IPSL-CM5B-LR" "MIROC-ESM-CHEM" "MIROC-ESM" "MIROC5" "MRI-CGCM3" "NorESM1-M")
YEAR_BLOCKS=("2021_2025" "2026_2030" "2031_2035" "2036_2040" "2041_2045" "2046_2050" "2051_2055" "2056_2060" "2061_2065" "2066_2070" "2071_2075" "2076_2080" "2081_2085" "2086_2090" "2091_2095" "2096_2099")
RCPS=("rcp45" "rcp85")


# Create a temp directory
tmpdir=$(mktemp -d)

for model in ${MODELS[@]}; do

    for rcp in ${RCPS[@]}; do

        for block in ${YEAR_BLOCKS[@]}; do

            bfirst=${block:0:4}  # first year in year block
            blast=${block:5:4}   # last year in year block
            sourcedir="${sourcepdir}/${model}"  # subset files for current model
            destdir="${destpdir}/${model}"      # output files for current model
            mkdir -p ${destdir}/yearly     # output dir for yearly files

            # Arbitrarily using the huss file as a template for other 
            # filenames. Most models have "r1i1p1" but at least one has 
            # "r6i1p1", so instead of a condition based on the model name, 
            # we're detecting the correct string with bash's glob expressions 
            # (i.e. the "r?i1p1" part). This could be a problem if more than 
            # one file matches the expression.
            huss_basename=$(basename $(echo ${sourcedir}/macav2metdata_huss_${model}_r?i1p1_${rcp}_${block}_${watershed}_daily.nc))
            # Prefix for yearly files; we'll add "${year}.nc" to the end 
            # when creating yearly files
            huss_y_prefix="${huss_basename/_${block}_${watershed}_daily.nc/_${watershed}_daily}"

            # Set variables for the input and output files using the huss 
            # filename as a template replacing "huss" with the target 
            # variable name
            subset_vas="${sourcedir}/${huss_basename/huss/vas}"
            subset_uas="${sourcedir}/${huss_basename/huss/uas}"
            envision_wind="${destdir}/${huss_basename/huss/wind}"
            rm ${envision_wind}

            # Do all the conversions
            #    subset vas and uas --> envision wind

            # Calculated variable: wind (wind_speed)
            # Make a temporary copy of the vas file. (It will be modified.)
            cp ${subset_vas} ${tmpdir}/vas.nc
            # Append the uas variable (eastward_wind) into the vas file
            ncks -A ${subset_uas} ${tmpdir}/vas.nc
            # Calculate wind speed. The output file from this command will have
            # three variables: northward_wind, eastward_wind, and wind_speed.
            ncap2 -s 'lon=lon-360; wind_speed=sqrt(northward_wind^2 + eastward_wind^2)' ${tmpdir}/vas.nc ${tmpdir}/wind.nc
            # Extract the wind_speed variable into its own file
            ncks -3 -v wind_speed ${tmpdir}/wind.nc ${envision_wind}
            # Get rid of temporary files
            rm ${tmpdir}/vas.nc ${tmpdir}/wind.nc
            # Update metadata. (The wind_speed variable inherited the
            # metadata from the vas variable, northward_wind.)
            ncatted -a comments,wind_speed,o,c,'Surface (10m) wind speed' -a long_name,wind_speed,o,c,'Wind Speed' -a standard_name,wind_speed,o,c,'wind_speed' ${envision_wind}

            # Split all the envision files into yearly files
            # Doing a C-style for loop here because "{${bfirst}..${blast}}"
            # doesn't work like it would with literals.
            for ((y=${bfirst}; y<=${blast}; y++)); do
                infile="${destdir}/${huss_basename/huss/wind}"
                outfile="${destdir}/yearly/${huss_y_prefix/huss/wind}_${y}.nc"
                rm ${outfile}
                # Year Bounds
                # Get the first and last day of the given year in the form
                # of "first,last" where first and last are floating point
                # numbers representing "days since 1900-01-01".
                yb=$(bash yearbounds.sh ${y})

                ncks -d time,${yb} ${infile} ${outfile}

            done  # end of years inside year block loop

        done  # end of year blocks loop

    done  # end of rcps loop

done  # end of models loop

# Cleanup: remove the temp directory we created at the beginning
rm -rf ${tmpdir}

