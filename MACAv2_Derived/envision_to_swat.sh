#!/bin/bash

watershed="$1"
centroids_file="${watershed}_centroids.txt"
lons_file="${watershed}_lons.txt"
lats_file="${watershed}_lats.txt"
var_map="vars.txt"

MODELS=("bcc-csm1-1-m" "bcc-csm1-1" "BNU-ESM" "CanESM2" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "GFDL-ESM2G" "GFDL-ESM2M" "HadGEM2-CC365" "HadGEM2-ES365" "inmcm4" "IPSL-CM5A-LR" "IPSL-CM5A-MR" "IPSL-CM5B-LR" "MIROC-ESM-CHEM" "MIROC-ESM" "MIROC5" "MRI-CGCM3" "NorESM1-M")
YEAR_BLOCKS=("2021_2025" "2026_2030" "2031_2035" "2036_2040" "2041_2045" "2046_2050" "2051_2055" "2056_2060" "2061_2065" "2066_2070" "2071_2075" "2076_2080" "2081_2085" "2086_2090" "2091_2095" "2096_2099")
RCPS=("rcp45" "rcp85")

sourcepdir="/data/public/datasets/MACA/MACAv2_Derived/${watershed}_Envision"
destpdir="/data/public/datasets/MACA/MACAv2_Derived/${watershed}_SWAT"


function closest_point {
    point=$1
    file=$2
    #echo "point: ${point}, file: ${file}" 1>&2
    awk -v point=${point} 'function abs(v) {return v < 0 ? -v : v} abs($1 - point) < 0.02083 {closest=$1} END{print closest}' ${file}
}


tmpdir=$(mktemp -d)

while read lon lat; do
    #blah
    closest_lon=$(closest_point ${lon} ${lons_file})
    closest_lat=$(closest_point ${lat} ${lats_file})
    lonrounded=$(printf "%.3f\n" ${closest_lon})
    latrounded=$(printf "%.3f\n" ${closest_lat})
    eb=$(echo "${lonrounded} + 0.001" | bc)  # east boundary
    wb=$(echo "${lonrounded} - 0.001" | bc)  # west boundary
    nb=$(echo "${latrounded} + 0.001" | bc)  # north boundary
    sb=$(echo "${latrounded} - 0.001" | bc)  # south boundary
    echo "${lon},${lat} --> ${lonrounded},${latrounded}"
    for model in ${MODELS[@]}; do
        sourcedir="${sourcepdir}/${model}"
        destdir="${destpdir}/${model}"
        mkdir -p ${destdir}
        for rcp in ${RCPS[@]}; do
            while read fvar ncvar; do
                outfile="${destdir}/${latrounded}_${lonrounded}_${fvar}_${rcp}.txt"
                echo "20210101" > ${outfile}
                for block in ${YEAR_BLOCKS[@]}; do
                    #infile=$(echo ${sourcedir}/macav2metdata_${fvar}_${model}_r?i1p1_${rcp}_${watershed}_daily_${year}.nc)
                    infile=$(echo ${sourcedir}/macav2metdata_${fvar}_${model}_r?i1p1_${rcp}_${block}_${watershed}_daily.nc)
                    temp_netcdf="${tmpdir}/temp.nc"
                    ncks -d lat,${sb},${nb} -d lon,${wb},${eb} ${infile} ${temp_netcdf} || { echo ${infile}; exit; }
                    ncdump -v ${ncvar} ${temp_netcdf} | sed -n '/^ '"${ncvar}"' =/,/\;/{s/ *'"${ncvar}"' = *//g;p}' | tr -d '\n ' | tr ',;' '\n\n' >> ${outfile}
                    rm ${temp_netcdf}
                done  # end year blocks loop
            done < ${var_map}  # end variables loop
        done  # end rcp loop
    done  # end model loop
done < ${centroids_file}  # end centroids loop

rm -rf ${tmpdir}
