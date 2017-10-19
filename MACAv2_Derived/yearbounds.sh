#!/bin/bash
#
# Print the first (Jan 1) and last (Dec 31) days of a given year as the
# number of days since a reference date.
#
# The output will be two floating point numbers rounded to one decimal point
# and separated with a comma. (This is intended to be used in a command
# to create a time subset of a netCDF file.)
#
# This script requires the target year to be passed in as a command
# line argument.
#

# the target year
y="$1"

# Calculate the reference time in Unix time
reftime=$(date --date='1900-01-01' +%s)

# Calculate January 1st and December 31st of the given year in Unix time
s1=$(date --date="${y}-01-01" +%s)
s2=$(date --date="${y}-12-31" +%s)

# Calculate difference from reference time and convert from seconds to days
d1=$(echo "(${s1} - ${reftime})/86400" | bc) 
d2=$(echo "(${s2} - ${reftime})/86400" | bc)

# Print the numbers in floating point format
printf "%.1f,%.1f\n" ${d1} ${d2}
