#!/bin/bash

for v in {pr,sph,srad,tmmx,tmmn,vs}; do
    for y in {1979..2017}; do
        url="http://www.northwestknowledge.net/metdata/data/${v}_${y}.nc"
        wget ${url}
    done
done
