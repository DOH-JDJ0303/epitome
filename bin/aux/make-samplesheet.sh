#!/bin/bash

dir=$1

echo "taxa,assembly" > ref-samplesheet.csv

taxa=$(aws s3 ls $dir | grep -v 'pipeline_info' | sed 's/.*PRE//g' | tr -d ' \t/')
for t in $taxa
do
    dir2="${dir%/}/${t}/"
    segment=$(aws s3 ls $dir2 | sed 's/.*PRE//g' | tr -d ' \t/')
    for s in $segment
    do
        dir3="${dir2%/}/${s}/consensus/"
        consensus=$(aws s3 ls $dir3 | awk '{print $4}')
        for c in $consensus
        do
            file="${dir3%/}/${c}"
            bn=$(basename ${file%.fa})
            echo "${bn},${file}" >> ref-samplesheet.csv
        done      
    done
done

aws s3 cp ref-samplesheet.csv ${dir%/}/