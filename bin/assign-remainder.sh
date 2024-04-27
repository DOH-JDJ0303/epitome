#!/bin/bash

top=$1
remainder=$2
clusters=$3
threshold=$4
threads=$5

[ -f reps.txt ] && rm reps.txt
[ -f reps.fa ] && rm reps.fa
for c in $(cat $clusters | cut -f 4 -d ',' | grep -v 'cluster' | sort | uniq)
do
    rep=$(cat $clusters | tr ',' '\t' | awk -v OFS=',' -v c=$c '$4 == c {print $1,$4}' | sed -n 1p)
    echo ${rep} >> reps.csv
    cat ${top} | paste - - | awk -v s=">$(echo ${rep} | cut -f 1 -d ',')" '$1 == s {print}' | tr '\t' '\n' >> reps.fa 
done

mash sketch -p $threads -o reps -i reps.fa
mash sketch -p $threads -o remainder -i ${remainder}

mash dist -p $threads reps.msh remainder.msh | awk -v OFS=',' -v t=${threshold} '$3 < t {print $1,$2,$3}' > remainder-mash.csv
cat ${remainder} | grep '>' | tr -d '>' > remainder-list.csv