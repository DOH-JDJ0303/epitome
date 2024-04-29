#!/bin/bash

# Sequence headers must be in this format '>EPI_ISL_212995|EPI715434|A/Indiana/18/2016|NP|5'

SEQS=$1

BN=$(basename ${SEQS%.*})

cat ${SEQS} \
    | sed 's/>.*$/@&@/g' \
    | tr -d '\n' \
    | tr '@' '\n' \
    | tail -n +2  \
    | paste - - \
    | tr '|' '\t' \
    | awk -v bn="${BN}" '{print $1"\n"$6 >> bn"_seg-"$5"_"$4".fa" }'