#!/bin/bash
version="1.1"

# input-qc.sh
# Author: Jared Johnson, jared.johnson@doh.wa.gov

set -o pipefail

# get version info
if [ "$1" == "version" ]; then echo "${version}" && exit 0; fi

# input
FASTA=$1
PREFIX=$2
LENGTH=$3
LENGTH_THRESHOLD=$4
AMB_THRESHOLD=$5
MAX_CLUSTER=$6

# make each sequence a single line
cat ${FASTA} | \
    sed 's/>.*$/@&@/g' | \
    tr -d '\n' | \
    tr '@' '\n' | \
    grep -v '>' | \
    tail -n +2 | \
    awk '{print toupper($0)}' \
    > seqs

# check for unexpected characters
bad_chars=$(cat seqs | tr -d '\-ATCGRYSWKMBDHVN\n')
if [[ $(echo "${bad_chars}" | wc -c) > 1 ]]
then
    echo "Error: the input contains illegal characters: ${bad_chars}"
    exit 1
fi 

#---- FILTER 1: CONDENSE REPLICATES ----#
cat seqs | sort | uniq > f1

#---- FILTER 2: REMOVE SEQUENCES WITH AMBIGUOUS BASES ----#
paste f1 <(cat f1 | sed -E 's/[^N-]//g') | \
    awk -v amb=${AMB_THRESHOLD} '{if(length($2)/length($1) <= amb){print $1 >> "f2"} else {print $1 >> "amb_seqs"}}'
if [ ! -s f2 ]
then
    echo "Error: All sequences had ambiguous bases!" && exit 1
fi

#---- FILTER 3: REMOVE SEQUENCES DIFFERING BY MORE THAN A PERCENTAGE OF THE EXPECTED LENGTH ----#
cat f2 | \
    awk -v exp_len=${LENGTH} -v len_thresh=${LENGTH_THRESHOLD} 'length($1) > exp_len*(1-len_thresh) && length($1) < exp_len*(1+len_thresh) {print $0}' \
    > f3
if [ ! -s f3 ]
then
    echo "Error: All sequences differ from the expected length: ${LENGTH}!" && exit 1
fi

#---- FILTER 4: REMOVE OUTLIERS BASED ON LENGTH ----#
# function for filtering outliers based on sequence length
filter_outliers () {
    # input
    local seqs=$1
    local out=$2

    # calculate sequence length & GC
    cat $seqs | \
        awk '{print $1, length($1)}' |  \
        tr -d 'AT' | \
        awk -v OFS='\t' '{print $2, 100*length($1)/$2}' \
        > m_len

    # calculate the mean and stdev of length
    mean_len=$(cat m_len | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
    sd_len=$(cat m_len | awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR); }')

    # filter outliers
    paste $seqs m_len | \
        awk -v ml=${mean_len} -v sl=${sd_len} '$2 >= ml-(3*sl) || $2 <= ml+(3*sl) {print $0}' \
        > ${out}

    # clean up
    rm m_len
}

# filter outliers
filter_outliers f3 f4
if [ ! -s f4 ]
then
    echo "Error: All sequences were considered outliers. This should not happen. Check your sequences and parameters." && exit 1
fi

# output summary of filtered samples
echo "total,filter1,filter2,filter3,filter4" > ${PREFIX}-qc-summary.csv
echo "$(cat seqs | wc -l),$(cat f1 | wc -l),$(cat f2 | wc -l),$(cat f3 | wc -l),$(cat f4 | wc -l)" >> ${PREFIX}-qc-summary.csv

# output sequences
## shuffle and number sequences
cat f4 | shuf | awk -v OFS='\t' '{print ">"NR, $1}' > shufd
## get seq count
n_clean=$(cat shufd | wc -l) 
## parition sequences - this is necessary for large datasets
cat shufd | sed -n "1,${MAX_CLUSTER}p" | tr '\t' '\n' > ${PREFIX}.top.fa
cat shufd | sed -n "$((MAX_CLUSTER+1)),\$p" | tr '\t' '\n' > ${PREFIX}.remainder.fa
# clean up
rm seqs f1 f2 f3 f4 shufd
