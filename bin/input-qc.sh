#!/bin/bash
version="1.0"

# input-qc.sh
# Author: Jared Johnson, jared.johnson@doh.wa.gov

set -o pipefail

# get version info
if [ "$1" == "version" ]; then echo "${version}" && exit 0; fi

# input
fasta=$1
prefix=$2
expected_length=$3
length_threshold=$4
max_cluster=$5

# make each sequence a single line
cat ${fasta} | sed 's/>.*$/@&@/g' | tr -d '\n' | tr '@' '\n' | grep -v '>' | tail -n +2 | awk '{print toupper($0)}' > seqs

#---- FILTER 1: CONDENSE REPLICATES ----#
cat seqs | sort | uniq > f1

#---- FILTER 2: REMOVE SEQUENCES WITH AMBIGUOUS BASES ----#
cat f1 | grep -E '^[ATCG]+$' > f2
if [ ! -s f2 ]
then
    echo "Error: All sequences had ambiguous bases!" && exit 1
fi

#---- FILTER 3: REMOVE SEQUENCES DIFFERING BY MORE THAN A PERCENTAGE OF THE EXPECTED LENGTH ----#
cat f2 | awk -v exp_len=${expected_length} -v len_thresh=${length_threshold} 'length($1) > exp_len*(1-len_thresh) && length($1) < exp_len*(1+len_thresh) {print $0}' > f3
if [ ! -s f3 ]
then
    echo "Error: All sequences differ from the expected length: ${expected_length}!" && exit 1
fi

#---- FILTER 4: REMOVE OUTLIERS BASED ON LENGTH AND GC CONTENT ----#
# function for filtering outliers based on sequence length and GC content
filter_outliers () {
    # input
    local seqs=$1
    local out=$2

    # calculate sequence length & GC
    cat $seqs | awk '{print $1, length($1)}' |  tr -d 'AT' | awk -v OFS='\t' '{print $2, 100*length($1)/$2}' > length_gc

    # calculate the mean and stdev of length and GC
    mean_len=$(cat length_gc | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
    sd_len=$(cat length_gc | awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR); }')
    mean_gc=$(cat length_gc | awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }')
    sd_gc=$(cat length_gc | awk '{delta = $2 - avg; avg += delta / NR; mean2 += delta * ($2 - avg); } END { print sqrt(mean2 / NR); }')

    # filter sequences based on length and GC content
    paste seqs length_gc | awk -v ml=${mean_len} -v sl=${sd_len} -v mg=${mean_gc} -v sg=${sd_gc} '$2 >= ml-(3*sl) && $2 <= ml+(3*sl) && $3 >= mg-(3*sg) && $3 <= mg+(3*sg) {print $0}' > ${out}

    # clean up
    rm length_gc
}

# filter outliers
filter_outliers f3 f4
if [ ! -s f4 ]
then
    echo "Error: All sequences were considered outliers. This should not happen. Check your sequences and parameters." && exit 1
fi

# output summary of filtered samples
echo "total,filter1,filter2,filter3,filter4" > ${prefix}-qc-summary.csv
echo "$(cat seqs | wc -l),$(cat f1 | wc -l),$(cat f2 | wc -l),$(cat f3 | wc -l),$(cat f4 | wc -l)" >> ${prefix}-qc-summary.csv

# output sequences
## shuffle and number sequences
cat f4 | shuf | awk -v OFS='\t' '{print ">"NR, $1}' > shufd
## get seq count
n_clean=$(cat shufd | wc -l) 
## parition sequences
cat shufd | sed -n "1,${max_cluster}p" | tr '\t' '\n' > ${prefix}.top.fa
cat shufd | sed -n "$((max_cluster+1)),\$p" | tr '\t' '\n' > ${prefix}.remainder.fa
# clean up
rm seqs f1 f2 f3 f4 shufd