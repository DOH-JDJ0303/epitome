#!/bin/bash
set -o pipefail

# input
fasta=$1

# make each sequence a single line
cat ${fasta} | sed 's/>.*$/@&@/g' | tr -d '\n' | tr '@' '\n' | grep -v '>' | tail -n +2 | awk '{print toupper($0)}' > seqs

# remove sequences with ambiguous bases, and condense replicates
cat seqs | sort | uniq | grep -vE 'R|Y|M|K|S|W|H|B|V|D|N' > f1

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
## round 1
filter_outliers f1 f2
## round 2
filter_outliers f2 f3

# output summary of filtered samples
echo "total,filter1,filter2,filter3" > input-qc-summary.csv
echo "$(cat seqs | wc -l),$(cat f1 | wc -l),$(cat f2 | wc -l),$(cat f3 | wc -l)" >> input-qc-summary.csv
# output cleaned sequences & clean up
cat f3 | awk -v OFS='\n' '{print ">"NR, $1 > NR".fa"}'
rm seqs f1 f2 f3