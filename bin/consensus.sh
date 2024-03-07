#!/bin/bash

# input
name=$1
aln=$2

# function for selecting the most common base - ties are chosen at random
select_base () {
    local site=$1
    cat seq.txt | sed 's/./& /g' | cut -f ${site} -d ' ' | sort | uniq -c | sed 's/^[ ]*//g' | tr ' ' ',' > base-counts.txt
    max_count=$(cat base-counts.txt | cut -f 1 -d ',' | sort -nr | sed -n 1p)
    most_common=$(cat base-counts.txt | tr ',' '\t' | awk -v max=${max_count} '$1 == max {print $2}' | shuf -n 1)
    line="${site}\t$(cat base-counts.txt | tr '\n' ';')\t${max_count}\t${most_common}"
    echo -e ${line} | tee -a site-list.txt
    rm base-counts.txt
}
## export the function for use with parallel
export -f select_base

# format sequences to single lines, remove headers, and convert all bases to uppercase
cat ${aln} | sed 's/>.*$/@&@/g' | tr -d '\n' | tr '@' '\n' | grep -v '>' | tail -n +2 | awk '{print toupper($0)}' > seq.txt
# get the alignment length
seq_len=$(cat seq.txt | sed -n 1p | wc -c) && echo "Alignment Length: ${seq_len}"

# remove file from previous run and then select most common base at each site
rm site-list.txt 2>&1 /dev/null | true
#seq 1 ${seq_len} | parallel -j0 select_base

for i in $(seq ${seq_len})
do
    select_base ${i}
done

# save files
## header
echo ">${name}" > ${name}.fa
## sort sites and collapse to single line.
cat site-list.txt | sort -nk 1 | cut -f 4 | tr -d '\n -' | sed 's/$/\n/g' >> ${name}.fa
## save sequence length 
echo ${seq_len} > length.csv
