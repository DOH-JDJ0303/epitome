#!/bin/bash
version="1.0"

# consensus.sh
# Author: Jared Johnson, jared.johnson@doh.wa.gov

set -o pipefail

# get version info
if [ "$1" == "version" ]; then echo "${version}" && exit 0; fi

# input
name=$1
aln=$2
threads=$3

# format sequences to single lines, remove headers, and convert all bases to uppercase
cat ${aln} | sed 's/>.*$/@&@/g' | tr -d '\n' | tr '@' '\n' | grep -v '>' | tail -n +2 | awk '{print toupper($0)}' > seq.txt
# get the alignment length
seq_len=$(cat seq.txt | sed -n 1p | wc -c) && echo "Alignment Length: ${seq_len}"

# remove file from previous run and then select most common base at each site
rm site-list.txt 2>&1 /dev/null | true

# function for selecting the most common base - ties are chosen at random
select_base () {
    local site=$1
    cat seq.txt | sed 's/./& /g' | cut -f ${site} -d ' ' | sort | uniq -c | sed 's/^[ ]*//g' | tr ' ' ',' > ${site}-bc
    max_count=$(cat ${site}-bc | cut -f 1 -d ',' | sort -nr | sed -n 1p)
    most_common=$(cat ${site}-bc | tr ',' '\t' | awk -v max=${max_count} '$1 == max {print $2}' | shuf -n 1)
    line="${site}\t$(cat ${site}-bc | tr '\n' ';')\t${max_count}\t${most_common}"
    echo -e ${line} >> site-list.txt
    rm ${site}-bc
}
## export the function for use with xargs
export -f select_base
# run select_base in parallel
seq ${seq_len} | xargs -d \\n -n 1 -P ${threads} bash -c 'select_base "$@"' _

# save files
## header
echo ">${name}" > ${name}.fa
## sort sites and collapse to single line.
cat site-list.txt | sort -nk 1 | cut -f 4 | tr -d '\n -' | sed 's/$/\n/g' >> ${name}.fa
