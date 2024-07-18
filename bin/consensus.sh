#!/bin/bash
version="2.0"

# consensus.sh
# Author: Jared Johnson, jared.johnson@doh.wa.gov

set -o pipefail

# get version info
if [ "$1" == "version" ]; then echo "${version}" && exit 0; fi

# input
name=$1
aln=$2

# format sequences to single lines, remove headers, and convert all bases to uppercase
cat ${aln} | sed 's/>.*$/@&@/g' | tr -d '\n' | tr '@' '\n' | grep -v '>' | tail -n +2 | awk '{print toupper($0)}' > seq.txt

# check for unexpected characters
bad_chars=$(cat seq.txt | tr -d '\-ATCG\n' | wc -c)
if [[ ${bad_chars} > 0 ]]
then
    echo "Error: the alignment contains characters other than '-', 'A', 'T', 'C', 'G'."
    exit 1
fi 

# get the alignment length & sequence count
seq_len=$(cat seq.txt | sed -n 1p | tr -d '\t\r\n ' | wc -c) && echo "Alignment Length: ${seq_len}"
seq_count=$(cat seq.txt | wc -l) && echo "Number of Sequences: ${seq_count}"

# count the occurance of each base and report the most supported call - choose ties at random
echo "site,-,A,T,C,G,max,call" > site-list.txt
cat seq.txt \
    | sed 's/./&\t/g' \
    | awk -v OFS=',' -v rc=${seq_count} '{for (i=1; i<=NF; i++){if(NR == 1){nd[i] = 0; na[i] = 0; nt[i] = 0; nc[i] = 0; ng[i] = 0}; if($i == "-"){nd[i]++}; if($i == "A"){na[i]++}; if($i == "T"){nt[i]++}; if($i == "C"){nc[i]++}; if($i == "G"){ng[i]++}; if(NR == rc){print i,nd[i],na[i],nt[i],nc[i],ng[i]}}}' \
    | awk -v FS=',' -v OFS=',' '{max=$2; if(max < $3){max=$3} else if(max < $4){max=$4} else if(max < $5){max=$5} else if(max < $6){max=$6}; print $0,max }' \
    | awk -v FS=',' -v OFS=',' 'BEGIN{srand()}{delete call; if ($2 == $7){call[length(call)] = "-"}; if ($3 == $7){call[length(call)] = "A"}; if ($4 == $7){call[length(call)] = "T"}; if ($5 == $7){call[length(call)] = "C"}; if ($6 == $7){call[length(call)] = "G"}; r=int(rand()*length(call)); print $0,call[r]}' \
    >> site-list.txt

# save files
echo ">${name}" > ${name}.fa
cat site-list.txt | tail -n +2 | cut -f 8 -d ',' | tr -d '\n -' | sed 's/$/\n/g' >> ${name}.fa