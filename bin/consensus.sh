#!/bin/bash
version="2.0"

# consensus.sh
# Author: Jared Johnson, jared.johnson@doh.wa.gov

set -o pipefail

# get version info
if [ "$1" == "version" ]; then echo "${version}" && exit 0; fi

# input
NAME=$1
ALN=$2

# format sequences to single lines, remove headers, and convert all bases to uppercase
cat ${ALN} | sed 's/>.*$/@&@/g' | tr -d '\n' | tr '@' '\n' | grep -v '>' | tail -n +2 | awk '{print toupper($0)}' > seq.txt

# check for unexpected characters
bad_chars=$(cat seq.txt | tr -d '\-ATCGRYSWKMBDHVN\n')
if [[ $(echo "${bad_chars}" | wc -c) > 1 ]]
then
    echo "Error: the alignment contains illegal characters: ${bad_chars}"
    exit 1
fi 

# get the alignment length & sequence count
seq_len=$(cat seq.txt | sed -n 1p | tr -d '\t\r\n ' | wc -c) && echo "Alignment Length: ${seq_len}"
seq_count=$(cat seq.txt | wc -l) && echo "Number of Sequences: ${seq_count}"

# count the occurance of each base and report the most supported call - choose ties at random
echo "site,-,A,T,C,G,R,Y,S,W,K,M,B,D,H,V,N,max,call" > site-list.txt
cat seq.txt \
    | sed 's/./&\t/g' \
    | awk -v OFS=',' -v rc=${seq_count} '{for (i=1; i<=NF; i++){
                                                                if(NR == 1){nx[i] = 0; 
                                                                            na[i] = 0; 
                                                                            nt[i] = 0; 
                                                                            nc[i] = 0; 
                                                                            ng[i] = 0;
                                                                            nr[i] = 0;
                                                                            ny[i] = 0;
                                                                            ns[i] = 0;
                                                                            nw[i] = 0;
                                                                            nk[i] = 0;
                                                                            nm[i] = 0
                                                                            nb[i] = 0;
                                                                            nd[i] = 0;
                                                                            nh[i] = 0;
                                                                            nv[i] = 0;
                                                                            nn[i] = 0}; 
                                                                if($i == "-"){nx[i]++}; 
                                                                if($i == "A"){na[i]++}; 
                                                                if($i == "T"){nt[i]++}; 
                                                                if($i == "C"){nc[i]++}; 
                                                                if($i == "G"){ng[i]++};
                                                                if($i == "R"){nr[i]++; na[i]++; ng[i]++};
                                                                if($i == "Y"){ny[i]++; nc[i]++; nt[i]++};
                                                                if($i == "S"){ns[i]++; ng[i]++; nc[i]++};
                                                                if($i == "W"){nw[i]++; na[i]++; nt[i]++};
                                                                if($i == "K"){nk[i]++; ng[i]++; nt[i]++};
                                                                if($i == "M"){nm[i]++; na[i]++; nc[i]++};
                                                                if($i == "B"){nb[i]++; nc[i]++; ng[i]++; nt[i]++};
                                                                if($i == "D"){nd[i]++; na[i]++; ng[i]++; nt[i]++};
                                                                if($i == "H"){nh[i]++; na[i]++; nc[i]++; nt[i]++};
                                                                if($i == "V"){nv[i]++; na[i]++; nc[i]++; ng[i]++};
                                                                if($i == "N"){nn[i]++; na[i]++; nt[i]++; nc[i]++; ng[i]++};
                                                                if(NR == rc){print i,nx[i],na[i],nt[i],nc[i],ng[i],nr[i],ny[i],ns[i],nw[i],nk[i],nm[i],nb[i],nd[i],nh[i],nv[i],nn[i]}}}' \
    | awk -v FS=',' -v OFS=',' '{max=$2;
                                 if(max < $3){max=$3}
                                 if(max < $4){max=$4} 
                                 if(max < $5){max=$5} 
                                 if(max < $6){max=$6}; 
                                 print $0,max }' \
    | awk -v FS=',' -v OFS=',' 'BEGIN{srand()}
                                {delete call; 
                                if ($2 == $18){call[length(call)] = "-"}; 
                                if ($3 == $18){call[length(call)] = "A"}; 
                                if ($4 == $18){call[length(call)] = "T"}; 
                                if ($5 == $18){call[length(call)] = "C"}; 
                                if ($6 == $18){call[length(call)] = "G"}; 
                                r=int(rand()*length(call));
                                print $0,call[r]}' \
    >> site-list.txt

# save files
echo ">${NAME}" > ${NAME}.fa
cat site-list.txt | tail -n +2 | cut -f 19 -d ',' | tr -d '\n -' | sed 's/$/\n/g' >> ${NAME}.fa