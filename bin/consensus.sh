#!/bin/bash
set -euo pipefail

version="2.1"

# consensus.sh
# Author: Jared Johnson, jared.johnson@doh.wa.gov

# Show version
if [[ "${1:-}" == "version" ]]; then
    echo "${version}"
    exit 0
fi

# Input args
NAME=$1
ALN=$2

# Detect if input is gzipped and use appropriate command
if [[ "${ALN}" == *.gz ]]; then
    READ_CMD="gunzip -c"
else
    READ_CMD="cat"
fi

# Reformat alignment: single-line sequences, uppercase, no headers
${READ_CMD} "${ALN}" \
    | sed 's/>.*$/@&@/g' \
    | tr -d '\n' \
    | tr '@' '\n' \
    | grep -v '>' \
    | tail -n +2 \
    | awk '{print toupper($0)}' > seq.txt

# Check for illegal characters
bad_chars=$(tr -d '\-ATCGRYSWKMBDHVN\n' < seq.txt)
if [[ -n "${bad_chars}" ]]; then
    echo "Error: the alignment contains illegal characters: ${bad_chars}"
    exit 1
fi

# Alignment stats
seq_len=$(sed -n 1p seq.txt | tr -d '\t\r\n ' | wc -c)
seq_count=$(wc -l < seq.txt)
echo "Alignment Length: ${seq_len}"
echo "Number of Sequences: ${seq_count}"

# Generate site-wise consensus calls
echo "site,-,A,T,C,G,R,Y,S,W,K,M,B,D,H,V,N,max,call" > site-list.txt

gawk -v OFS=',' -v rc=${seq_count} '
{
    for (i = 1; i <= length($0); i++) {
        base = substr($0, i, 1)
        a[i]["-"] += (base == "-")
        a[i]["A"] += (base == "A")
        a[i]["T"] += (base == "T")
        a[i]["C"] += (base == "C")
        a[i]["G"] += (base == "G")
        a[i]["R"] += (base == "R")
        a[i]["Y"] += (base == "Y")
        a[i]["S"] += (base == "S")
        a[i]["W"] += (base == "W")
        a[i]["K"] += (base == "K")
        a[i]["M"] += (base == "M")
        a[i]["B"] += (base == "B")
        a[i]["D"] += (base == "D")
        a[i]["H"] += (base == "H")
        a[i]["V"] += (base == "V")
        a[i]["N"] += (base == "N")
    }
}
END {
    for (i = 1; i <= length(a); i++) {
        dash = a[i]["-"] + 0
        A = a[i]["A"] + a[i]["R"] + a[i]["W"] + a[i]["M"] + a[i]["D"] + a[i]["H"] + a[i]["V"] + a[i]["N"]
        T = a[i]["T"] + a[i]["Y"] + a[i]["W"] + a[i]["K"] + a[i]["B"] + a[i]["D"] + a[i]["H"] + a[i]["N"]
        C = a[i]["C"] + a[i]["Y"] + a[i]["S"] + a[i]["M"] + a[i]["B"] + a[i]["H"] + a[i]["V"] + a[i]["N"]
        G = a[i]["G"] + a[i]["R"] + a[i]["S"] + a[i]["K"] + a[i]["B"] + a[i]["D"] + a[i]["V"] + a[i]["N"]
        max = dash
        if (A > max) max = A
        if (T > max) max = T
        if (C > max) max = C
        if (G > max) max = G

        # Randomly resolve ties
        split("", calls)
        if (dash == max) calls[length(calls)+1] = "-"
        if (A == max) calls[length(calls)+1] = "A"
        if (T == max) calls[length(calls)+1] = "T"
        if (C == max) calls[length(calls)+1] = "C"
        if (G == max) calls[length(calls)+1] = "G"
        call = calls[int(rand()*length(calls)) + 1]

        print i, dash, A, T, C, G, a[i]["R"]+0, a[i]["Y"]+0, a[i]["S"]+0, a[i]["W"]+0, a[i]["K"]+0, a[i]["M"]+0, a[i]["B"]+0, a[i]["D"]+0, a[i]["H"]+0, a[i]["V"]+0, a[i]["N"]+0, max, call
    }
}' seq.txt >> site-list.txt

# Output consensus
echo ">${NAME}" > "${NAME}.fa"
cut -d',' -f19 site-list.txt | tail -n +2 | tr -d '\n -' | sed 's/$/\n/' >> "${NAME}.fa"