#!/bin/bash

cat $1 | sed 's/>.*$/@&@/g' | tr -d '\n' | tr '@' '\n' | grep -v '>' | tail -n +2 | awk '{print toupper($0)}' | sort | uniq | grep -vE 'R|Y|M|K|S|W|H|B|V|D|N' | awk '{print ">"NR"\n"$0}' > no-dups.fa