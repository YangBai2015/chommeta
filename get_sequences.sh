#!/bin/bash

seq_list_file=$1
fasta_file=$2

for seq in $(cat $seq_list_file)
do

    awk "/>$seq$/ {flag=1;print;next} />/ {flag=0} flag {print}" $fasta_file

done

