#!/bin/bash

module load bbc2/pigz/pigz-2.7

DirName='raw_data'

FileName="$1"
read_pair=${FileName#*_R}   # Remove everything up to and including _R
read_pair=${read_pair%%_*}  # Remove everything from next _ onwards

pigz -dc ${DirName}/${FileName}  | \
awk -v readn="$read_pair" 'BEGIN{FS=OFS=":"} /^@AV234602/ {$8=""; sub(/::/, " " readn ":", $0) } { print }' | \
pigz > ${DirName}_mod/${FileName}
