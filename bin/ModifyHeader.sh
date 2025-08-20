#!/bin/bash

module load bbc2/pigz/pigz-2.7

FileName="$1"
OutName="$2"
read_pair=${BaseName#*_[R,I]}   # Remove everything up to and including _R or _I
read_pair=${read_pair%%_*}  # Remove everything from next _ onwards

pigz -dc ${FileName}  | \
awk -v readn="$read_pair" 'BEGIN{FS=OFS=":"} /^@AV234602/ {$8=""; sub(/::/, " " readn ":", $0) } { print }' | \
pigz > ${OutName}
