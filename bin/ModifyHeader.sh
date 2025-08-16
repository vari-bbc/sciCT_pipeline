#!/bin/bash

module load bbc2/pigz/pigz-2.7

FileName="$1"
BaseName=$(basename "${FileName}")
DirName=$(dirname "${FileName}")
read_pair=${BaseName#*_[R,I]}   # Remove everything up to and including _R or _I
read_pair=${read_pair%%_*}  # Remove everything from next _ onwards

mkdir -p ${DirName}_mod
pigz -dc ${FileName}  | \
awk -v readn="$read_pair" 'BEGIN{FS=OFS=":"} /^@AV234602/ {$8=""; sub(/::/, " " readn ":", $0) } { print }' | \
pigz > ${DirName}_mod/${BaseName}
