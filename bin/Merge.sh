#!/bin/bash

DirName1="$1"
DirName2="$2"
mkdir -p merged
for f in $( ls ${DirName1}/*.fastq.gz) ; do
    name=$(basename ${f})
    cat ${f} ${DirName2}/${name} > merged/${name}
done