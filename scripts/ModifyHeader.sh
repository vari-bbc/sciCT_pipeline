#!/bin/bash

module load bbc2/pigz/pigz-2.7

DirName='PR001746_JAND_sciCT'

mkdir -p ${DirName}_mod

# for f in $(ls ${DirName}/*_001.fastq.gz) ; do
#     FileName=$(basename ${f})   
#     pigz -dc $f \
#     | sed -u '1~4 s/^\(@[^:]\+\(:[^:]\+\)\{6\}\):[^ ]\+\( .*\)/\1\3/' \
#     | pigz > ${DirName}_mod/${FileName}   
# done




for n in 1 2 ; do 
    for f in $(ls ${DirName}/*${n}_001.fastq.gz) ; do
        FileName=$(basename ${f})
        pigz -dc ${f}  | \
        awk -v readn="$n" 'BEGIN{FS=OFS=":"} /^@AV234602/ {$8=""; sub(/::/, " " readn ":", $0) } { print }' | \
        pigz > ${DirName}_mod/${FileName}
    done
done

# pigz -dc ${DirName}/${DirName}_JAND_sciCandT/

# pigz -dc PR001626_JAND_sciCandT_2/sciCUTTAG_S1_L001_R1_001.fastq.gz | \
#   awk 'BEGIN{FS=OFS=":"} /^@AV234602/ {$8=""; sub(/::/, " 1:", $0) } { print }' | \
#   pigz > PR001626_JAND_sciCandT_2_mod/sciCUTTAG_S1_L001_R1_001.fastq.gz