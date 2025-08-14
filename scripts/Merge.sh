#!/bin/bash

# for f in $( ls PR001746_JAND_sciCT_mod/*.fastq.gz) ; do
#     echo "cat PR001746_JAND_sciCT_mod/${f} PR001746_JAND_sciCT_2_mod/${f} > PR001746_JAND_sciCT_merged_mod/${f}"
# done


cat PR001746_JAND_sciCT_mod/DJCoCnT060625_S2_L002_R1_001.fastq.gz PR001746_JAND_sciCT_2_mod/DJCoCnT060625_L000_R1_001.fastq.gz \
  > PR001746_JAND_sciCT_merged_mod/DJCoCnT060625_L002_R1_001.fastq.gz

cat PR001746_JAND_sciCT_mod/DJCoCnT060625_S2_L002_R2_001.fastq.gz PR001746_JAND_sciCT_2_mod/DJCoCnT060625_L000_R2_001.fastq.gz \
  > PR001746_JAND_sciCT_merged_mod/DJCoCnT060625_L002_R2_001.fastq.gz

cat PR001746_JAND_sciCT_mod/DJCoCnT060625_UMI_S2_L002_I1_001.fastq.gz PR001746_JAND_sciCT_2_mod/DJCoCnT060625_UMI_S1_L001_I1_001.fastq.gz \
  PR001746_JAND_sciCT_2_mod/DJCoCnT060625_UMI_S1_L002_I1_001.fastq.gz > PR001746_JAND_sciCT_merged_mod/DJCoCnT060625_UMI_L000_I1_001.fastq.gz

cat PR001746_JAND_sciCT_mod/DJCoCnT060625_UMI_S2_L002_I2_001.fastq.gz PR001746_JAND_sciCT_2_mod/DJCoCnT060625_UMI_S1_L001_I2_001.fastq.gz \
  PR001746_JAND_sciCT_2_mod/DJCoCnT060625_UMI_S1_L002_I2_001.fastq.gz > PR001746_JAND_sciCT_merged_mod/DJCoCnT060625_UMI_L000_I2_001.fastq.gz