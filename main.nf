#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.primer_annot   = "index_files/Primer_Annotation.csv"
params.tn5_annot  = "index_files/Tn5_Barcode_Annotation.xlsx"
params.fastq1     = "raw_data/sciCHIP_L000_R1_001.fastq.gz"
params.fastq2     = "raw_data/sciCHIP_L000_R2_001.fastq.gz"
params.umi1     = "raw_data/sciCHIP_UMI_S1_L001_I1_001.fastq.gz"
params.umi2     = "raw_data/sciCHIP_UMI_S1_L001_I2_001.fastq.gz"

params.output_dir = "results"

process PROCESS_PRIMER_ANNOT {
    tag "$primer_annot"
    input:
    path primer_annot

    output:
    path "${primer_annot.baseName}_2.csv"

    script:
    """
    AddPrimerAnnotationCols.py -i $primer_annot -o "${primer_annot.baseName}_2.csv"
    dos2unix "${primer_annot.baseName}_2.csv"
    """
}

process PROCESS_TN5_ANNOT {
    tag "$tn5_annot"
    input:
    path tn5_annot

    output:
    path "${tn5_annot.baseName}_2.csv"

    script:
    """
    Check_Tn5_file.py -i $tn5_annot -o "${tn5_annot.baseName}_2.csv"
    dos2unix "${tn5_annot.baseName}_2.csv"
    """
}

process PROCESS_FASTQ {
    tag "$fastq.baseName"
    input:
    path fastq

    output:
    path "${fastq.baseName}_mod.fastq.gz"

    script:
    """
    bash ModifyHeader.sh $fastq "${fastq.baseName}_mod.fastq.gz"
    """
}

process RUN_DEMUX {

    input:
    path primer_annot2   
    path tn5_annot2  
    val processed_fastqs 

    output: 
    path "demux_out/*.fq.gz", emit: demux_ch

    script:
    """
    mkdir -p demux_out
    sciCTextract --forward-mode \
      --outdir demux_out \
      --Primer_Barcode $primer_annot2 \
      --Tn5_Barcode $tn5_annot2 \
      ${processed_fastqs.join(" ")}
    """
}

workflow {
    primer_annot_ch  = Channel.of(file(params.primer_annot))
    tn5_annot_ch = Channel.of(file(params.tn5_annot))

    fastqs = Channel.of(
        file(params.fastq1),
        file(params.fastq2),
        file(params.umi1),
        file(params.umi2)
    )

    primer_annot_out  = PROCESS_PRIMER_ANNOT(primer_annot_ch)
    tn5_annot_out = PROCESS_TN5_ANNOT(tn5_annot_ch)
    fastq_out     = PROCESS_FASTQ(fastqs)

    RUN_DEMUX(primer_annot_out, tn5_annot_out, fastq_out.collect())
}
