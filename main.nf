#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.meta_csv   = "metadata.csv"
params.meta_xlsx  = "metadata.xlsx"
params.fastq1     = "reads_1.fastq.gz"
params.fastq2     = "reads_2.fastq.gz"
params.fastq3     = "reads_3.fastq.gz"
params.fastq4     = "reads_4.fastq.gz"

params.output_dir = "results"

process PROCESS_METADATA_CSV {
    tag "$meta_csv"
    input:
    path meta_csv from file(params.meta_csv)

    output:
    path "metadata1_2.csv"

    script:
    """
    python script1.py $meta_csv metadata1_2.csv
    """
}

process PROCESS_METADATA_XLSX {
    tag "$meta_xlsx"
    input:
    path meta_xlsx from file(params.meta_xlsx)

    output:
    path "metadata2_2.csv"

    script:
    """
    python script2.py $meta_xlsx metadata2_2.csv
    """
}

process PROCESS_FASTQ {
    tag "$fastq"
    input:
    path fastq from fastqs

    output:
    path "processed_${fastq.getName()}"

    script:
    """
    bash modify_fastq.sh $fastq processed_${fastq.getName()}
    """
}

process RUN_EXTRACTION {
    publishDir params.output_dir, mode: 'copy'

    input:
    path meta_csv2   from PROCESS_METADATA_CSV.out
    path meta_xlsx2  from PROCESS_METADATA_XLSX.out
    path processed_fastqs from PROCESS_FASTQ.out.collect()

    script:
    """
    conda activate my_extract_env
    python extract.py $meta_csv2 $meta_xlsx2 ${processed_fastqs.join(" ")}
    """
}

workflow {
    meta_csv_ch  = Channel.of(file(params.meta_csv))
    meta_xlsx_ch = Channel.of(file(params.meta_xlsx))

    fastqs = Channel.of(
        file(params.fastq1),
        file(params.fastq2),
        file(params.fastq3),
        file(params.fastq4)
    )

    meta_csv_out  = PROCESS_METADATA_CSV(meta_csv_ch)
    meta_xlsx_out = PROCESS_METADATA_XLSX(meta_xlsx_ch)
    fastq_out     = PROCESS_FASTQ(fastqs)

    RUN_EXTRACTION(meta_csv_out, meta_xlsx_out, fastq_out)
}
