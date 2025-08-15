#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.primer_annot   = "index_files/Primer_Annotation.csv"
params.meta_xlsx  = "metadata.xlsx"
params.fastq1     = "reads_1.fastq.gz"
params.fastq2     = "reads_2.fastq.gz"
params.fastq3     = "reads_3.fastq.gz"
params.fastq4     = "reads_4.fastq.gz"

params.output_dir = "results"

process PROCESS_PRIMER_ANNOT {
    tag "$primer_annot"
    input:
    path primer_annot from file(params.primer_annot)

    output:
    path "${primer_annot.baseName}_2.csv"

    script:
    """
    conda activate sciCT_env
    python AddPrimerAnnotationCols.py -i $primer_annot -o "${primer_annot.baseName}_2.csv"
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
    path primer_annot2   from PROCESS_PRIMER_ANNOT.out
    path meta_xlsx2  from PROCESS_METADATA_XLSX.out
    path processed_fastqs from PROCESS_FASTQ.out.collect()

    script:
    """
    conda activate my_extract_env
    python extract.py $primer_annot2 $meta_xlsx2 ${processed_fastqs.join(" ")}
    """
}

workflow {
    primer_annot_ch  = Channel.of(file(params.primer_annot))
    meta_xlsx_ch = Channel.of(file(params.meta_xlsx))

    fastqs = Channel.of(
        file(params.fastq1),
        file(params.fastq2),
        file(params.fastq3),
        file(params.fastq4)
    )

    primer_annot_out  = PROCESS_PRIMER_ANNOT(primer_annot_ch)
    meta_xlsx_out = PROCESS_METADATA_XLSX(meta_xlsx_ch)
    fastq_out     = PROCESS_FASTQ(fastqs)

    RUN_EXTRACTION(primer_annot_out, meta_xlsx_out, fastq_out)
}
