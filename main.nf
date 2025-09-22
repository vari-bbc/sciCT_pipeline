#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.primer_annot   = "index_files/Primer_Annotation.csv"
params.tn5_annot  = "index_files/Tn5_Barcode_Annotation.xlsx"
params.fastq1     = "raw_data/sciCHIP_L000_R1_001.fastq.gz"
params.fastq2     = "raw_data/sciCHIP_L000_R2_001.fastq.gz"
params.umi1     = "raw_data/sciCHIP_UMI_S1_L001_I1_001.fastq.gz"
params.umi2     = "raw_data/sciCHIP_UMI_S1_L001_I2_001.fastq.gz"
params.min_reads = params.min_reads ?: 10000
params.bt2_index = "/varidata/research/projects/bbc/versioned_references/latest/data/hg38_gencode/indexes/bowtie2/hg38_gencode"

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

process EXPECTED_SAMPLES {
    tag 'expected-samples'
    input:
        path primer_annot2
        path tn5_annot2
    output:
        path 'expected_samples.txt'
    """
        Expected_Samples.py -t "$tn5_annot2" -p "$primer_annot2" -o expected_samples.txt
    """
}

process PROCESS_FASTQ {
    tag "$fastq.baseName"
    input:
    tuple val(idx), path(fastq)

    output:
    tuple val(idx), path("${fastq.baseName}_mod.fastq.gz")

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

process COUNT_READS {
    tag "$sid"
    input:
        tuple val(sid), path(r1), path(r2)
    output:
        tuple val(sid), path(r1), path(r2), val(total_reads)
    """
    n1=\$(zcat "$r1" | awk 'END{print NR/4}')
    n2=\$(zcat "$r2" | awk 'END{print NR/4}')
    echo \$(( ${n1:-0} + ${n2:-0} )) > reads.txt
    """

}

process WRITE_QC {
    tag 'qc-report'
    input:
        val missing_sids
        val lows
    output:
        path 'demux_qc.tsv'
    """
    {
        echo -e "status\tsample_id\treads"
        for s in ${missing_sids[@]:-}; do echo -e "MISSING\t\$s\t0"; done
    } > demux_qc.tsv

    # append low-read lines
    while read -r line; do
        sid=\$(echo "\$line" | cut -f1)
        cnt=\$(echo "\$line" | cut -f2)
        echo -e "LOW_READS\t\$sid\t\$cnt"
    done < <(printf "%s\n" ${lows[@]})
    """
}

process CUTADAPT {
    tag "$sid"
    input:
        tuple val(sid), path(r1), path(r2)
    output:
        tuple val(sid),
            path("${sid}.trim.R1.fq.gz"),
            path("${sid}.trim.R2.fq.gz")
    """
    cutadapt -j $task.cpus -m 20 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -Z \
        -o ${sid}.trim.R1.fq.gz \
        -p ${sid}.trim.R2.fq.gz \
        "$r1" "$r2"
    """
}

process BOWTIE2 {
    tag "$sid"
    input:
        tuple val(sid), path(r1), path(r2)
    output:
        tuple val(sid), path("${sid}.sorted.bam")
    """
    bowtie2 -x "$params.bt2_index" --very-sensitive-local --soft-clipped-unmapped-tlen --no-mixed --no-discordant --dovetail --phred33 -I 10 -X 1000 \
        -1 "$r1" -2 "$r2" \
        | samtools sort -o ${sid}.sorted.bam
    """
}

workflow {
    primer_annot_ch  = Channel.of(file(params.primer_annot))
    tn5_annot_ch = Channel.of(file(params.tn5_annot))

    fastqs = Channel.of(
        tuple( 0, file(params.fastq1)),
        tuple( 1, file(params.fastq2)),
        tuple( 2, file(params.umi1)),
        tuple( 3, file(params.umi2))
    )

    primer_annot_out  = PROCESS_PRIMER_ANNOT(primer_annot_ch)
    tn5_annot_out = PROCESS_TN5_ANNOT(tn5_annot_ch)
    expected_samples_ch = EXPECTED_SAMPLES(primer_annot_out, tn5_annot_out)
        .flatMap { file('expected_samples.txt').splitText() }
        .map { it.trim() }
        .filter { it }
        .map { sid -> tuple(sid) }     

    ordered_fastqs = PROCESS_FASTQ(fastqs)
        .toList()
        .map { items -> items.sort{ a, b -> a[0] <=>b [0] } }
        .map { items -> items.collect { it[ 1 ] } }

    RUN_DEMUX(primer_annot_out, tn5_annot_out, ordered_fastqs)

    dmx_gen_out = demux_ch
        .filter { f -> !(f.name =~ /Undetermined/i) }
        .map { f ->
            def m = (f.baseName =~ /^(.+)_R([12])$/)
            assert m, "Unexpected demux filename: ${f.name}"
            def sid  = m[0][1]
            def mate = m[0][2] as int
            tuple(sid, mate, f)
    }
    paired_gen_out = dmx_gen_out
        .groupTuple()
        .map { sid, items ->
            def r1 = items.find { it[1] == 1 }?.getAt(2)
            def r2 = items.find { it[1] == 2 }?.getAt(2)
            assert r1 && r2 : "Missing mate for ${sid}"
            tuple(sid, r1, r2)
    }

    expected_list = expected_samples_ch.map{ it[0] }.collect()
    found_list    = paired_gen_out.map{ sid, r1, r2 -> sid }.collect()

    missing_ids = expected_list.combine(found_list)
        .map { exp, got -> exp - got }
    
    read_counts = COUNT_READS(paired_gen_out)

    valid_pairs = read_counts.filter { sid, r1, r2, n -> n >= params.min_reads }
                            .map    { sid, r1, r2, n -> tuple(sid, r1, r2) }

    low_pairs   = read_counts.filter { sid, r1, r2, n -> n <  params.min_reads }
                            .map    { sid, r1, r2, n -> tuple(sid, n) }

    missing_vec = missing_ids.map { arr -> arr as String[] } 
    low_vec     = low_pairs.collect().map { lst -> lst.collect{ it[0] + "\t" + it[1] } as String[] }

    qc_report = WRITE_QC(missing_vec, low_vec)

    trimmed_pairs = CUTADAPT(valid_pairs)
    aligned = BOWTIE2(trimmed_pairs)




}
