nextflow.enable.dsl=2

log.info """\
    O N C O T R A N S C R I P T - T N B C   
    =====================================================
    Samplesheet  : ${params.samplesheet}
    Genome FASTA : ${params.genome_fasta}
    GTF GFF File : ${params.gtf}
    STAR Index   : ${params.star_index}
    Outputs Dir  : ${params.outdir}
    =====================================================
"""
.stripIndent()

process FASTQC {
    tag "QC on ${sample_id}"
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}", emit: qc_files

    script:
    """
    fastqc --threads ${task.cpus} ${reads[0]} ${reads[1]}
    """
}

process TRIMGALORE {
    tag "Trimming ${sample_id}"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*.fq.gz"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.fq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt", emit: reports

    script:
    """
    trim_galore \\
        --paired \\
        --cores ${task.cpus} \\
        --fastqc \\
        --gzip \\
        ${reads[0]} \\
        ${reads[1]}
    """
}

process STAR_ALIGN {
    tag "Aligning ${sample_id}"
    publishDir "${params.outdir}/aligned", mode: 'copy', pattern: "*.sortedByCoord.out.bam"

    input:
    tuple val(sample_id), path(reads)
    path star_index

    output:
    tuple val(sample_id), path("*.sortedByCoord.out.bam"), emit: bam
    path "*.Log.final.out", emit: log_final

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --quantMode TranscriptomeSAM GeneCounts
    """
}

process FEATURE_COUNTS {
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path bams
    path gtf

    output:
    path "merged_gene_counts.txt", emit: matrix
    path "merged_gene_counts.txt.summary", emit: summary

    script:
    """
    featureCounts \\
        -p \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -t exon \\
        -g gene_id \\
        -o merged_gene_counts.txt \\
        ${bams}
    """
}

process MULTIQC {
    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    path reports

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc . --name "multiqc_report.html"
    """
}

workflow {
    ch_samplesheet = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true, strip:true)
        .map { row -> tuple(row.sample_id, [file(row.fastq_1), file(row.fastq_2)]) }

    ch_star_index = file(params.star_index)
    ch_gtf        = file(params.gtf)

    FASTQC(ch_samplesheet)
    TRIMGALORE(ch_samplesheet)

    STAR_ALIGN(TRIMGALORE.out.trimmed_reads, ch_star_index)

    ch_all_bams = STAR_ALIGN.out.bam.map { it[1] }.collect()
    FEATURE_COUNTS(ch_all_bams, ch_gtf)

    ch_qc_reports = FASTQC.out.qc_files.collect()
        .mix(TRIMGALORE.out.reports.collect())
        .mix(STAR_ALIGN.out.log_final.collect())
        .mix(FEATURE_COUNTS.out.summary)
        .collect()

    MULTIQC(ch_qc_reports)
}
