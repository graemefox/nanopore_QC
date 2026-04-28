#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// set up and create an output directory
outdir = file(params.outdir)
outdir.mkdir()

nextflow_version="v.0.1"

params.bam    = false
params.fastq   = false
params.my_param = false

if (!params.bam && !params.fastq) {
    error "ERROR: Please provide either --bams or --fastq"
}
if (params.bam && params.fastq) {
    error "ERROR: Please provide either --bams or --fastq, not both"
}

log.info """\

        ================================================================
        Oxford Nanopore Data QC - Nextflow P I P E L I N E     ${nextflow_version}
        ================================================================

        INPUTS
        ================================================================
        bams                     : ${params.bam}
        fastq                    : ${params.fastq}
        ref_dir                  : ${params.ref_dir}
        threads                  : ${params.threads}
        outdir                   : ${params.outdir}

        """
        .stripIndent()

process GET_VERSIONS {
    input:

    publishDir("${params.outdir}", mode: 'copy')

    output:
        path "versions.txt", emit: versions

    script:
        """
        samtools --version | head -n 1 >> versions.txt
        NanoPlot -v >> versions.txt
        nextflow -v >> versions.txt
        nextflow run epi2me-labs/wf-alignment | grep '^Launching' | awk '{for (i=1; i<=NF; i++) if (\$i ~ /^revision:/) print "epi2me-labs/wf-alignment " \$i " " \$(i+1)}' >> versions.txt
        """
}

process CAT_BAMS {
    input:
        path(bams)
        val(threads)

    output:
       path "merged.sorted.bam", emit: merged_bam

    script:
        """
        samtools cat -o merged.bam ${bams}
        samtools sort -@${threads} -o merged.sorted.bam merged.bam
        """
}

process INDEX_MERGED_BAM {
    input:
        path(input_bam)
        val(threads)

    output:
        path "*.bai", emit: indexed_bam
        tuple path(input_bam), path("*.bai"), emit: indexed_bam_tuple

    script:
        """
        samtools \
        index \
        -@${threads} \
        ${input_bam}
        """
}

process NANOPLOT_BAM {
    input:
        path(input_bam)
        val(threads)
        path(indexed_bam)

    publishDir("${params.outdir}/NanoPlot",      mode: 'copy', pattern: "NanoPlot-report.html")
    publishDir("${params.outdir}/NanoPlot",      mode: 'copy', pattern: "NanoStats.txt")
    publishDir("${params.outdir}/NanoPlot",      mode: 'copy', pattern: "NanoPlot-data.tsv.gz")
    publishDir("${params.outdir}/NanoPlot/pngs", mode: 'copy', pattern: "*.png")

    output:
        path "NanoPlot-report.html", emit: nanoplot_report
        path "NanoStats.txt",        emit: nanoplot_stats
        path "NanoPlot-data.tsv.gz", emit: nanoplot_data
        path "*.png"

    script:
        """
        NanoPlot \
        -t ${threads} \
        --bam ${input_bam} \
        --raw
        """
}

process NANOPLOT_FASTQ {
    input:
        path(fastq_files)
        val(threads)

    publishDir("${params.outdir}/NanoPlot",      mode: 'copy', pattern: "NanoPlot-report.html")
    publishDir("${params.outdir}/NanoPlot",      mode: 'copy', pattern: "NanoStats.txt")
    publishDir("${params.outdir}/NanoPlot",      mode: 'copy', pattern: "NanoPlot-data.tsv.gz")
    publishDir("${params.outdir}/NanoPlot/pngs", mode: 'copy', pattern: "*.png")

    output:
        path "NanoPlot-report.html", emit: nanoplot_report
        path "NanoStats.txt",        emit: nanoplot_stats
        path "NanoPlot-data.tsv.gz", emit: nanoplot_data
        path "*.png"

    script:
        """
        NanoPlot \
        -t ${threads} \
        --fastq ${fastq_files} \
        --raw
        """
}

process WF_ALIGNMENT_BAM {
    input:
        path(merged_bam)
        path(ref_dir)
        val(threads)

    publishDir("${params.outdir}/wf-alignment",      mode: 'copy', pattern: "wf-alignment-report.html")
    publishDir("${params.outdir}/wf-alignment/data", mode: 'copy', pattern: "**.{hist,tsv,json}")

    output:
        path("wf-alignment-report.html")
        path("**.{hist,tsv,json}")

    script:
        """
        samtools fastq -@ ${threads} -0 reads.fastq.gz -c 6 ${merged_bam}

        nextflow run epi2me-labs/wf-alignment -r master \
        -profile standard \
        -ansi-log false \
        -work-dir ./nf-work \
        --references ${ref_dir} \
        --fastq reads.fastq.gz \
        --threads ${threads} \
        --depth_coverage false \
        --out_dir .
        """
}

process WF_ALIGNMENT_FASTQ {
    input:
        path(fastq_files)
        path(ref_dir)
        val(threads)

    publishDir("${params.outdir}/wf-alignment",      mode: 'copy', pattern: "wf-alignment-report.html")
    publishDir("${params.outdir}/wf-alignment/data", mode: 'copy', pattern: "**.{hist,tsv,json}")

    output:
        path("wf-alignment-report.html")
        path("**.{hist,tsv,json}")

    script:
        """
        mkdir -p fastq_input
        for f in ${fastq_files}; do ln -s \$(realpath \$f) fastq_input/; done

        nextflow run epi2me-labs/wf-alignment -r master \
        -profile standard \
        -ansi-log false \
        -work-dir ./nf-work \
        --references ${ref_dir} \
        --fastq fastq_input \
        --threads ${threads} \
        --depth_coverage false \
        --out_dir .
        """
}

///////////////////////////
// MAIN WORKFLOW SECTION //
///////////////////////////

workflow {

    Channel.fromPath(params.ref_dir, checkIfExists: true).set { ref_dir }
    Channel.from(params.threads).set { threads }

    GET_VERSIONS()

    if (params.bam) {

        Channel.fromPath(params.bam, checkIfExists: true).collect().set { bams }

        CAT_BAMS_CH        = CAT_BAMS(bams, threads)
        INDEX_MERGED_BAM_CH = INDEX_MERGED_BAM(CAT_BAMS_CH.merged_bam, threads)

        NANOPLOT_BAM(CAT_BAMS_CH.merged_bam, threads, INDEX_MERGED_BAM_CH.indexed_bam)
        WF_ALIGNMENT_BAM(CAT_BAMS_CH.merged_bam, ref_dir, threads)

    } else {

        Channel.fromPath(params.fastq, checkIfExists: true).collect().set { fastq_files }

        NANOPLOT_FASTQ(fastq_files, threads)
        WF_ALIGNMENT_FASTQ(fastq_files, ref_dir, threads)

    }
}
