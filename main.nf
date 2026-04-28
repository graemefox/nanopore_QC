#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// set up and create an output directory
outdir = file(params.outdir)
outdir.mkdir()

nextflow_version="v.0.1"

params.bam   = false
params.fastq = false

if (!params.bam && !params.fastq) {
    error "ERROR: Please provide either --bam or --fastq"
}
if (params.bam && params.fastq) {
    error "ERROR: Please provide either --bam or --fastq, not both"
}

log.info """\

        ================================================================
        Oxford Nanopore Data QC - Nextflow P I P E L I N E     ${nextflow_version}
        ================================================================

        INPUTS
        ================================================================
        bam                      : ${params.bam}
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

process BAM_TO_FASTQ {
    input:
        path(bams)

    output:
        path "reads.fastq.gz", emit: reads_fastq

    script:
        """
        samtools cat ${bams} | samtools fastq -0 - | gzip -c > reads.fastq.gz
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
        BAM_TO_FASTQ_CH = BAM_TO_FASTQ(bams)
        NANOPLOT_FASTQ(BAM_TO_FASTQ_CH.reads_fastq, threads)
        WF_ALIGNMENT_FASTQ(BAM_TO_FASTQ_CH.reads_fastq, ref_dir, threads)

    } else {

        Channel.fromPath(params.fastq, checkIfExists: true).collect().set { fastq_files }
        NANOPLOT_FASTQ(fastq_files, threads)
        WF_ALIGNMENT_FASTQ(fastq_files, ref_dir, threads)

    }
}
