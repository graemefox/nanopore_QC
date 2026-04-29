#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// set up and create an output directory
outdir = file(params.outdir)
outdir.mkdir()

nextflow_version="v.0.1"

params.bam                  = false
params.fastq                = false
params.downsample           = false   // fraction to keep e.g. 0.1 for 10%. false = no downsampling
params.EPI2ME_profile = 'standard'

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
        downsample               : ${params.downsample}
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

// Strip alignment info per BAM in parallel — only if the BAM has been aligned.
// @SQ lines in the header indicate alignment to a reference; absent = already unmapped.
process RESET_BAM {
    input:
        path(bam)

    output:
        path("${bam.baseName}.reset.bam"), emit: reset_bam

    script:
        """
        if samtools view -H ${bam} | grep -q "^@SQ"; then
            echo "Aligned BAM detected — resetting ${bam}"
            samtools reset -x OA,XA --no-PG ${bam} -o ${bam.baseName}.reset.bam
        else
            echo "Unaligned BAM — skipping reset for ${bam}"
            ln -s \$(realpath ${bam}) ${bam.baseName}.reset.bam
        fi
        """
}

// Downsample per reset BAM in parallel using a read-name hash.
// -s seed.fraction: safe on coordinate-sorted BAMs (hashes read name, not position).
process DOWNSAMPLE_BAM {
    input:
        path(bam)

    output:
        path("${bam.baseName}.ds.bam"), emit: downsampled_bam

    script:
        def fraction = (params.downsample * 100).toInteger()
        """
        samtools view -s 42.${fraction} -b ${bam} -o ${bam.baseName}.ds.bam
        """
}

// Downsample per FASTQ in parallel using seqtk sample (seed 42, deterministic).
process DOWNSAMPLE_FASTQ {
    input:
        path(fastq)

    output:
        path("${fastq.simpleName}.ds.fq.gz"), emit: downsampled_fastq

    script:
        """
        seqtk sample -s 42 ${fastq} ${params.downsample} | gzip > ${fastq.simpleName}.ds.fq.gz
        """
}

// NanoPlot from BAM input (full data — no downsampling)
process NANOPLOT_BAM {
    input:
        path(bams)
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
        --ubam ${bams} \
        --raw
        """
}

// NanoPlot from FASTQ input
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

// wf-alignment from BAM input 
process WF_ALIGNMENT_BAM {
    input:
        path(bams)
        path(ref_dir)
        val(threads)

    publishDir("${params.outdir}/wf-alignment",      mode: 'copy', pattern: "wf-alignment-report.html")
    publishDir("${params.outdir}/wf-alignment/data", mode: 'copy', pattern: "**.{hist,tsv,json}")

    output:
        path("wf-alignment-report.html")
        path("**.{hist,tsv,json}")

    script:
        """

        nextflow run epi2me-labs/wf-alignment -r master \
        -profile ${params.EPI2ME_profile} \
        -ansi-log false \
        -work-dir ./nf-work \
        --references ${ref_dir} \
        --bam ${bams} \
        --threads ${threads} \
        --depth_coverage false \
        --out_dir .
        """
}

// wf-alignment from FASTQ input
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
        nextflow run epi2me-labs/wf-alignment -r master \
        -profile ${params.EPI2ME_profile} \
        -ansi-log false \
        -work-dir ./nf-work \
        --references ${ref_dir} \
        --fastq ${fastq_files} \
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

    GET_VERSIONS_CH = GET_VERSIONS()

    if (params.bam) {
// Working with BAMs:
        Channel.fromPath(params.bam, checkIfExists: true).set { bams_ch }

        // Reset all BAMs in parallel (strips alignment info from pre-aligned BAMs)
        RESET_BAM_CH = RESET_BAM(bams_ch)

        // NanoPlot always gets full (non-downsampled) reset BAMs for accurate stats
        NANOPLOT_BAM_CH = NANOPLOT_BAM(RESET_BAM_CH.reset_bam, threads)

        // wf-alignment gets downsampled BAMs if --downsample is set, otherwise full
        if (params.downsample) {
        // Yes, we're downsampling
            DOWNSAMPLE_BAM_CH = DOWNSAMPLE_BAM(RESET_BAM_CH.reset_bam)

            WF_ALIGNMENT_BAM_CH = WF_ALIGNMENT_BAM(DOWNSAMPLE_BAM_CH.downsampled_bam, ref_dir, threads)

        } else {
        // No downsampling
              WF_ALIGNMENT_BAM_CH = WF_ALIGNMENT_BAM(RESET_BAM_CH.reset_bam, ref_dir, threads)
        }

    } else {
//        Working with FastQs....
        Channel.fromPath(params.fastq, checkIfExists: true).set { fastq_ch }

        // NanoPlot always gets the full (non-downsampled) data for accurate stats
        NANOPLOT_FASTQ_CH = NANOPLOT_FASTQ(fastq_ch, threads)

        // wf-alignment gets downsampled FASTQs if --downsample is set, otherwise full
        if (params.downsample) {
        // Yes, we're downsampling
            DOWNSAMPLE_FASTQ_CH = DOWNSAMPLE_FASTQ(fastq_ch)
            WF_ALIGNMENT_FASTQ_CH = WF_ALIGNMENT_FASTQ(DOWNSAMPLE_FASTQ_CH.downsampled_fastq, ref_dir, threads)

        } else {
        // No downsampling
            WF_ALIGNMENT_FASTQ_CH = WF_ALIGNMENT_FASTQ(fastq_files, ref_dir, threads)
        }
    }
}
