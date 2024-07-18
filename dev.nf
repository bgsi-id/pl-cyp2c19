#!/usr/bin/env nextflow

// Define parameters
params.input_dir = params.input_dir ?: ''
params.ref = 'pl-cyp2c19/static/GRCh38.cyp2c19.fa'
params.threads = "12"
params.model = "r1041_e82_400bps_hac_v430"
params.bed = "pl-cyp2c19/static/region.bed"
params.output = "pl-cyp2c19_output"


include { SAMTOOLS_FASTQ }                       from '../modules/nf-core/samtools/fastq/main.nf'
include { SAMTOOLS_SORT }                        from '../modules/nf-core/samtools/sort/main.nf'
include { MINIMAP2_ALIGN }                       from '../modules/nf-core/minimap2/align/main.nf'

// Workflow definition
// Workflow definition
workflow ALIGN_WORKFLOW {
    def inputFiles = Channel.fromPath("${params.input_dir}/*")

    // Process files based on their extensions
    inputFiles.view { file ->
        def fileName = file.toString()
        if (fileName.endsWith('.fastq.gz') || fileName.endsWith('.fastq')) {
            def fastqFilesChannel = Channel.fromPath(file)
            MINIMAP2_ALIGN(input: fastqFilesChannel, index: params.ref) | SAMTOOLS_SORT
        } else if (fileName.endsWith('.bam')) {
            def bamFilesChannel = Channel.fromPath(file)
            bamFilesChannel | SAMTOOLS_FASTQ | MINIMAP2_ALIGN(index: params.ref) | SAMTOOLS_SORT
        }
    }

    // Check if no files were found and exit if true
    inputFiles.subscribe { file ->
        if (!file) {
            println "No fastq, fastq.gz or bam files found in ${params.input_dir}"
            System.exit(1)
        }
    }
}