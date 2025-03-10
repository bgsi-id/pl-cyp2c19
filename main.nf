#!/usr/bin/env nextflow

// Define parameters
params.input_dir = params.input_dir ?: ''
params.ref = 'pl-cyp2c19/static/GRCh38.cyp2c19.fa'
params.threads = "12"
params.model = "r1041_e82_400bps_hac_v430"
params.bed = "pl-cyp2c19/static/region.bed"
params.output = "pl-cyp2c19_output"

// Extract the directory name to use as the sample name

input_ch = Channel.fromPath("${params.input_dir}/*", type:'dir')

// Processes
process alignmentProcess {
    publishDir "$params.output/$input_dir"

    input:
    path input_dir
    path ref

    output:
    tuple path("${input_dir}.bam"), path("${input_dir}.bam.bai")

    script:
    """
    if ls ${input_dir}/*.fastq.gz 1> /dev/null 2>&1; then
        zcat ${input_dir}/*.fastq.gz | minimap2 -ax map-ont ${ref} - | samtools sort -o "${input_dir}.bam" - && samtools index "${input_dir}.bam"
    elif ls ${input_dir}/*.fastq 1> /dev/null 2>&1; then
        samtools fastq ${input_dir}/*.fastq | minimap2 -ax map-ont ${ref} - | samtools sort -o "${input_dir}.bam" - && samtools index "${input_dir}.bam"
    elif ls ${input_dir}/*.bam 1> /dev/null 2>&1; then
        samtools fastq ${input_dir}/*.bam | minimap2 -ax map-ont ${ref} - | samtools sort -o "${input_dir}.bam" - && samtools index "${input_dir}.bam"
    else
        echo "No fastq, fastq.gz or bam files found in ${input_dir}"
        exit 1
    fi
    """
}

process index {
    input:
    path ref

    output:
    path "GRCh38.cyp2c19.fa.fai"

    script:
    """
    samtools faidx ${ref}
    """
}

process variantCalling {
    publishDir "${params.output}/${bam.simpleName}"
    input:
    path(bam)
    path(bam_index)
    path(ref)
    path(ref_index)
    path(model)

    output:
    path "clair3/*"

    script:
    """
    run_clair3.sh \
        --bam_fn=${bam} \
        --ref_fn=${ref} \
        --platform="ont" \
        --threads=${params.threads} \
        --model_path=${model} \
        --include_all_ctgs \
        --var_pct_full=1 \
        --var_pct_phasing=1 \
        --gvcf \
        --output=clair3
    """
}

process liftoverToHg38 {
    publishDir "${params.output}/${bam.simpleName}"
    input:
    path bam
    path vcf
    path gvcf

    output:
    tuple path("${bam.simpleName}.hg38.bam"), path("${bam.simpleName}.hg38.vcf"), path("${bam.simpleName}.hg38.gvcf"), path("${bam.simpleName}.hg38.bam.bai")

    script:
    """
    echo "Starting liftover for BAM: ${bam} VCF: ${vcf} GVCF: ${gvcf}"

    samtools view -H ${bam} | sed 's/SN:chr10:94757681-94855547/SN:chr10/' > header.txt && \
    samtools view ${bam} | sed 's/chr10:94757681-94855547/chr10/' - | \
    awk -v shift=94757680 '\$3 == "chr10" {OFS="\\t"; \$4 = \$4 + shift} {print \$0}' > body.txt && \
    cat body.txt >> header.txt && \
    samtools view -b header.txt > "${bam.simpleName}.hg38.bam" && \
    samtools index "${bam.simpleName}.hg38.bam" && \
    rm -f header.txt body.txt

    zcat ${vcf} | awk 'BEGIN{OFS="\\t"} {if (\$1=="chr10:94757681-94855547") {\$1="chr10"; \$2+=94757680} print}' - > "${bam.simpleName}.hg38.vcf"
    zcat ${gvcf} | awk 'BEGIN{OFS="\\t"} {if (\$1=="chr10:94757681-94855547") {\$1="chr10"; \$2+=94757680} print}' - > "${bam.simpleName}.hg38.gvcf"

    echo "Completed liftover for BAM: ${bam.simpleName}.hg38.bam VCF: ${bam.simpleName}.hg38.vcf GVCF: ${bam.simpleName}.hg38.gvcf"
    """
}

process rename_contig {
    publishDir "${params.output}/${vcf_hg38.simpleName}"

    input:
    path vcf_hg38

    output:
    path("${vcf_hg38.simpleName}.vcf")

    script:
    """
    sed -E 's/(##contig=<ID=chr[^,:]+):[^,]+(,length=[0-9]+>)/\1\2/' ${vcf_hg38} > ${vcf_hg38.simpleName}.vcf
    """
}

process calculating_bamCoverage {
    publishDir "${params.output}/${bam_hg38.simpleName}"
    input:
    path bam_hg38
    path bed

    output:
    path "${bam_hg38.simpleName}.depth.txt"

    script:
    """
    bedtools coverage -a ${bed} -b ${bam_hg38} -d > "${bam_hg38.simpleName}.depth.txt"
    """
}

process creating_igvReport {
    publishDir "${params.output}/${bam_hg38.simpleName}"
    input:
    path bed
    path vcf_hg38
    path bam_hg38
    path bam_hg38_index

    output:
    path "${bam_hg38.simpleName}.igv.html"

    script:
    """
    create_report ${bed} \
        --genome hg38 \
        --flanking 1000 \
        --tracks ${vcf_hg38} ${bam_hg38} \
        --output "${bam_hg38.simpleName}.igv.html"
    """
}

process running_pharmcat {
    publishDir "${params.output}/${vcf_hg38.simpleName}"
    input:
    path vcf_hg38

    output:
    path "pharmcat_out/*"

    script:
    """
    mv ${vcf_hg38} ${vcf_hg38.baseName}.2.vcf
    cp ${vcf_hg38.baseName}.2.vcf ${vcf_hg38.baseName}.vcf
    rm ${vcf_hg38.baseName}.2.vcf
    pharmcat_pipeline ${vcf_hg38} -o pharmcat_out
    """
    //docker run --rm -v ./:/pharmcat/data pgkb/pharmcat ./pharmcat_pipeline /pharmcat/data/${vcf_hg38} -o /pharmcat/data/pharmcat_out
}

// Workflow definition
workflow {
    // Define channels
    ref_ch = Channel.fromPath(params.ref).first()
    model_ch = Channel.fromPath(params.model).first()
    bed_ch = Channel.fromPath(params.bed).first()

    // Alignment process based on file type
    bam_ch = alignmentProcess(input_ch, ref_ch)

    // Split BAM and BAI files
    bam_file_ch = bam_ch.map { it[0] }
    bam_index_ch = bam_ch.map { it[1] }

    // Index reference genome
    ref_index_ch = index(ref_ch)

    // Perform variant calling
    variant_calling_ch = variantCalling(bam_file_ch, bam_index_ch, ref_ch, ref_index_ch, model_ch)

    // Extract VCF and GVCF files from variant calling output
    vcf_file_ch = variant_calling_ch.map { it -> it.findAll { file -> file.name.endsWith('merge_output.vcf.gz') } }
    gvcf_file_ch = variant_calling_ch.map { it -> it.findAll { file -> file.name.endsWith('merge_output.gvcf.gz') } }

    // Liftover to hg38
    hg38_files_ch = liftoverToHg38(bam_file_ch, vcf_file_ch, gvcf_file_ch)
    bam_hg38_ch = hg38_files_ch.map { it -> it.findAll { file -> file.name.endsWith('.hg38.bam') } }
    vcf_hg38_ch = hg38_files_ch.map { it -> it.findAll { file -> file.name.endsWith('.hg38.vcf') } }
    bam_hg38_index_ch = hg38_files_ch.map { it -> it.findAll { file -> file.name.endsWith('.hg38.bam.bai') } }
    // gvcf_hg38_ch = hg38_files_ch.map { it -> it.findAll { file -> file.name.endsWith('.hg38.gvcf') } }

    //Rename contig
    rename_contig(vcf_hg38_ch)

    // Calculate coverage depth from BAM file
    calculating_bamCoverage(bam_hg38_ch, bed_ch)

    // Generating IGV report
    creating_igvReport(bed_ch, vcf_hg38_ch, bam_hg38_ch, bam_hg38_index_ch)

    // Generate pharmcat report
    running_pharmcat(vcf_hg38_ch)
}
