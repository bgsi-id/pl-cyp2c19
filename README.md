# pgkb-hwgs

## Introduction

**pgkb-hwgs** is an alternative guideline to generate pharmacogenomics report from a targeted CYP2C19 gene panel VCF.

The pipeline can currently perform the following:

- Phase VCF based on reference
- Generate pharmacogenomics report using PharmCAT

## Steps

- Separate multiallelic variants to different rows in VCF using bcftools.
- Phase the input VCF with filtered 1000Genomes VCF as reference using Eagle. The reference VCF is obtained from (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) and filtered on region chr10:94761288-9485291.
- Run the PharmCAT installed in docker.
  
## Prerequisites
Install docker.

Install conda.

## Installation

Download and install Eagle_v2.4.1 (https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz)

Install PharmCAT with docker (https://pharmcat.org/using/PharmCAT-in-Docker/)

Install bcftools in a Conda environment

## Commands

Input Directory Structure

The following is the structure of the `input_directory` directory:

```
example_input_directory/
├── sample1.vcf.gz
├── sample2.vcf.gz
├── sample3.vcf.gz
└──...
```

1. Declare input and activate environment

```
INPUT_DIR="/home/$USER/input_vcf_dir"
VCF_REF="$HOME/Eagle_v2.4.1/ALL_chr10_phase3_cyp2c19.vcf.gz"
GENETIC_MAP="$HOME/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
OUT_DIR="$HOME/output_pharmcat_hwgs"

# Set up memory usage limit 
ulimit -v 8097152

# Create the output directory
mkdir -p $OUT_DIR

# Activate environment installed with bcftools
conda activate bcftools
```

2. Phasing and pharmcat

```
# Run the bcftools, Eagle, and pharmCAT

for i in $(ls $INPUT_DIR/*vcf.gz); do \
    id=$(basename $i .vcf.gz) ;\
    mkdir -p "$OUT_DIR"/"$id" ;\
    #handle multiallelic rows in vcf
    bcftools norm -m -any $i -o $OUT_DIR/"$id"/"$id"_mono.vcf.gz ;\
    
    #index input vcf
    tabix $OUT_DIR/"$id"/"$id"_mono.vcf.gz ;\
    
    #phasing input vcf based on available reference
    $HOME/Eagle_v2.4.1/eagle --vcfTarget $OUT_DIR/"$id"/"$id"_mono.vcf.gz \
        --outPrefix="$OUT_DIR"/"$id"/"$id" \
        --geneticMapFile $GENETIC_MAP \
        --chrom chr10 \
        --vcfRef $VCF_REF \
        --numThreads 8 \
        --outputUnphased ;\
    
    mv "$OUT_DIR"/"$id"/"$id".unphased.vcf.gz "$OUT_DIR"/"$id"/"$id".vcf.gz;\
    rm $OUT_DIR/"$id"/*_mono.vcf.gz*  ;\

    #running pharmcat
    docker run --rm -v $OUT_DIR/"$id":/pharmcat/data pgkb/pharmcat ./pharmcat_pipeline /pharmcat/data/"$id".vcf.gz -o /pharmcat/data/"$id"/ ;\
    rm $OUT_DIR/*.vcf.gz ;\
    done
```



## Disclaimer

> [!WARNING]
> Unless otherwise specified, the repository utilizes the MIT license, and code is provided "as is" without warranty of any kind, express or implied.
> The user recognizes they are using the pipeline at their own risk.


