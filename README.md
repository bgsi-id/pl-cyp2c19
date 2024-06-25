# pl-cyp2c19

## Introduction

**pl-cyp2c19** is a workflow for detecting variants and generating pharmacogenomics (PGx) report for a nanopore CYP2C19 gene panel.

The pipeline can currently perform the following:

- Align nanopore reads to CYP2C19 reference gene 
- Variant calling using Clair3
- Adjust coordinate to genome build 38
- Check for missing important sites in CYP2C19
- Generate variant report using IGV Report
- Generate pharmacogenomics report using PharmCAT

## Installation

Download source and reference

```
wget https://github.com/bgsi-id/pl-cyp2c19/archive/refs/heads/dev.zip
unzip dev.zip
rm dev.zip
```

Install docker and nextflow

```
bash pl-cyp2c19-dev/setup.sh
```

Download Clair3 Model

```
## see 'https://github.com/nanoporetech/rerio/tree/master/clair3_models' for full model
wget https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_v430.tar.gz
tar -xvf r1041_e82_400bps_hac_v430.tar.gz
rm r1041_e82_400bps_hac_v430.tar.gz
```

## Usage

Input Directory Structure

The following is the structure of the `input_directory` directory:

```
example_input_directory/
├── sample1/
│   ├── sample1_file1.fastq
│   └── sample1_file2.fastq
├── sample2/
│   ├── sample2_file1.fastq.gz
│   └── sample2_file2.fastq.gz
└── sample3/
    ├── sample3_file1.bam
    └── sample3_file2.bam
```

Then, run the script in command line:

```
nextflow run bgsi-id/pl-cyp2c19 \ 
    -r 1.0.0 \
    -input_dir path/to/input_directory \ 
    --ref path/to/reference.fa \
    --threads 12 \
    --model path/to/clair3_model \
    --bed path/to/region.bed  \
    --output output_folder

# FOR EXAMPLE
nextflow run bgsi-id/pl-cyp2c19 \ 
    -r 1.0.0 \
    --input_dir example_input_directory \
    --ref pl-cyp2c19-dev/static/GRCh38.cyp2c19.fa \
    --threads 12 \
    --model r1041_e82_400bps_hac_v430 \
    --bed pl-cyp2c19-dev/static/region.bed  \
    --output pl-cyp2c19_output
```

## Disclaimer

> [!WARNING]
> Unless otherwise specified, the repository utilizes the MIT license, and code is provided "as is" without warranty of any kind, express or implied.
> The user recognizes they are using the pipeline at their own risk.


