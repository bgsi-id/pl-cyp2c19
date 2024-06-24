# pl-cyp2c19
## Introduction
Align nanopore reads to reference (CYP2C19 gene, from genome build 38), variant calling using Clair3, and adjust to make it compatible with Pharmcat, check for missing important sites in CYP2C19.

## Installation
Download Source code
```
wget https://github.com/bgsi-id/pl-cyp2c19/archive/refs/heads/RT.zip
unzip RT.zip
```

Install docker and nextflow
```
bash pl-cyp2c19-RT/prequisites_installation.sh
```

Download Clair3 Model
```
## see 'https://github.com/nanoporetech/rerio/tree/master/clair3_models' for full model
wget https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_v430.tar.gz
tar -xvf r1041_e82_400bps_hac_v430.tar.gz
```

Build docker image
```
docker build -f Dockerfile --target image1 -t pl-cyp2c19-image:latest .
docker build -f Dockerfile --target pharmcat_image -t pgkb/pharmcat:latest .
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
    ├── sample3_file1.fastq
    └── sample3_file2.fastq.gz
```

Then, run the script in command line:
```
nextflow run pl-cyp2c19-RT/pl-cyp2c19.nf --input_dir path/to/input_directory --ref path/to/reference.fa --threads 12 --model path/to/clair3_model --bed path/to/region.bed  --output output_folder

# For example:
nextflow run pl-cyp2c19-RT/pl-cyp2c19.nf --input_dir example_input_directory --ref pl-cyp2c19/static/GRCh38.cyp2c19.fa --threads 12 --model r1041_e82_400bps_hac_v430 --bed pl-cyp2c19/static/region.bed  --output pl-cyp2c19_output
```



