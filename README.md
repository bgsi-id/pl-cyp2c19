# pl-cyp2c19
## Introduction
Align nanopore reads to reference (CYP2C19 gene, from genome build 38), variant calling using Clair3, and adjust to make it compatible with Pharmcat, check for missing important sites in CYP2C19.

## Installation
Install Dependencies
```
conda config --add channels defaults && \
conda config --add channels bioconda && \
conda config --add channels conda-forge && \
conda create -n pl-cyp2c19 -c bioconda clair3 samtools minimap2 python=3.9.0 -y
conda activate pl-cyp2c19
```

Download Reference
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip -d GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna chr10:94757681-94855547 > GRCh38.cyp2c19.fa 
samtools faidx GRCh38.cyp2c19.fa
```

Download Clair3 Model
```
## see 'https://github.com/nanoporetech/rerio/tree/master/clair3_models' for full model
wget https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_v430.tar.gz
tar -xvf r1041_e82_400bps_sup_v410.tar.gz
```
## Usage
Run Pipeline
```
bash run.sh \
    -s sample \
    -f /dir/to/fastq \
    -r /reference.fa \
    -t 48 \
    -m /dir/to/model \
    -o /dir/output

## example
# bash run.sh \
#     -s barcode06 \
#     -f /data/rspon/Batch_1/fastq_pass/barcode06 \
#     -r /data/refs/GRCh38.cyp2c19.fa \
#     -t 48 \
#     -m /data/models/r1041_e82_400bps_sup_v410 \
#     -o output-barcode06
```

Run Pharmcat
```
docker run -it -v `pwd`:/data pgkb/pharmcat /pharmcat/pharmcat_pipeline /data/sample.hg38.vcf -o /data/pcat

## example
# docker run -it -v `pwd`:/data pgkb/pharmcat /pharmcat/pharmcat_pipeline /data/output-barcode06/barcode06.hg38.vcf -o /data/pcat
```



