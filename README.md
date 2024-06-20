# pl-cyp2c19
## Introduction
Align nanopore reads to reference (CYP2C19 gene, from genome build 38), variant calling using Clair3, and adjust to make it compatible with Pharmcat, check for missing important sites in CYP2C19.

## Prerequisites Installation
Install conda
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
conda init
conda activate
```

Install docker
```
sudo apt update
sudo apt install apt-transport-https ca-certificates curl software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu focal stable"
apt-cache policy docker-ce
sudo apt install docker-ce
sudo usermod -aG docker ${USER}
sudo chown ${USER}:${USER} /var/run/docker.sock
#check docker permission
docker image list
```


## Installation
Download Source code
```
wget https://github.com/bgsi-id/pl-cyp2c19/archive/refs/heads/dev.zip
unzip dev.zip
```

Install Dependencies
```
conda config --add channels defaults && \
conda config --add channels bioconda && \
conda config --add channels conda-forge && \
conda create -n pl-cyp2c19 -c bioconda clair3 samtools minimap2 bedtools python=3.9.0 -y
conda activate pl-cyp2c19
```

Download Clair3 Model
```
## see 'https://github.com/nanoporetech/rerio/tree/master/clair3_models' for full model
wget https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_v430.tar.gz
tar -xvf r1041_e82_400bps_hac_v410.tar.gz
```

Pull pharmcat docker image
```
docker pull pgkb/pharmcat
```

## Usage
Run All Sample
```
bash pl-cyp2c19/run_analysis.sh \
    /dir/to/folder \
    /dir/output
```

Run Single Sample
```
bash pl-cyp2c19/run.sh \
    -s sample \
    -f /dir/to/fastq \
    -r /reference.fa \
    -b /region.bed
    -t 48 \
    -m /dir/to/model \
    -o /dir/output

## example
bash pl-cyp2c19/run.sh \
    -s barcode06 \
    -f /data/rspon/Batch_1/fastq_pass/barcode06 \
    -r pl-cyp2c19/static/GRCh38.cyp2c19.fa \
    -b pl-cyp2c19/static/region.bed
    -t 48 \
    -m r1041_e82_400bps_sup_v410 \
    -o output-barcode06
```

Run Multiple Samples

Input Directory Structure

The following is the structure of the `input_directory` directory:
input_directory/
├── sample1/
│   ├── sample1_file1.fastq
│   └── sample1_file2.fastq
├── sample2/
│   ├── sample2_file1.fastq.gz
│   └── sample2_file2.fastq.gz
└── sample3/
    ├── sample3_file1.fastq
    └── sample3_file2.fastq.gz

Then, run the script in command line:
```
bash pl-cyp2c19/analysis_run.sh input_folder output_folder
```


Run Pharmcat Manually
```
docker run -it -v `pwd`:/data pgkb/pharmcat /pharmcat/pharmcat_pipeline /data/sample.hg38.vcf -o /data/pcat

## example
docker run -it -v `pwd`:/data pgkb/pharmcat /pharmcat/pharmcat_pipeline /data/output-barcode06/barcode06.hg38.vcf -o /data/pcat
```



