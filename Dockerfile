# Stage 1: Build environment with continuumio/miniconda3
FROM continuumio/miniconda3:latest as builder

# Install necessary tools and dependencies
RUN apt-get update && apt-get install -y wget unzip && \
    pip install --upgrade pip

# Set up Conda environment
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --add channels r && \
    conda create -n pl-cyp2c19 -c bioconda clair3 samtools minimap2 bedtools python=3.9.0 -y

# Activate Conda environment
RUN echo "conda activate pl-cyp2c19" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# Install additional Python packages
RUN pip install igv-reports

# Install pysam
RUN conda install -n pl-cyp2c19 -c bioconda pysam -y

# Stage 2: Final image based on pgkb/pharmcat:latest
FROM pgkb/pharmcat:latest as pharmcat_image