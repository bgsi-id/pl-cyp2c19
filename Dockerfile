# Stage 1: Build environment with continuumio/miniconda3
FROM continuumio/miniconda3:latest as image1

# Install necessary tools and dependencies
RUN apt-get update && apt-get install -y wget unzip && \
    pip install --upgrade pip

# Set up Conda environment
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --add channels r && \
    conda create -n pl-cyp2c19 -c bioconda clair3 samtools minimap2 bedtools python=3.9.0 -y

# Activate Conda environment and install additional packages
SHELL ["conda", "run", "-n", "pl-cyp2c19", "/bin/bash", "-c"]

# Install additional Python packages
RUN pip install igv-reports

# Install pysam
RUN conda install -n pl-cyp2c19 -c bioconda pysam -y

# Set environment variables for Conda
ENV PATH="/opt/conda/envs/pl-cyp2c19/bin:/opt/conda/bin:$PATH"
ENV CONDA_DEFAULT_ENV=pl-cyp2c19

# Ensure conda environment is activated when running commands
RUN echo "conda activate pl-cyp2c19" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# Default command
CMD ["bash"]

# Stage 2: Final image based on pgkb/pharmcat:latest
FROM pgkb/pharmcat:latest as pharmcat_image
