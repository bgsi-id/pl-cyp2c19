#!/bin/bash

# Check if the input and output directories are provided
if [ -z "$1" ] || [ -z "$2" ]; then
  echo "Usage: $0 <input_folder> <output_folder>"
  exit 1
fi

# Get absolute paths of the input and output directories
INPUT_ROOT=$(readlink -f "$1")
OUTPUT_ROOT=$(readlink -f "$2")

# Create the output root directory if it doesn't exist
mkdir -p "$OUTPUT_ROOT"

# Activate the Conda environment
source /home/rtriwijaya/miniforge3/etc/profile.d/conda.sh
conda activate pl-cyp2c19

# Iterate over each folder in the input root directory
for INPUT_FOLDER in "$INPUT_ROOT"/*; do
  if [ -d "$INPUT_FOLDER" ]; then
    INPUT_FOLDER_NAME=$(basename "$INPUT_FOLDER")
    OUTPUT_FOLDER="$OUTPUT_ROOT/$INPUT_FOLDER_NAME"
    
    # Create the output directory for the current input folder
    mkdir -p "$OUTPUT_FOLDER"
    
    # Run the command
    bash pl-cyp2c19/run.sh \
      -s "$INPUT_FOLDER_NAME" \
      -f "$INPUT_FOLDER" \
      -r pl-cyp2c19/static/GRCh38.cyp2c19.fa \
      -b pl-cyp2c19/static/region.bed \
      -m reference/r1041_e82_400bps_hac_v430 \
      -t 48 \
      -o "$OUTPUT_FOLDER"
  fi
done
