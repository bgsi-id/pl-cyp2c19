while getopts ":s:t:f:r:b:m:o:" opt; do
  case $opt in
    s)
      SAMPLE=$OPTARG
      ;;
    t)
      T=$OPTARG
      ;;
    f)
      FASTQ=$OPTARG
      ;;
    r)
      REF=$OPTARG
      ;;
    b)
      BED=$OPTARG
      ;;
    m)
      MODEL=$OPTARG
      ;;
    o)
      OUTPUT=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

echo "Checking environment ..."
required_tools=("minimap2" "samtools" "run_clair3.sh" "bedtools" "create_report")

for tool in "${required_tools[@]}"; do
    which "$tool" > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        echo "- $tool is already installed."
    else
        echo "- $tool is not installed."
        exit
    fi
done

mkdir $OUTPUT

# Step 1
echo "(1/4) Mapping reads to reference"
if ls $FASTQ/*.fastq.gz 1> /dev/null 2>&1; then
    # Handle fastq.gz files
    zcat $FASTQ/*.fastq.gz | minimap2 -ax map-ont $REF - | samtools sort -o $OUTPUT/$SAMPLE.bam - && samtools index $OUTPUT/$SAMPLE.bam
elif ls $FASTQ/*.fastq 1> /dev/null 2>&1; then
    # Handle fastq files
    samtools fastq $FASTQ/*.fastq | minimap2 -ax map-ont $REF - | samtools sort -o $OUTPUT/$SAMPLE.bam - && samtools index $OUTPUT/$SAMPLE.bam
#elif ls $FASTQ_DIR/*.bam 1> /dev/null 2>&1; then
#    # Handle BAM files by converting to FASTQ
#    samtools fastq $FASTQ_DIR/*.bam | minimap2 -ax map-ont $REF - | samtools sort -o $OUTPUT_DIR/$SAMPLE.bam -
else
    echo "No fastq, fastq.gz, or bam files found in $FASTQ_DIR"
exit 1
fi

# Check the exit status of the previous command
if [ $? -ne 0 ]; then
    echo "Error in mapping reads to reference. Exiting."
    exit 1
fi

# Step 2
echo "(2/5) Variant calling"
run_clair3.sh \
  --bam_fn=$OUTPUT/$SAMPLE.bam \
  --ref_fn=$REF \
  --platform="ont" \
  --threads=$T \
  --model_path=$MODEL \
  --include_all_ctgs \
  --var_pct_full=1 \
  --var_pct_phasing=1 \
  --gvcf \
  --output=$OUTPUT/clair3

# Check the exit status of the previous command
if [ $? -ne 0 ]; then
    echo "Error in variant calling. Exiting."
    exit 1
fi

# Step 3
echo "(3/5) Liftover to hg38"

rm -f header.txt body.txt && \
samtools view -H $OUTPUT/$SAMPLE.bam | sed 's/SN:chr10:94757681-94855547/SN:chr10/' > header.txt && \
samtools view $OUTPUT/$SAMPLE.bam | sed 's/chr10:94757681-94855547/chr10/' - | \
awk -v shift=94757680 '$3 == "chr10" {OFS="\t"; $4 = $4 + shift} {print $0}' > body.txt && \
cat body.txt >> header.txt && \
samtools view -b header.txt > $OUTPUT/$SAMPLE.hg38.bam && \
samtools index $OUTPUT/$SAMPLE.hg38.bam && \
rm -f header.txt body.txt

zcat $OUTPUT/clair3/merge_output.vcf.gz | awk 'BEGIN{OFS="\t"} {if ($1=="chr10:94757681-94855547") {$1="chr10"; $2+=94757680} print}' - > $OUTPUT/$SAMPLE.hg38.vcf
zcat $OUTPUT/clair3/merge_output.gvcf.gz | awk 'BEGIN{OFS="\t"} {if ($1=="chr10:94757681-94855547") {$1="chr10"; $2+=94757680} print}' - > $OUTPUT/$SAMPLE.hg38.gvcf

# Check the exit status of the previous command
if [ $? -ne 0 ]; then
    echo "Error in liftover to hg38. Exiting."
    exit 1
fi

# Step 4
echo "(4/5) Checking variants"
bedtools coverage -a $BED -b $OUTPUT/$SAMPLE.hg38.bam > $OUTPUT/$SAMPLE.depth.txt -d

echo "(5/5) Generating report"

CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
create_report $BED \
    --genome hg38 \
    --flanking 1000 \
    --tracks $OUTPUT/$SAMPLE.hg38.vcf $OUTPUT/$SAMPLE.hg38.bam \
    --output $OUTPUT/$SAMPLE.igv.html

echo "Complete!"