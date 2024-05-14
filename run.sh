while getopts ":s:t:f:r:m:o:" opt; do
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
required_tools=("minimap2" "samtools" "run_clair3.sh")

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
zcat $FASTQ/*.fastq.gz | minimap2 -ax map-ont $REF - | samtools sort - > $OUTPUT/$SAMPLE.bam && samtools index $OUTPUT/$SAMPLE.bam

# Check the exit status of the previous command
if [ $? -ne 0 ]; then
    echo "Error in mapping reads to reference. Exiting."
    exit 1
fi

# Step 2
echo "(2/4) Variant calling"
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
echo "(3/4) Liftover to hg38"

rm -f header.txt body.txt && \
samtools view -H $OUTPUT/$SAMPLE.bam | sed 's/SN:chr10:94761900-94855547/SN:chr10/' > header.txt && \
samtools view $OUTPUT/$SAMPLE.bam | sed 's/chr10:94761900-94855547/chr10/' - | \
awk -v shift=94761899 '$3 == "chr10" {OFS="\t"; $4 = $4 + shift} {print $0}' > body.txt && \
cat body.txt >> header.txt && \
samtools view -b header.txt > $OUTPUT/$SAMPLE.hg38.bam && \
samtools index $OUTPUT/$SAMPLE.hg38.bam && \
rm -f header.txt body.txt

zcat $OUTPUT/clair3/merge_output.vcf.gz | awk 'BEGIN{OFS="\t"} {if ($1=="chr10:94761900-94855547") {$1="chr10"; $2+=94761899} print}' - > $OUTPUT/$SAMPLE.hg38.vcf
zcat $OUTPUT/clair3/merge_output.gvcf.gz | awk 'BEGIN{OFS="\t"} {if ($1=="chr10:94761900-94855547") {$1="chr10"; $2+=94761899} print}' - > $OUTPUT/$SAMPLE.hg38.gvcf

# Check the exit status of the previous command
if [ $? -ne 0 ]; then
    echo "Error in liftover to hg38. Exiting."
    exit 1
fi

# Step 4
echo "(4/4) Checking variants"
echo "CHR\tPOS\tDP" > $OUTPUT/$SAMPLE.depth.txt
# Loop through the provided positions
for position in \
  "94761900" \
  "94761288" \
  "94762656" \
  "94762693" \
  "94762706" \
  "94762712" \
  "94762715" \
  "94762755" \
  "94762760" \
  "94762788" \
  "94762856" \
  "94775106" \
  "94775121" \
  "94775160" \
  "94775185" \
  "94775367" \
  "94775416" \
  "94775423" \
  "94775453" \
  "94775489" \
  "94775507" \
  "94780574" \
  "94780579" \
  "94780653" \
  "94781858" \
  "94781859" \
  "94781944" \
  "94781999" \
  "94842861" \
  "94842866" \
  "94842879" \
  "94842995" \
  "94849995" \
  "94852738" \
  "94852765" \
  "94852785" \
  "94852914"

do
  # Extract depth from VCF for the given position
  depth=$(grep $position $OUTPUT/$SAMPLE.hg38.gvcf | awk '{print $10}' | awk -F':' '{print $3}')
  #Print result
  if [ -n "$depth" ]; then
    echo -e "chr10\t$position\t$depth" >> $OUTPUT/$SAMPLE.depth.txt
  else
    echo -e "chr10\t$position\tNo depth" >> $OUTPUT/$SAMPLE.depth.txt
  fi
done

if [ $? -ne 0 ]; then
    echo "Error in liftover to hg38. Exiting."
    exit 1
fi

echo "Complete!"