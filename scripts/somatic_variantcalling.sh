#!/bin/bash
# ==============================================================
# SOMATIC VARIANT CALLING PIPELINE (Tumor vs Normal)
# Author: Karthik R
# Samples: SRR27718988 (Tumor) and SRR27718966 (Normal)
# ==============================================================

# Use 8 threads for faster performance
THREADS=8

# ==============================================================
# STEP 1: QUALITY CONTROL
# ==============================================================
echo "[Step 1] Running FastQC..."
fastqc SRR27718988_1.fastq.gz SRR27718988_2.fastq.gz SRR27718966_1.fastq.gz SRR27718966_2.fastq.gz -t $THREADS -o QC_reports

echo "[Step 2] Generating MultiQC report..."
multiqc QC_reports -o QC_reports

# ==============================================================
# STEP 2: REFERENCE INDEXING
# ==============================================================
echo "[Step 3] Indexing reference genome..."
bwa index reference/GRCh38.p14.fa
samtools faidx reference/GRCh38.p14.fa

# ==============================================================
# STEP 3: ALIGNMENT
# ==============================================================
echo "[Step 4] Aligning reads..."
bwa mem -t $THREADS reference/GRCh38.p14.fa SRR27718988_1.fastq.gz SRR27718988_2.fastq.gz | \
samtools sort -@ $THREADS -o SRR27718988.sorted.bam
samtools index SRR27718988.sorted.bam

bwa mem -t $THREADS reference/GRCh38.p14.fa SRR27718966_1.fastq.gz SRR27718966_2.fastq.gz | \
samtools sort -@ $THREADS -o SRR27718966.sorted.bam
samtools index SRR27718966.sorted.bam

# ==============================================================
# STEP 4: ADD READ GROUPS
# ==============================================================
echo "[Step 5] Adding read groups..."
gatk AddOrReplaceReadGroups \
  -I SRR27718988.sorted.bam \
  -O SRR27718988.RG.bam \
  -RGID SRR27718988 -RGLB Lib1 -RGPL ILLUMINA -RGPU Unit1 -RGSM SRR27718988

gatk AddOrReplaceReadGroups \
  -I SRR27718966.sorted.bam \
  -O SRR27718966.RG.bam \
  -RGID SRR27718966 -RGLB Lib1 -RGPL ILLUMINA -RGPU Unit1 -RGSM SRR27718966

# ==============================================================
# STEP 5: MARK DUPLICATES
# ==============================================================
echo "[Step 6] Removing duplicates..."
gatk MarkDuplicates \
  -I SRR27718988.RG.bam \
  -O SRR27718988.markdup.bam \
  -M SRR27718988.metrics.txt \
  --CREATE_INDEX true

gatk MarkDuplicates \
  -I SRR27718966.RG.bam \
  -O SRR27718966.markdup.bam \
  -M SRR27718966.metrics.txt \
  --CREATE_INDEX true

# ==============================================================
# STEP 6: SOMATIC VARIANT CALLING
# ==============================================================
echo "[Step 7] Running Mutect2..."
gatk Mutect2 \
  -R reference/GRCh38.p14.fa \
  -I SRR27718988.markdup.bam -tumor SRR27718988 \
  -I SRR27718966.markdup.bam -normal SRR27718966 \
  -O SRR27718988_vs_SRR27718966.vcf.gz

# ==============================================================
# STEP 7: FILTER VARIANTS
# ==============================================================
echo "[Step 8] Filtering variants..."
gatk FilterMutectCalls \
  -V SRR27718988_vs_SRR27718966.vcf.gz \
  -O SRR27718988_vs_SRR27718966.filtered.vcf.gz

# ==============================================================
# STEP 8: ANNOTATION
# ==============================================================
echo "[Step 9] Annotating variants with VEP..."
vep -i SRR27718988_vs_SRR27718966.filtered.vcf.gz \
    -o SRR27718988_vs_SRR27718966.annotated.vep.txt \
    --assembly GRCh38 --everything --fork $THREADS

# ==============================================================
# COMPLETION MESSAGE
# ==============================================================
echo "SOMATIC VARIANT CALLING COMPLETED SUCCESSFULLY!"
echo "Outputs:"
echo "  - QC Reports: QC_reports/"
echo "  - BAM Files: *.bam"
echo "  - Raw VCF: SRR27718988_vs_SRR27718966.vcf.gz"
echo "  - Filtered VCF: SRR27718988_vs_SRR27718966.filtered.vcf.gz"
echo "  - Annotation: SRR27718988_vs_SRR27718966.annotated.vep.txt"
