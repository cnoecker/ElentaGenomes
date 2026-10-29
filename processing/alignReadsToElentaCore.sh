#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --ntasks=6
#SBATCH --mem=24g

## Script to align MGS reads to E lenta core genes

FILE=im_samples_list.txt
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $FILE | cut -f1)
echo $SAMPLE

#SAMPLE=xxx
R1file=immunomicrobiomeMGS/${SAMPLE}_L1_1.fastq.gz
R2file=immunomicrobiomeMGS/${SAMPLE}_L1_2.fastq.gz

hisat_path=../../sftwr/hisat2-2.2.1/
module load python3
module load samtools

outdir=im_align_out/
mkdir -p $outdir

## $hisat_path/hisat2-build 
$hisat_path/hisat2 -x panaroo_out_fixed/hisat_db/Elenta_core_genes -1 $R1file -2 $R2file  \
  --no-spliced-alignment -p 6 --no-unal --sensitive > $outdir/${SAMPLE}_ElentaCore.sam
  
samtools view -bS $outdir/${SAMPLE}_ElentaCore.sam | samtools sort > $outdir/${SAMPLE}_ElentaCore.bam
samtools index $outdir/${SAMPLE}_ElentaCore.bam
samtools coverage $outdir/${SAMPLE}_ElentaCore.bam  -o $outdir/${SAMPLE}_ElentaCore.cov
samtools mpileup $outdir/${SAMPLE}_ElentaCore.bam > $outdir/${SAMPLE}_pileup_out.txt
rm $outdir/${SAMPLE}_ElentaCore.sam

