#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l scratch=50G
#$ -l h_rt=48:0:0
#$ -pe smp 6
#$ -t 1-4

######### Nanopore isolate assembly pipeline
######################################################################################################################################################################################
# Bacterial assembly Pipeline 
# Instructions:
# - Obtain the hybrid_assembly_manifest.tsv and modify as appropriate
# - Adjust the lines above as appropriate:
#   -t 1-10 refers to running samples 1 to 10 in the manifest, adjust to be 1 to the number of total samples in your run
#   -pe smp 12 refers to the number of cores used ($NSLOTS), increase if you want your jobs to finish faster, but be aware it may take longer to start
#   -mem_free=8G refers to the amount of RAM to be used, the actual RAM usage is the number of processors * mem_free = 8x12 = 96Gb
#   -scratch=200G refers to the amount of space required on the node running the code, 200Gb should be sufficient for most applications
#   -h_rt=48:0:0 is the maximum run time on a per sample basis, in this case 48 hours
# - Do not modify anything below this line unless you know what you are doing
######################################################################################################################################################################################


######################################################################################################################################################################################
### Setting up for Run
export PATH=$PATH:$QB3/shared_resources/bin/
#Extract Parameters from manifest one line at a time assuming user settinsg on are lines 6 to 12
# Also removes the last character (carriage return newline)
for i in $(seq 6 14); do export $(sed "${i}q;d" hybrid_assembly_manifest.tsv | sed 's/\#.*//g' | sed 's/\r//'); done 
echo "ReadDir is "$ReadDir
echo "WorkDir is "$WorkDir
echo "OutDir is "$OutDir
echo "BUSCO_DATA_PATH is "$BUSCO_DATA_PATH
echo "REF_GENOME_COMPARE is "$REF_GENOME_COMPARE
echo "MEDAKA_CONSENSUS_MODEL is "$MEDAKA_CONSENSUS_MODEL
echo "MEDAKA_VARIANT_MODEL is "$MEDAKA_VARIANT_MODEL
echo "NANOPORE_FILTER is "$NANOPORE_FILTER
echo "NANOPORE_FIRST_CUTOFF is "$NANOPORE_FIRST_CUTOFF

#FILTER=true

mkdir -p $OutDir
mkdir -p $OutDir/assembly_logs
WorkDir=$( mktemp -d -t -p /scratch $WorkDir )
JobID=$(expr $SGE_TASK_ID + 18) # offset the SGE_TASK_ID to get to the correct line in the manifest
SampleID=$(sed "${JobID}q;d" hybrid_assembly_manifest.tsv | cut -f 1)
echo "SampleID is "$SampleID
R1_L1=$(sed "${JobID}q;d" hybrid_assembly_manifest.tsv | cut -f 2 | sed 's/\r$//')
R2_L1=$(sed "${JobID}q;d" hybrid_assembly_manifest.tsv | cut -f 3 | sed 's/\r$//')
nanopore1=$(sed "${JobID}q;d" hybrid_assembly_manifest.tsv | cut -f 4 | sed 's/\r$//')
nanopore2=$(sed "${JobID}q;d" hybrid_assembly_manifest.tsv | cut -f 5 | sed 's/\r$//')
#LibNumber=$(if [ -f "$ReadDir/$R1_L2" ]; then echo 2; else echo 1; fi) # to know if there are multiple libraries
#echo "LibNumber is "$LibNumber

MEMU=96 #Amount of RAM requested (ie should be equal to pe smp x mem_free)

mkdir -p $WorkDir/reads
cp $ReadDir/$SampleID/$R1_L1 $WorkDir/reads/${SampleID}_R1_L1.fastq.gz
cp $ReadDir/$SampleID/$R2_L1 $WorkDir/reads/${SampleID}_R2_L1.fastq.gz
if [  -f "$ReadDir/$SampleID/$nanopore1" ]; then cp $ReadDir/$SampleID/$nanopore1 $WorkDir/reads/${SampleID}_nanopore1.fastq; fi
if [ -f "$ReadDir/$SampleID/$nanopore2"  ]; then cp $ReadDir/$SampleID/$nanopore2 $WorkDir/reads/${SampleID}_nanopore2.fastq; fi
if [ -f "$ReadDir/$SampleID/$nanopore2"  ]; then cat $WorkDir/reads/${SampleID}_nanopore2.fastq >> $WorkDir/reads/${SampleID}_nanopore1.fastq; fi


pigz -dc $WorkDir/reads/${SampleID}_R1_L1.fastq.gz > $WorkDir/${SampleID}_1.fastq
#if [ "$LibNumber" == "2" ]; then pigz -dc $WorkDir/reads/${SampleID}_R1_L2.fastq.gz >> $WorkDir/${SampleID}_1.fastq; fi

pigz -dc $WorkDir/reads/${SampleID}_R2_L1.fastq.gz > $WorkDir/${SampleID}_2.fastq
#if [ "$LibNumber" == "2" ]; then pigz -dc $WorkDir/reads/${SampleID}_R2_L2.fastq.gz >> $WorkDir/${SampleID}_2.fastq; fi


######################################################################################################################################################################################
### Set up run log
exec > ${OutDir}/assembly_logs/${SampleID}.log 2> ${OutDir}/assembly_logs/${SampleID}.err
echo $(date) ..... Running analysis for $SampleID on $(hostname)
echo SampleID: $SampleID
echo SGE_TASK_ID: $SGE_TASK_ID
echo Temporary Directory: $WorkDir
echo R1_L1: $R1_L1
echo R2_L1: $R2_L1
echo nanopore1: $nanopore1
echo nanopore2: $nanopore2
#echo LibNumber: $LibNumber

mkdir -p $OutDir/$SampleID

### QC 


## Make plots with Nanoplot
mkdir -p $OutDir/$SampleID/nanoplot_prefilt
conda activate Nanoplot_1.40.0
NanoPlot -t $NSLOTS --verbose -o $OutDir/$SampleID/nanoplot_prefilt --raw -f pdf --plots kde --N50 --fastq $WorkDir/reads/${SampleID}_nanopore1.fastq

if [[ $NANOPORE_FILTER == "TRUE" ]]; then
	# Use Filtlong to filter - let's try skipping this for now
	mkdir -p $OutDir/$SampleID/filtlong/

	filtlong --min_length 500 --verbose $WorkDir/reads/${SampleID}_nanopore1.fastq > $OutDir/$SampleID/filtlong/$SampleID".fastq"

	# Run Nanoplot again
	mkdir -p $OutDir/$SampleID/nanoplot_postfilt
	NanoPlot -t $NSLOTS --verbose -o $OutDir/$SampleID/nanoplot_postfilt/ --raw -f pdf --plots kde --N50 --fastq $OutDir/$SampleID/filtlong/$SampleID".fastq"

	# Update WorkDir version for subsequent steps
	cp $OutDir/$SampleID/filtlong/$SampleID".fastq" $WorkDir/reads/${SampleID}_nanopore1.fastq
	#READDIR=$OutDir/$SampleID/filtlong/
	conda deactivate
fi


## Run assembly with Flye
conda activate flye_2.9.2
#flye --nano-raw $OUTDIR/filtlong/$SAMPLEID"_filt.fastq" -g 3.6m -o $OUTDIR/flye_assembly/ -t $NSLOTS
flye --nano-corr $WorkDir/reads/${SampleID}_nanopore1.fastq -g 3.6m -o $OutDir/$SampleID/flye_assembly/ -t $NSLOTS --plasmids
conda deactivate

## Check stats from flye to decide whether to switch to Illumina-first
flye_largest_frag=$(grep "Largest frg:" $OutDir/$SampleID/flye_assembly/flye.log | cut -f 3)
echo "Largest fragment - Flye: "$flye_largest_frag

if [[ $NANOPORE_FIRST_CUTOFF -gt $flye_largest_frag ]]; then
	# Switch to Unicycler
    echo "Low-quality nanopore assembly - Switching to Unicycler"
    conda activate unicycler_0.5.0

	mkdir -p $OutDir/$SampleID/unicycler/
	unicycler -1 $WorkDir/${SampleID}_1.fastq -2 $WorkDir/${SampleID}_2.fastq -l $WorkDir/reads/${SampleID}_nanopore1.fastq -t $NSLOTS -o $OutDir/$SampleID/unicycler/
	#cp -r $WorkDir/${SampleID}"_unicycler"/* $OutDir/${SampleID}

	conda deactivate
	cp $OutDir/$SampleID/unicycler/assembly.fasta $OutDir/${SampleID}_final.fasta
    
else
	# Carry on with Flye assembly
	# # Option here to add Trycycler to get better consensus - only worthwhile for complete genomes and 25x+ coverage according to paper
# 
	# ## Run long-read polishing - Medaka. Check whether model is the correct selection for your data!!!
	conda activate medaka_1.8.0
#medaka_consensus -i $OUTDIR/filtlong/$SAMPLEID"_filt.fastq" -d $OUTDIR/flye_assembly/assembly.fasta -t $NSLOTS -o $OUTDIR/medaka
	medaka_consensus -i $WorkDir/reads/${SampleID}_nanopore1.fastq -d $OutDir/$SampleID/flye_assembly/assembly.fasta -t $NSLOTS -o $OutDir/$SampleID/medaka \
		-m $MEDAKA_CONSENSUS_MODEL -f

	conda deactivate

	# Circularize if possible - doesn't actually use spades. install circlator into base conda environment (pip3 install circlator)
	$QB3/shared_resources/Conda/bin/circlator fixstart $OutDir/$SampleID/medaka/consensus.fasta $OutDir/$SampleID/$SampleID"_fixedstart"

	# Align short reads to medaka genome
	bwa index $OutDir/$SampleID/$SampleID"_fixedstart.fasta"
	bwa mem -t $NSLOTS -a $OutDir/$SampleID/$SampleID"_fixedstart.fasta" $WorkDir/${SampleID}_1.fastq > $WorkDir/${SampleID}_shortreads_1.sam
	bwa mem -t $NSLOTS -a $OutDir/$SampleID/$SampleID"_fixedstart.fasta" $WorkDir/${SampleID}_2.fastq > $WorkDir/${SampleID}_shortreads_2.sam
	polypolish_insert_filter.py --in1 $WorkDir/${SampleID}_shortreads_1.sam --in2 $WorkDir/${SampleID}_shortreads_2.sam --out1 $WorkDir/${SampleID}_shortreads_filt1.sam --out2 $WorkDir/${SampleID}_shortreads_filt2.sam
	# Short-read polishing
	polypolish $OutDir/$SampleID/$SampleID"_fixedstart.fasta" $WorkDir/${SampleID}_shortreads_filt1.sam $WorkDir/${SampleID}_shortreads_filt2.sam > $OutDir/${SampleID}_polypolished.fasta
	
	mv $OutDir/${SampleID}_polypolished.fasta $OutDir/${SampleID}_final.fasta

fi


## Align reads to final genome
# Long reads
minimap2 -L -N 0 -a $OutDir/$SampleID"_final.fasta" $WorkDir/reads/${SampleID}_nanopore1.fastq > $OutDir/$SampleID/$SampleID"_longreads.sam"
module load samtools # Get Wynton-compatible version
	
samtools view -S -b $OutDir/$SampleID/$SampleID"_longreads.sam" | samtools sort > $OutDir/$SampleID/$SampleID"_longreads.bam"
rm $OutDir/$SampleID/$SampleID"_longreads.sam"
samtools index $OutDir/$SampleID/$SampleID"_longreads.bam" $OutDir/$SampleID/$SampleID"_longreads.bai"
samtools coverage $OutDir/$SampleID/$SampleID"_longreads.bam" -o $OutDir/$SampleID/$SampleID"_longreads.cov"

# Short reads
#Skip multimaps or get them all (-a)
hisat2-build $OutDir/$SampleID"_final.fasta" $SampleID

hisat2 -x $SampleID \
     -1 $ReadDir/$SampleID/$R1_L1 \
     -2 $ReadDir/$SampleID/$R2_L1 \
     -X 1000 -a  \
     --very-sensitive\
     --new-summary --summary-file $OutDir/$SampleID/$SampleID"_shortread_summary.tsv" \
     --threads $NSLOTS \
     -S $OutDir/$SampleID/$SampleID"_shortreads.sam"

samtools view -S -b $OutDir/$SampleID/$SampleID"_shortreads.sam" | samtools sort > $OutDir/$SampleID/$SampleID"_shortreads.bam"
rm $OutDir/$SampleID/$SampleID"_shortreads.sam"
samtools index $OutDir/$SampleID/$SampleID"_shortreads.bam" $OutDir/$SampleID/$SampleID"_shortreads.bai"
samtools coverage $OutDir/$SampleID/$SampleID"_shortreads.bam" -o $OutDir/$SampleID/$SampleID"_shortreads.cov"



# ## Annotate with prokka
# conda activate prokka_1.13.3
# prokka --cpus $NSLOTS --outdir $OUTDIR/prokka_out/ --prefix $SAMPLEID $OutDir/${SampleID}_final.fasta
# conda deactivate

#### Run Quast - QC

conda activate quast-5.0.2
quast.py $OutDir/${SampleID}_final.fasta --threads $NSLOTS -o $OutDir/$SampleID/quast/
conda deactivate

## Run BUSCO - More QC
##Busco only takes relative paths
conda activate busco_5.3.2

mkdir -p $SampleID # Have it save to current directory, then move to the correct place
busco -i $OutDir/${SampleID}_final.fasta -o $SampleID/busco_out -m genome --auto-lineage-prok --offline --download_path $BUSCO_DATA_PATH --cpu $NSLOTS
mv $SampleID/busco_out/ $OutDir/$SampleID/
rmdir $SampleID

conda deactivate

#Bakta
conda activate bakta_v1.7.0
mkdir -p $OutDir/${SampleID}/bakta_annot/

bakta --db $QB3/shared_resources/databases/bakta_db/db/ --prefix $SampleID --output $OutDir/${SampleID}/bakta_annot/ \
--prefix $SampleID $OutDir/${SampleID}_final.fasta

conda deactivate

conda activate eggnog-mapper_2.1.9

mkdir -p $OutDir/$SampleID/eggnog/

emapper.py -i $OutDir/${SampleID}/bakta_annot/${SampleID}.faa -o $SampleID --output_dir $OutDir/$SampleID/eggnog/ --tax_scope bacteria \
	--data_dir $QB3/shared_resources/databases/eggnog_mapper/ --cpu $NSLOTS
	
conda deactivate

## Interproscan
#module load openjdk
 #shared_resources/sftwrshare/interproscan-5.54-87.0/interproscan.sh -i $OutDir/${SampleID}/bakta_annot/${SampleID}.faa -pa

## Annotate taxonomy with GTDB-tk
conda activate gtdbtk_2.1.1
echo "$OutDir/${SampleID}_final.fasta	$SampleID" > $WorkDir/tmp_${SampleID}_batch
gtdbtk classify_wf --batchfile $WorkDir/tmp_${SampleID}_batch --out_dir $OutDir/$SampleID/gtdbtk \
	--cpus $NSLOTS --extension fasta

conda deactivate



