#!/bin/bash
# Working directory
WORK_DIR=~/variant-calling/data

# Create the directory and cd into it
cd ${WORK_DIR}

########################################################################################################################
# Requirements:
#   Java (version 8)
#   GATK (version 3.3)
########################################################################################################################

java -version
java -jar ${GATK} --help
java -jar ${PICARD}


###################################################
### Overall workflow of variant-calling with GATK #
###################################################

# GATK overall workflow:
# 1. Prepare GATK input data (pre-processing)
# 2. Local realignment around indels
# 3. Base quality recalibration
# 4. Call variants
# 5. Recalibrate variants
# 6. Post-processing


#####################
# Extract resources #
#####################
# Download resources

# Known indels
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz.tbi

# Known SNPs
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz.tbi


# Choose variable names
KNOWN_INDELS=${WORK_DIR}/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz
KNOWN_SNP=${WORK_DIR}/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz
REF_GENOME=Homo_sapiens.Chr20.fa


#######################################
### Prepare reference genome for GATK #
#######################################

# Index reference
# Command: samtools faidx
# Input: reference genome (.fa)
# Output: indexed genome (.fa.fai) (no need to mention it)
samtools faidx ${REF_GENOME}

# Create a sequence dictionnary for the reference genome
# Command: CreateSequenceDictionary (PICARDtools)
# Input: reference genome (.fa)
# Sequence dictionnary (.dict)
java -jar ${PICARD} CreateSequenceDictionary \
	R=${REF_GENOME} \
	O=${REF_GENOME/.fa/.dict}


#############################
### Prepare GATK input data #
#############################

for FILE_NAME in "HG02024" "HG02025" "HG02026"
do

# Mark Duplicate reads
# Command: MarkDuplicates (PICARDtools)
# Input: alignment (.bam)
# Ouput: alignment with duplicates marked (.bam)
java -jar ${PICARD} MarkDuplicates \
	REMOVE_DUPLICATES=FALSE \
	I=${FILE_NAME}.bam \
	O=${FILE_NAME}.marked_dups.bam \
	M=${FILE_NAME}.dup_metrics.txt

# Make sure that index already obtained. Else, do it now
# Command: samtools index / BuilBamIndex (PICARDtools)
# Input: alignment (.bam)
java -jar ${PICARD} BuildBamIndex \
	INPUT=${FILE_NAME}.marked_dups.bam


#####################################
### Local realignemnt around indels #
#####################################

# Find regions that need to be realigned
# Command: gatk RealignerTargetCreator
# Input: preprocessed alignment (.bam) + compressed known indels (.vcf.gz) + reference genome (.fa)
# Output: list of intervals (.list / .txt)
#gatk RealignerTargetCreator \
java -jar ${GATK} -T RealignerTargetCreator \
	-R ${REF_GENOME} \
	-I ${FILE_NAME}.marked_dups.bam \
	--known ${KNOWN_INDELS} \
	-o ${FILE_NAME}.target_intervals.list 

# Perform local realignement
# Command: gatk IndelRealigner
# Input: preprocessed alignment (.bam) + compressed known indels (.vcf.gz) + list of intervals (.list / .txt) + reference genome (.fa)
# Output: realigned alignment (.bam)
#gatk IndelRealigner \
java -jar ${GATK} -T IndelRealigner \
  -R ${REF_GENOME} \
	-known ${KNOWN_INDELS} \
	-I ${FILE_NAME}.marked_dups.bam \
	-targetIntervals ${FILE_NAME}.target_intervals.list \
	-o ${FILE_NAME}.realigned_reads.bam

################################
### Base quality recalibration #
################################

# Generate a base recalibration table to analyse patterns of covariation in the dataset
# Command: gatk BaseRecalibrator
# Input: realigned alignment (.bam) + compressed known indels (.vcf.gz) + compressed known snps (.vcf.gz) + reference genome (.fa)
# Output: base recalibration table (.table / .txt)
java -jar ${GATK} -T BaseRecalibrator \
	-R ${REF_GENOME} \
	-knownSites ${KNOWN_SNP} \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov ContextCovariate \
	-I ${FILE_NAME}.realigned_reads.bam \
	-o ${FILE_NAME}.recal_data.table

# Perform base quality recalibration
# Command: gatk PrintReads + BQSR option
# Input: realigned alignment (.bam) + reference genome (.fa) + base quality recalbration table (.table / .txt)
# Output: base quality recalibrated alignement (.bam)
java -jar ${GATK} -T PrintReads \
	-R ${REF_GENOME} \
	-I ${FILE_NAME}.realigned_reads.bam \
	-BQSR ${FILE_NAME}.recal_data.table \
	-o ${FILE_NAME}.recal_reads.bam


###################
### Call variants #
###################

# Perform variant calling
# Command: gatk HaplotypeCaller
# Input: base quality recalibrated alignement (.bam) + reference genome (.fa)
# Output: Genomic variant calling file (.g.vcf)
java -jar ${GATK} -T HaplotypeCaller \
                  -R ${REF_GENOME} \
                  -I ${FILE_NAME}.recal_reads.bam \
                  -o ${FILE_NAME}.g.vcf  \
                  --genotyping_mode DISCOVERY \
                  -variant_index_type LINEAR \
                  -variant_index_parameter 128000 \
                  --emitRefConfidence GVCF
done

# Perform variant calling
# Command: gatk GenotypeGVCFs
# Input : genomic variant calling files (.g.vcf) + reference genome (.fa)
# Output: Variant calling file (.vcf)
java -jar ${GATK} \
   -T GenotypeGVCFs \
   -R ${REF_GENOME} \
   --variant ${FILE_NAME}.g.vcf \
   -o ${FILE_NAME}.vcf
 