#!bin/bash
#Installation du dossier de travail
WORK_DIR =~/TP-NGS/data
mkdir -p ${WORK_DIR}
cd ${WORK_DIR}

#################################################################
## Téléchargement des et préparation des séquences nécessaires ##
#################################################################
#Téléchargement de la séquence de référence
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz -O Homo_sapiens.Chr20.fa.gz
#Extraction
gunzip Homo_sapiens.Chr20.fa.gz
REF_GENOME=Homo_sapiens.Chr20.fa
#Indexation
bwa index ${REF_GENOME}

FTP_SEQ_FOLDER=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3 # Ftp folder from 1000Genomes project
# Download index file containing sequencing runs information
wget ${FTP_SEQ_FOLDER}/20130502.phase3.analysis.sequence.index -O 20130502.phase3.index

# Référence de l'index pour GAKT
samtools faidx ${REF_GENOME}
# Dictionnaire de séquence pour la référence
java -jar ${PICARD} CreateSequenceDictionary \
	R=${REF_GENOME} \
	O=${REF_GENOME/.fa/.dict}

# Known indels
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz -o Indels.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz.tbi -o Indels.vcf.gz.tbi
# Known SNPs
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz -o SNPs.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz.tbi -o SPNs.vcf.gz.tbi
KNOWN_INDELS=${WORK_DIR}/Indels.vcf.gz
KNOWN_SNP=${WORK_DIR}/SNPs.vcf.gz

# Récupération de Pedigree
PEDIGREE=20130606_g1k.ped
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/${PEDIGREE} -O ${PEDIGREE}
# Modify PEDIGREE to keep only first 6 columns and no header
cut -f 1-6 ${PEDIGREE} |sed '1,1d' > ${PEDIGREE}.txt


##################################################################################
## Récupération, adressage et détection de variants sur les séquences d'intérêt ##
##################################################################################
for INDIVIDUAL in "HG02024" "HG02025" "HG02026" 
do
#Préparation de la liste
echo ""> ${INDIVIDUAL}.bamlist

#Sélection des runs à étudier
grep ${INDIVIDUAL} 20130502.phase3.index | grep "exome" | grep 'PAIRED' | grep "Pond-" | grep -v 'Solexa' | grep -v 'from blood' | grep -v '_1.filt.fastq.gz' | grep -v '_2.filt.fastq.gz' | sed 's/\t/,/g' > ${INDIVIDUAL}.index

#Pour chaque séquence, alignement, tri, ajout des read groups et indexation
cat ${INDIVIDUAL}.index | while IFS="," read FASTQ_FILE MD5 RUN_ID STUDY_ID STUDY_NAME CENTER_NAME SUBMISSION_ID SUBMISSION_DATE SAMPLE_ID SAMPLE_NAME POPULATION EXPERIMENT_ID INSTRUMENT_PLATFORM INSTRUMENT_MODEL LIBRARY_NAME RUN_NAME RUN_BLOCK_NAME INSERT_SIZE LIBRARY_LAYOUT PAIRED_FASTQ WITHDRAWN WITHDRAWN_DATE COMMENT READ_COUNT BASE_COUNT ANALYSIS_GROUP
do
    # Variables definition
    FASTQ_FILE_1=${FASTQ_FILE/.filt.fastq.gz/_1.filt.fastq.gz} # Path of the fasta file in the FTP folder
    FASTQ_FILE_2=${FASTQ_FILE/.filt.fastq.gz/_2.filt.fastq.gz} # Path of the fasta file in the FTP folder (pairing file)

    # Download paired sequencing reads for the ${SAMPLE_NAME}
    wget ${FTP_SEQ_FOLDER}/${FASTQ_FILE_1} -O ${SAMPLE_NAME}_${RUN_ID}_1.filt.fastq.gz
    wget ${FTP_SEQ_FOLDER}/${FASTQ_FILE_2} -O ${SAMPLE_NAME}_${RUN_ID}_2.filt.fastq.gz

    # Map, filter, and sort the paired reads of the sequencing run against the reference genome
    bwa mem -M -t 2 ${REF_GENOME}  ${SAMPLE_NAME}_${RUN_ID}_1.filt.fastq.gz ${SAMPLE_NAME}_${RUN_ID}_2.filt.fastq.gz | samtools view -b -f 3 | samtools sort > ${SAMPLE_NAME}_${RUN_ID}.sorted.bam

    # Add Read group
    java -jar ${PICARD} AddOrReplaceReadGroups I=${SAMPLE_NAME}_${RUN_ID}.sorted.bam O=${SAMPLE_NAME}_${RUN_ID}.sorted.RG.bam \
                                         RGID=${RUN_ID} RGLB=${LIBRARY_NAME} RGPL=${INSTRUMENT_PLATFORM} \
                                         RGPU=${RUN_NAME} RGSM=${SAMPLE_NAME} RGPI=${INSERT_SIZE}

    # Append the file name (.bam) to the list of alignments that will be merged
    echo ${SAMPLE_NAME}_${RUN_ID}.sorted.RG.bam >> ${SAMPLE_NAME}.bamlist
done

# Merge the list of alignments into a single file
samtools merge -f -b ${INDIVIDUAL}.bamlist ${INDIVIDUAL}.bam

# Index the alignment
samtools index ${INDIVIDUAL}.bam > ${INDIVIDUAL}.bam.bai

# Mark Duplicate reads
java -jar ${PICARD} MarkDuplicates \
	REMOVE_DUPLICATES=FALSE \
	I=${INDIVIDUAL}.bam \
	O=${INDIVIDUAL}.marked_dups.bam \
	M=${INDIVIDUAL}.dup_metrics.txt

# Obtain index
java -jar ${PICARD} BuildBamIndex \
	INPUT=${INDIVIDUAL}.marked_dups.bam


######################################
## Local realignement around indels ##
######################################

# Find regions that need to be realigned
java -jar ${GATK} -T RealignerTargetCreator \
	-R ${REF_GENOME} \
	-I ${INDIVIDUAL}.marked_dups.bam \
	--known ${KNOWN_INDELS} \
	-o ${INDIVIDUAL}.target_intervals.list 

# Perform local realignement
java -jar ${GATK} -T IndelRealigner \
  -R ${REF_GENOME} \
	-known ${KNOWN_INDELS} \
	-I ${INDIVIDUAL}.marked_dups.bam \
	-targetIntervals ${INDIVIDUAL}.target_intervals.list \
	-o ${INDIVIDUAL}.realigned_reads.bam


################################
## Base quality recalibration ##
################################

# Generate a base recalibration table to analyse patterns of covariation in the dataset
java -jar ${GATK} -T BaseRecalibrator \
	-R ${REF_GENOME} \
	-knownSites ${KNOWN_SNP} \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov ContextCovariate \
	-I ${INDIVIDUAL}.realigned_reads.bam \
	-o ${INDIVIDUAL}.recal_data.table

# Perform base quality recalibration
java -jar ${GATK} -T PrintReads \
	-R ${REF_GENOME} \
	-I ${INDIVIDUAL}.realigned_reads.bam \
	-BQSR ${INDIVIDUAL}.recal_data.table \
	-o ${INDIVIDUAL}.recal_reads.bam


###################
## Call variants ##
###################

# Perform variant calling
java -jar ${GATK} -T HaplotypeCaller \
                  -R ${REF_GENOME} \
                  -I ${INDIVIDUAL}.recal_reads.bam \
                  -o ${INDIVIDUAL}.g.vcf  \
                  --genotyping_mode DISCOVERY \
                  -variant_index_type LINEAR \
                  -variant_index_parameter 128000 \
                  --emitRefConfidence GVCF
done

# Perform joint variant calling
java -jar ${GATK} \
   -T GenotypeGVCFs \
   -R ${REF_GENOME} \
   --variant HG02024.g.vcf \
   --variant HG02025.g.vcf \
   --variant HG02026.g.vcf \
 	 -o trio.vcf

#It is here impossible to recalibrate variants becauuse we deal with exomes

########################################
## Post-processing: Analysis of trios ##
########################################

# Phase trios
java -jar ${GATK} -T PhaseByTransmission \
	-R ${REF_GENOME} \
	--variant trio.vcf \
	-ped ${PEDIGREE}.txt \
	-o trio.phased.vcf

# Evaluate variants by computing control metrics
java -jar ${GATK} -T VariantEval \
	-R ${REF_GENOME} \
	--eval:set1 trio.vcf \
	--eval:set2 trio.phased.vcf \
	-o trio.phased.VE.txt

# Tabulate the number of sites which overlap and share alleles
 java -jar GenomeAnalysisTK.jar \
   -T GenotypeConcordance \
   -R ${REF_GENOME} \
   -eval trio.vcf \
   -comp trio.phased.vcf \
   -o triotabulate.txt