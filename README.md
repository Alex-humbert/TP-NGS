# TP-NGS
This repository is meant to enable anyone to reproduce what we did during our practical work. We wanted to develop a tool to detect genetic variants, using data from the 1000 Genomes project. Our program enabled us to work on three individuals but could be used for greater scale studies with a few rearrangements. Have a good readthrough !


# 0. Downloading the files
The commands are as follow :
```
git clone https://github.com/Alex-humbert/TP-NGS.git
cd TP-NGS/
```


# 1. Requirements:
Java (version 8)
FastQC (version 0.11.7)
BWA-MEM (version 0.7.17-r1194-dirty)
SAMtools (version 1.9)
IGV (version 2.4.14)

The installation.sh file is useful only if you don't already have the tools preinstalled on your machine. It also contains Picard and GATK tools, that were useful for our further work.


# 2. Mapping
The first steps are designed in downloading all data necessary. This contains the reference genome from the 1000 Genomes project, and here three indivual genomes. Those need to be sometimes unzipped and indexed. Indexation is realised by bwa index. We also download the index of all reads of all individuals in the project.

We then sort this index in order to only keep reads corresponding to our individuals, the sampling method used (here, blood), the machine type used, and last but not least their position. We only want to look at exons, because the model to study variants is way easier this way.

Once the index is filtered, the aim is to get all reads corresponding to exons in the individuals we study. So we go through ou index, accessing the corresponding reads in the fasta sequence, and for each pair of reads (positive and negative), map the read against the reference, and sort them by position order. This is made by bwa mem, and samtools suite (view and sort).

Read groups are then constructed by reanalysing the sorted bam file using Gatk's AddOrReplaceReadGroups; they are groups of sequences (from 1 to... a lot !) that origin from the same run. This limits the number of separate runs and thus the size of the files. Those read groups are then merged into a single bam file, which can contain non-sequenced zones anyway. This is samtools merge.

The final step of this subpart is the indexation of our sequences, using samtools index.




# 3. Variant Calling
This is the real core part of our job. Nevertheless, it starts once again with the acquisition and preparation of data. It begins this time with the downloading of known Indels' and SNPs' lists.

Every sequence used have to be indexed ; although we already did that, the code was more modular that way and we do that again. A sequence dictionnary (an exhaustive list of smaller sequences) is also required and hereby constructed.

If we have duplicates, then variants we be called several times, so we need to mark them using Picard's MarkDuplicates. Yet we won't merge them all ; it would cause the loss of the ability to detect heterozygotic sites.

Now come two steps of quality enhancement. The first one uses indels to realign locally the sequences around them. It is made using two functions: initially, one needs to identify the loci where a realignment could benefit ; and then the realignment itself should be made. This is proposed by the Gatk suite.

Then quite a resembling action is realised using SNPs. This allows the sequences to be modified according to the most credible version. It mostly is the case when facing single false positive variants, for instance due to sequencing failures. The use of known SNPs is permitting to decipher differences to the reference genome.

Once all this is made, we get to the most important functions, which are HaplotypeCaller and GenotypeGVCFs, once again from Gatk's suite. Their combined action create a vcf file containing the list of detected variants in regard to the reference genome. Using the gvcf format will be easier afterwards for us to pile the three individual files !


# 4. Trio Analysis
Here again, a small bit of new data is required, and also does preparation of our preexistent data. We merge all three gvcf files in a single trio.vcf file containing the informations on all three individuals.

We also download the pedigree enabling us to confirm the familial relationship between our three subjects (a woman and her parents). It needs to be slightly modified to be used by Gatk, which explains why we cut its header and only keep the 6 first columns. And we ought to save time this way, so...

PhaseByTransmission will take into account the familial bonds between the individuals in the vcf, and phase the variants by resemblance.

VariantEval will create control metrics.

Those metrics are then used by GenomeAnalysisTK to control the concordance between parental genotypes and the daughter's genotype.
