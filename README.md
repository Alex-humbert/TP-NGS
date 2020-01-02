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
The first steps are 
