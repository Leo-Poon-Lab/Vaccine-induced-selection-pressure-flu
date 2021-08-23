# Vaccine induced selection pressure on seasonal influenza in mice 
*Authors*

---

This repository host the scripts for analyzing the nucleotide diversity (π) via deep-sequencing data in the study "Vaccine induced selection pressure on seasonal influenza in mice". The viral population diversity measurements including πN and πS were estimated using [SNPGenie](https://academic.oup.com/bioinformatics/article/31/22/3709/241742). 

Details:

1. `./scripts/0.resequencing.sh`: The raw Illumina fastq reads were aligned to the reference using BWA-MEM2 (v2.0). Primer trimming and quality trimming were performed using [iVar](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1618-7). Variants using for population diversity estimation were called by `samtools mpileup` and `ivar variants`.
2. `./scripts/nu_diversity.R`: Building GTF, VCF and other files needed for input for SNPGenie.

---