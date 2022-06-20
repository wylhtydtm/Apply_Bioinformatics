# Applied_Bioinformatics

This is a repository for course 2022-Appl-Bioinfo provided by MEDBIOINFO.se to first-year PhD students

Data source: A real world dataset

41 Gbp of Illumina raw sequence reads (both metagenomes and metatranscriptomes) from 125 human oral swabs from positive or negative COVID PCR patients was published in Dec 2021. 

The publication: *A metagenomics workflow for SARS-CoV-2 identification, co-pathogen detection, and overall diversity* [[1]]

Installation:
sra-tools, snakemake, fastqc, flash2, kraken2, bracken,multiqc, R

To download fastq files from the NCBI SRA or EBI ENA archive's websites using the run identifiders(accession numbers)
```
module load sra-tools
srun --cpus-per-task=8 --time=00:30:00 xargs -a zliu_run_accessions.txt fastq-dump --readids --gzip \
--outdir ../data/sra_fastq/ --disable-multithreading --split-e  

```

To run the pipeline:

```
module load snakemake 
snakemake --snakefile Snakefile --cluster-config config.yaml --jobs 2
```



[1]: https://www.sciencedirect.com/science/article/pii/S1386653221002924


Course organisers: Pascal and Samuel 
