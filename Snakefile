configfile: "config.yaml"
with open('../analyses/zliu_run_accessions.txt') as f:
    samples= [sample for sample in f.read().split('\n') if len(sample) > 0]  # we dont want empty lines

#Alternatively;
#from numpy import loadtxt
#samples= loadtxt('../analyses/zliu_run_accessions.txt',dtype='str')
print(samples)

#Give the output of the final rule as the input for rule all, or output from each rule as checkpoint
rule all:
    input:
        expand(["fastqc/{sample}_{dir}_fastqc.html","merged_pairs/{sample}.flash.extendedFrags.fastq.gz","kraken2/kraken2.{sample}","kraken2/bracken.{sample}","kraken2/bracken.{sample}.report","kraken2/ziwei_sample_sub-set_multiqc_report_repeat.html","kraken2/ziwei_sample_sub-set_multiqc_report_repeat_data","kraken2/heatmap.pdf"],sample= samples, dir=[1,2])
# Perform fastqc for raw fastq

rule fastqc:
    input:
        fq1="../data/sra_fastq/{sample}_1.fastq.gz",
        fq2="../data/sra_fastq/{sample}_2.fastq.gz"
    output:
        "fastqc/{sample}_1_fastqc.html",
        "fastqc/{sample}_2_fastqc.html"
    params:
        outdir="fastqc"    
    shell:
        "module load fastqc;fastqc --outdir {params.outdir} --threads 4 --noextract {input.fq1} {input.fq2}"

#merge by flash2
rule flash2:
    input:
        fq1="../data/sra_fastq/{sample}_1.fastq.gz",
        fq2="../data/sra_fastq/{sample}_2.fastq.gz"
    output:
        "merged_pairs/{sample}.flash.extendedFrags.fastq.gz"
    params:
        outdir="merged_pairs"
    shell:
        "module load flash2;flash2 -z --output-directory {params.outdir} --threads 4 --output-prefix={wildcards.sample}.flash {input.fq1} {input.fq2} 2>&1 |tee ziwei_flash2.log"

#Kraken2 k-mer read taxonomic assignation
rule kraken2:
    input:
        fq1="../data/sra_fastq/{sample}_1.fastq.gz",
        fq2="../data/sra_fastq/{sample}_1.fastq.gz"
    output:
        out="kraken2/kraken2.{sample}",
        report="kraken2/kraken2.{sample}.report"
    threads:2
    resources: cpus=2, mem_mb=80000, time_min=20
    shell:
         "module load kraken2;kraken2 --threads {threads} --db /shared/projects/form_2022_19/kraken2/arch_bact_vir_hum_protoz_fung/ --paired --output {output.out} --report {output.report} --gzip-compressed {input.fq1} {input.fq2}"

#bracken to estimate species relative abundance in each sample.

rule bracken:
    input:
        "kraken2/kraken2.{sample}.report"
    output:
        out="kraken2/bracken.{sample}",
        report="kraken2/bracken.{sample}.report"
    resources: cpus=4,mem_mb=60000, time_min=20
    shell:
        "module load bracken;srun --job-name={wildcards.sample} bracken -d /shared/projects/form_2022_19/kraken2/arch_bact_vir_hum_protoz_fung/ -i {input}  -o {output.report} -w {output.out} -r 50 -l S -t 5 "

# simple integration: add kraken2 results to multiQC
rule multiqc:
    input:
        expand(["kraken2/kraken2.{sample}","kraken2/kraken2.{sample}.report","kraken2/bracken.{sample}","kraken2/bracken.{sample}.report"],sample=samples) 
    output:
        report="kraken2/ziwei_sample_sub-set_multiqc_report_repeat.html",
        data="kraken2/ziwei_sample_sub-set_multiqc_report_repeat_data"
    shell:
        "module load multiqc; multiqc --force --title 'ziwei_sample_sub-set' {input} -o kraken2/"

# Make heatmap using R script
rule makeheatmap:
    input:
        expand("kraken2/bracken.{sample}", sample=samples)
    output:
        "kraken2/heatmap.pdf"
    script:
        "makeheatmap.R"
        
