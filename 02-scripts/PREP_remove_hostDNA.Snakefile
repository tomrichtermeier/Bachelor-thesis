################################################################################
# Project: Bachelor thesis
# Part: Preparation of the data
# Step: Extract non-host DNA after nf-core/eager processing
# Date: 13.06.2023
# Author: Tom Richtermeier, Alexander Hübner
################################################################################

import os
import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = 1
os.environ["OMP_NUM_THREADS"] = 1

#### SAMPLES ###################################################################
eager_tbl = pd.read_csv("01-documentation/nfcore_eager_samplesheet.tsv", sep="\t")
SAMPLES = eager_tbl.groupby(['Sample_Name'])['Library_ID'].apply(list).to_dict()
################################################################################

rule all:
    input:
        "05-results/PREP_Nextflow_EAGER_noReads.tsv"

rule samtools_sort_by_name:
    output:
        pipe("tmp/eager_extract_unmapped/{lib}.nsorted.bam")
    message: "Sort the BAM file by name: {wildcards.lib}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 12,
        cores = 4
    params:
        bam = "04-analysis/eager/mapping/bwa/{lib}_PE.mapped.bam"
    threads: 4
    shell:
        """
        samtools sort -n -o {output} {params.bam}
        """
        
rule samtools_fixmate:
    input:
        "tmp/eager_extract_unmapped/{lib}.nsorted.bam"
    output:
        pipe("tmp/eager_extract_unmapped/{lib}.fixmate.bam")
    message: "Fix mate flags: {wildcards.lib}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 8,
        cores = 4
    threads: 4
    shell:
        """
        samtools fixmate -mcu -@ {threads} {input} {output}
        """

rule extract_unmapped_reads:
    input:
        "tmp/eager_extract_unmapped/{lib}.fixmate.bam"
    output:
        pe1 = temp("tmp/eager_extract_unmapped/{lib}_1.fastq.gz"),
        pe2 = temp("tmp/eager_extract_unmapped/{lib}_2.fastq.gz"),
        pe0 = temp("tmp/eager_extract_unmapped/{lib}_0.fastq.gz")
    message: "Extract all reads for which are not aligned in a proper pair and convert to fastq: {wildcards.lib}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 8,
        cores = 2
    threads: 2
    shell:
        """
        samtools view -uh -e '(flag.paired && (flag.unmap || flag.munmap)) || (!flag.paired && flag.unmap)' {input} | \
        samtools fastq -1 {output.pe1} \
                       -2 {output.pe2} \
                       -0 {output.pe0} -
        """
        
rule concat_fastqs:
    input:
        fqs = lambda wildcards: [f"tmp/eager_extract_unmapped/{lib}_{i}.fastq.gz" for lib in SAMPLES[wildcards.sample] for i in range(3)]
    output:
        pe1 = "03-data/processed_data/{sample}_1.fastq.gz",
        pe2 = "03-data/processed_data/{sample}_2.fastq.gz",
        pe0 = "03-data/processed_data/{sample}_0.fastq.gz"
    message: "Concatenate the FastQ files: {wildcards.sample}"
    params:
        tmpdir = "tmp/eager_extract_unmapped",
        outdir = "03-data/processed_data",
        pe1 = lambda wildcards: " ".join([f"tmp/eager_extract_unmapped/{lib}_1.fastq.gz" for lib in SAMPLES[wildcards.sample]]),
        pe2 = lambda wildcards: " ".join([f"tmp/eager_extract_unmapped/{lib}_2.fastq.gz" for lib in SAMPLES[wildcards.sample]]),
        pe0 = lambda wildcards: " ".join([f"tmp/eager_extract_unmapped/{lib}_0.fastq.gz" for lib in SAMPLES[wildcards.sample]])
    shell:
        """
        cat {params.pe1} > {output.pe1}
        cat {params.pe2} > {output.pe2}
        cat {params.pe0} > {output.pe0}
        """

rule count_reads:
    input:
        pe1 = "03-data/processed_data/{sample}_1.fastq.gz",
        pe2 = "03-data/processed_data/{sample}_2.fastq.gz",
        pe0 = "03-data/processed_data/{sample}_0.fastq.gz"
    output:
        temp("03-data/processed_data/{sample}.n")
    message: "Count the number of reads: {wildcards.sample}"
    conda: "ENVS_bioawk.yaml"
    resources:
        mem = 2
    shell:
        """
        reads_PE1=$(bioawk -c fastx 'END{{print NR}}' {input.pe1})
        reads_PE2=$(bioawk -c fastx 'END{{print NR}}' {input.pe2})
        reads_PE0=0
        echo -e "{wildcards.sample}\t${{reads_PE1}}\t${{reads_PE2}}\t${{reads_PE0}}" > {output}
        """

rule summarise_count_reads:
    input:
        expand("03-data/processed_data/{sample}.n", sample=SAMPLES)
    output:
        "05-results/PREP_Nextflow_EAGER_noReads.tsv"
    message: "Summarise the number of reads per sample"
    run:
        pd.concat([pd.read_csv(fn, sep="\t", header=None, names=['sample', 'R1', 'R2', 'R0'])
                   for fn in input]) \
            .sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False)
