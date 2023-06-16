################################################################################
# Project: Training project
# Part: Composition analysis
# Step: Taxonomic profiling with MetaPhlAn
#
# Dependent on:
#   - PREP_remove_hostDNA.Snakefile
#
# Alex Huebner, 18/04/23
################################################################################

import os

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
SAMPLES, = glob_wildcards("03-data/processed_data/{sample}_1.fastq.gz")
################################################################################

rule all:
    input:
        "05-results/COMP_MetaPhlAn4_species_profiles.tsv",
        "05-results/COMP_MetaPhlAn4_alignedreads.tsv"
    
rule metaphlan4_install_db:
    output:
        touch("04-analysis/metaphlan4/installed_database.done")
    message: "Install the CHOCOPhlAn database"
    conda: "ENVS_MetaPhlAn4.yaml"
    shell:
        "metaphlan --install"
        
rule concat_fq:
    output:
        pipe("tmp/metaphlan/{sample}.fastq.gz")
    message: "Concatenate paired-end seq data: {wildcards.sample}"
    group: "metaphlan"
    resources:
        mem = 2,
        cores = 1
    params:
        pe1 = "03-data/processed_data/{sample}_1.fastq.gz",
        pe2 = "03-data/processed_data/{sample}_2.fastq.gz"
    shell:
        "cat {params.pe1} {params.pe2} > {output}"
    
rule metaphlan:
    input:
        db = "04-analysis/metaphlan4/installed_database.done",
        fq = "tmp/metaphlan/{sample}.fastq.gz"
    output:
        sam = "04-analysis/metaphlan4/{sample}.metaphlan.sam.bz2",
        profile = "04-analysis/metaphlan4/{sample}.metaphlan.profile.txt"
    message: "Run MetaPhlAn4 with default settings for sample {wildcards.sample}"
    conda: "ENVS_MetaPhlAn4.yaml"
    group: "metaphlan"
    resources:
        mem = 24,
        cores = 8
    threads: 8
    shell:
        """
        metaphlan \
            {input.fq} \
            --input_type fastq \
            --force \
            --no_map \
            --index mpa_vOct22_CHOCOPhlAnSGB_202212 \
            --ignore_eukaryotes \
            -t rel_ab_w_read_stats \
            --sample_id {wildcards.sample} \
            --read_min_len 35 \
            -s {output.sam} \
            -o {output.profile} \
            --nproc {threads}
        """

rule merge_metaphlan:
    input:
        expand("04-analysis/metaphlan4/{sample}.metaphlan.profile.txt", sample=SAMPLES)
    output:
        "05-results/COMP_MetaPhlAn4_species_profiles.tsv"
    message: "Combine the MetaPhlAn results on species level"
    run:
        profiles = pd.concat([pd.read_csv(fn, sep="\t", skiprows=5) \
                              .rename({'#clade_name': 'clade_name'}, axis=1) \
                              .query('clade_name.str.contains("s__") and ~clade_name.str.contains("t__")') \
                              .assign(sample=os.path.basename(fn).replace(".metaphlan.profile.txt", ""))
                              for fn in input])

        profiles[['sample', 'clade_name', 'clade_taxid', 'relative_abundance',
                  'estimated_number_of_reads_from_the_clade']] \
        .sort_values(['sample', 'relative_abundance'], ascending=[True, False]) \
        .to_csv(output[0], sep="\t", index=False)


rule aligned_reads:
    input:
        "04-analysis/metaphlan4/{sample}.metaphlan.sam.bz2"
    output:
        temp("04-analysis/metaphlan4/{sample}.nreads.txt")
    message: "Determine the number of aligned reads for sample {wildcards.sample}"
    resources:
        mem = 4
    shell:
        """
        bzgrep -v "^@" {input} | wc -l > {output}
        """

rule summary_alignedreads:
    input:
        expand("04-analysis/metaphlan4/{sample}.nreads.txt", sample=SAMPLES)
    output:
        "05-results/COMP_MetaPhlAn4_alignedreads.tsv"
    message: "Summarise the number of aligned reads against the MetaPhlAn4 database"
    resources:
        mem = 4
    params:
        dir = "04-analysis/metaphlan4"
    run:
        nreads = [(os.path.basename(fn).replace(".nreads.txt", ""),
                  open(fn).readline().rstrip())
                  for fn in input]
        nreads_df = pd.DataFrame(nreads, columns=['sample', 'nreads'])
        nreads_df['nreads'] = nreads_df['nreads'].astype(int)
        nreads_df.sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False)
