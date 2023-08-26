################################################################################
# Taxonomic profiling of contigs with antiSMASH hits and ancient DNA damage
# analysis
#
# Alex Huebner, 16/08/23
################################################################################

from glob import glob
import os
import re

import pandas as pd
import pyfastx

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
ANTISMASH_DIR = "/mnt/archgen/users/richtermeier/bachelor_thesis/04-analysis/antismash"
ASSEMBLY_DIR = "/mnt/archgen/users/richtermeier/bachelor_thesis/04-analysis/snakemake/results_assembly"
SAMPLES = [os.path.basename(fn) for fn in glob(f"{ANTISMASH_DIR}/*")]
################################################################################

rule all:
    input:
        "05-results/ANNO_BGC_contig_taxclassification.tsv",
        "05-results/QUAL_pyDamage_results_BGCcontigs.tsv"

#### Identify the list of contigs with antiSMASH hits ##########################

rule extract_antismash_contigs:
    output:
        contigs = "tmp/mmseqs_pydamage/{sample}.antismash_contigs.txt",
        fasta = "tmp/mmseqs_pydamage/{sample}.antismash_contigs.fasta"
    message: "Extract the contigs that were identified by antiSMASH: {wildcards.sample}"
    resources:
        mem = 4,
        cores = 1
    params:
        antismash = lambda wildcards: f"{ANTISMASH_DIR}/{wildcards.sample}/{wildcards.sample}.tsv",
        fasta = lambda wildcards: f"{ASSEMBLY_DIR}/alignment/megahit/{wildcards.sample}-megahit.fasta.gz"
    threads: 1
    run:
        antismash = pd.read_csv(params.antismash, sep="\t")
        antismash['contig_num_id'] = antismash['Contig_ID'] \
            .str.extract(r'NODE_([0-9]+)_length_[0-9]+_cov_[0-9\.]+').astype(int)
        antismash = antismash.sort_values(['contig_num_id'])
        contigs = antismash['Contig_ID'].tolist()

        # Write list of contigs for subsetting BAM file
        extract_length = re.compile(r'NODE_[0-9]+_length_([0-9]+)_cov_[0-9\.]+')
        with open(output.contigs, "wt") as outfile:
            for c in contigs:
                length = extract_length.search(c).group(1)
                outfile.write(f"{c}:1-{length}\n")

        # Subset FastA file
        contigs = set(contigs)
        with open(output.fasta, "wt") as fastafile:
            for name, seq in pyfastx.Fastx(params.fasta):
                if name in contigs:
                    fastafile.write(f">{name}\n{seq}\n")

################################################################################

#### MMSeqs2 taxonomy ##########################################################

rule concat_bins:
    input:
        fastas = expand("tmp/mmseqs_pydamage/{sample}.antismash_contigs.fasta", sample=SAMPLES)
    output:
        temp("tmp/mmseqs_pydamage/mmseqs2/all_contigs.fasta")
    message: "Concatenate all contigs that were binned into a single FastA"
    resources:
        mem = 4,
        cores = 1
    run:
        with open(output[0], "wt") as outfile:
            for sample in SAMPLES:
                # Write to file if not empty
                for name, seq in pyfastx.Fasta(f"tmp/mmseqs_pydamage/{sample}.antismash_contigs.fasta", build_index=False):
                    outfile.write(f">{sample}:{name}\n{seq}\n")

rule createdb_bins:
    input:
        "tmp/mmseqs_pydamage/mmseqs2/all_contigs.fasta"
    output:
        "tmp/mmseqs_pydamage/mmseqs2/all_contigs.contigs"
    message: "Create database of contigs"
    container: "https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0"
    resources:
        mem = 16,
        cores = 1
    shell:
        "mmseqs createdb {input} {output}"

rule screen:
    input:
        "tmp/mmseqs_pydamage/mmseqs2/all_contigs.contigs"
    output:
        "tmp/mmseqs_pydamage/mmseqs2/all_contigs.mmseqs2_gtdb.index"
    message: "Assign taxonomy via the GTDB for contigs"
    container: "https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0"
    resources:
        mem = 750,
        cores = 24
    params:
        tmpdir = "tmp/mmseqs_pydamage/mmseqs2_tmpdir",
        prefix = "tmp/mmseqs_pydamage/mmseqs2/all_contigs.mmseqs2_gtdb",
        gtdb_db = "/mnt/archgen/users/huebner/refdbs/mmseqs2_gtdb_r207_db/mmseqs2_gtdb_r207_db"
    threads: 24
    shell:
        """
        mmseqs taxonomy \
            {input} \
            {params.gtdb_db} \
            {params.prefix} \
            {params.tmpdir} \
            -a \
            --tax-lineage 1 \
            --lca-ranks kingdom,phylum,class,order,family,genus,species \
            --majority 0.5 \
            --vote-mode 1 \
            --orf-filter 1 \
            --remove-tmp-files 1 \
            --threads {threads}
        """

rule create_tsv:
    input:
        contigs = "tmp/mmseqs_pydamage/mmseqs2/all_contigs.contigs",
        assignments = "tmp/mmseqs_pydamage/mmseqs2/all_contigs.mmseqs2_gtdb.index"
    output:
        temp("tmp/mmseqs_pydamage/mmseqs2/all_contigs.mmseqs2_gtdb.tsv")
    message: "Convert MMSeqs2 GTDB results to TSV"
    container: "https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0"
    resources:
        mem = 8,
        cores = 1
    params:
        prefix = "tmp/mmseqs_pydamage/mmseqs2/all_contigs.mmseqs2_gtdb"
    shell:
        """
        mmseqs createtsv {input.contigs} {params.prefix} {output}
        """

rule mmseqs2_annotatetsv:
    input:
        "tmp/mmseqs_pydamage/mmseqs2/all_contigs.mmseqs2_gtdb.tsv"
    output:
        "05-results/ANNO_BGC_contig_taxclassification.tsv"
    message: "Add header to MMSeqs2 table"
    run:
        res = pd.read_csv(input[0], sep="\t", header=None,
                          usecols=list(range(8)) + [9],
                          names=['contig', 'NCBItaxID', 'NCBIrank', 'NCBItaxName',
                                 'nFrags', 'retainedFrags', 'taxassignedFrags',
                                 'fractionAgreement', 'lineage'])
        res['sample'] = res['contig'].str.split(":").str[0]
        res['contig'] = res['contig'].str.split(":").str[1]
        res.iloc[:, [9] + list(range(9))] \
            .to_csv(output[0], sep="\t", index=False)

################################################################################

#### PyDamage analysis #########################################################

rule subset_bam_to_mag:
    input:
        "tmp/mmseqs_pydamage/{sample}.antismash_contigs.txt"
    output:
        temp("tmp/mmseqs_pydamage/{sample}.antismash.bam")
    message: "Subset BAM file for contigs of bin {wildcards.sample}"
    group: "prep_bamfile"
    container: "https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0"
    resources:
        mem = 2,
        cores = 1
    params:
        bam = lambda wildcards: f"{ASSEMBLY_DIR}/alignment/megahit/{wildcards.sample}.sorted.dedup.bam"
    shell:
        """
        cat {input} | xargs samtools view -o {output} -hu {params.bam}
        """

rule remove_extra_headerlines:
    input:
        contiglist = "tmp/mmseqs_pydamage/{sample}.antismash_contigs.txt",
        bam = "tmp/mmseqs_pydamage/{sample}.antismash.bam"
    output:
        bam = temp("tmp/mmseqs_pydamage/{sample}.antismash_clean.bam"),
        bai = temp("tmp/mmseqs_pydamage/{sample}.antismash_clean.bam.bai")
    message: "Remove contigs not belonging to bin to allow for faster iteration: {wildcards.sample}"
    group: "prep_bamfile_depth"
    resources:
        mem = 2,
        cores = 0
    params:
        assembler = "metaspades"
    wrapper:
        "file:/mnt/archgen/users/huebner/automatic_MAG_refinement/workflow/wrappers/remove_extra_headerlines"


rule pydamage:
    input:
        bam = "tmp/mmseqs_pydamage/{sample}.antismash_clean.bam",
        bai = "tmp/mmseqs_pydamage/{sample}.antismash_clean.bam.bai"
    output:
        "tmp/mmseqs_pydamage/{sample}.pydamage.tsv"
    message: "Quantify the ancient DNA damage with pyDamage: {wildcards.sample}"
    container: "https://depot.galaxyproject.org/singularity/pydamage:0.70--pyhdfd78af_0"
    resources:
        mem = 8,
        cores = 4
    params:
        tmpdir = "tmp/mmseqs_pydamage/pydamage_{sample}"
    threads: 4
    shell:
        """
        mkdir -p {params.tmpdir} && \
        pydamage --outdir {params.tmpdir} analyze -p {threads} -f {input.bam} && \
        mv {params.tmpdir}/pydamage_results.csv {output} && \
        rmdir {params.tmpdir}
        """

rule summarise_pydamage:
    input:
        expand("tmp/mmseqs_pydamage/{sample}.pydamage.tsv", sample=SAMPLES)
    output:
        "05-results/QUAL_pyDamage_results_BGCcontigs.tsv"
    message: "Summarise the pyDamage results into a single table"
    resources:
        mem = 4,
        cores = 1
    run:
        pd.concat([pd.read_csv(f"tmp/mmseqs_pydamage/{s}.pydamage.tsv", sep=",") \
                   .assign(sample=s)
                   for s in SAMPLES]) \
        .sort_values(['sample']) \
        .iloc[:, [36] + list(range(36))] \
        .to_csv(output[0], sep="\t", index=False, float_format="%.3f")

################################################################################
