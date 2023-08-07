################################################################################
# BGC discovery using antiSMASH v7
#
# Alex Huebner, 03/08/23
################################################################################

from glob import glob
import os

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
SAMPLES = {os.path.basename(fn).split("-")[0]: fn
           for fn in glob("/mnt/archgen/users/richtermeier/bachelor_thesis/genome_assembly_test/ancient_metagenome_assembly/results_assembly/alignment/megahit/*-megahit.fasta.gz")}
################################################################################

rule all:
    input:
        expand("04-analysis/antismash/{sample}/{sample}.tsv", sample=SAMPLES)

#### Filter FastA file #########################################################

rule filter_fasta:
    output:
        temp("tmp/bakta/{sample}.fasta")
    message: "Discard all contigs shorter 2 kb: {wildcards.sample}"
    container: "https://depot.galaxyproject.org/singularity/bioawk:1.0--h7132678_7"
    params:
        fasta = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        bioawk -c fastx '{{if (length($seq) >= 2000){{print ">" $name "\\n" $seq}}}}' {params.fasta} > {output}
        """

################################################################################

#### Annotate contigs using Bakta ##############################################

rule download_bakta_db:
    # ln -s /mnt/archgen/microbiome_paleobiotech/calcBGCecoevo/03-data/refdbs/bakta 03-data/refdbs/ && touch 03-data/refdbs/bakta/downloaded
    output:
        touch("03-data/refdbs/bakta_downloaded")
    message: "Download and prepare the reference DB for Bakta"
    conda: "ENVS_bakta.yaml"
    resources:
        mem = 4,
        cores = 1
    params:
        dir = "03-data/refdbs/bakta"
    shell:
        """
        bakta_db download --output {params.dir}
        """

rule bakta:
    # executed with --use-singularity --singularity-args "-B /mnt/archgen/microbiome_paleobiotech/calcBGCecoevo/03-data/refdbs/bakta/db"
    input:
        fasta = "tmp/bakta/{sample}.fasta",
        db = "03-data/refdbs/bakta_downloaded"
    output:
        "04-analysis/bakta/{sample}.gff3"
    message: "Annotate contigs using BAKTA: {wildcards.sample}"
    container: "/mnt/archgen/tools/singularity/containers/depot.galaxyproject.org-singularity-bakta-1.7.0--pyhdfd78af_0.img"
    resources:
        mem = 128,
        cores = 8
    params:
        prefix = lambda wildcards: wildcards.sample,
        outdir = "04-analysis/bakta",
        dbdir = "03-data/refdbs/bakta/db",
        extra = "--keep-contig-headers --meta --skip-trna --skip-tmrna --skip-rrna --skip-ncrna --skip-ncrna-region --skip-plot --skip-ori --skip-gap --skip-crispr"
    threads: 8
    shell:
        """
        bakta -p {params.prefix} \
            --db {params.dbdir} \
            --output {params.outdir} \
            {params.extra} \
            --threads {threads} \
            {input.fasta}
        """

################################################################################

#### AntiSMASH #################################################################

rule antismash:
    input:
        "04-analysis/bakta/{sample}.gff3"
    output:
        zip = "04-analysis/antismash/{sample}/{sample}.zip",
        gbk = "04-analysis/antismash/{sample}/{sample}.gbk"
    message: "Run antiSMASH: {wildcards.sample}"
    #container: "docker://antismash/standalone"
    singularity: "/mnt/archgen/users/huebner/containers/antismash_standalone_v7.0.0.sif"
    resources:
        mem = 50,
        cores = 8
    params:
        gbff = "04-analysis/bakta/{sample}.gbff",
        outdir = "04-analysis/antismash/{sample}"
    threads: 8
    log: "04-analysis/antismash/{sample}/antiSMASH.log"
    shell:
        """
        antismash \
            --clusterhmmer --tigrfam --asf --cc-mibig --cb-general --rre \
            -c {threads} \
            --genefinding-tool none \
            --cb-knownclusters \
            --output-dir {params.outdir} \
            --output-basename {wildcards.sample} \
            --minlength 3000 --allow-long-headers \
            --logfile {log} \
            {params.gbff}
        """

rule comBGC:
    input:
        "04-analysis/antismash/{sample}/{sample}.gbk"
    output:
        "04-analysis/antismash/{sample}/{sample}.tsv"
    message: "Parse antiSMASH output to TSV: {wildcards.sample}"
    conda: "ENVS_comBGC.yaml"
    resources:
        mem = 4,
        cores = 1
    params:
        outdir = "04-analysis/antismash/{sample}"
    threads: 1
    shell:
        """
        $HOME/.nextflow/assets/nf-core/funcscan/bin/comBGC.py -i {input} -o {params.outdir} && \
        mv {params.outdir}/combgc_summary.tsv {output}
        """

################################################################################
