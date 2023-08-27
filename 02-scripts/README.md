# This directory contains the scripts used in the bachelore_thesis project
## Data preparation
#### downloading data
This is a simple bash script for downloading the data from ENA using cURL.  
##### Scripts:  
- `curl_download_script.sh`
##### Commands:  
```
bash curl_download_script.sh
```
#### nextflow/eager
nexftflow/eager is used for multiple preprocessing steps. Those involve removing adapter sequences, quality checking, highlighting hostDNA etc.
##### Commands:  
```
nextflow run nf-core/eager -r 2.4.6 \
	-profile eva,archgen \
	--input '01-documentation/nfcore_eager_samplesheet.tsv' \
	--fasta '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.fa' \
	--fasta_index '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.fa.fai' \
	--bwa_index '/mnt/archgen/Reference_Genomes/Human/hs37d5' \
	--seq_dict '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.dict' \
	--skip_preseq \
	--skip_deduplication \
	--skip_damage_calculation \
	--skip_qualimap \
	--skip_collapse \
	--complexity_filter_poly_g \
	--bwaalno 1 --bwaalnl 32 \
	--outdir '04-analysis-eager' 
```
#### Removal of host DNA with snakemake
This script is used for removing the host DNA. Since out samples should only include DNA from microbes, all DNA from the DNA has to be removed. This is done using a Snakemake pipeline and its ENVS.
##### Scripts:
- `PREP_remove_hostDNA.Snakefile`
- `ENVS_bioawk.yaml`
- `ENVS_samtools.yaml`
##### Commands:  
```
snakemake -s 02-scripts/PREP_remove_hostDNA.Snakefile --use-conda --conda-prefix conda --profile sge_archgenq --cores 10 -j 8
```
## Composition analysis
#### MetaPhlAn4
MetaPhlAn4 is used for taxonomic classification. The Snakefile contains the configurations to run the program itself. In order to that it needs a environment with for the right package versions which first need to be created. This step requires the .yaml file, which contains the environment.
##### Scripts:  
- `COMP_MetaPhlAn4.Snakefile`
- `ENVS_MetaPhlAn4.yaml`  
##### Commands:  
Create environment:  
```
snakemake -s 02-scripts/COMP_MetaPhlAn4.Snakefile --use-conda --conda-prefix conda --profile sge_archgenq --cores 10 -j 8 --latency-wait 60 --conda-create-envs-only
```
Start program:  
```
snakemake -s 02-scripts/COMP_MetaPhlAn4.Snakefile --use-conda --conda-prefix conda --profile sge_archgenq --cores 10 -j 8 --latency-wait 60
```

## de novo assembly, binning and contig annotation
#### Snakemake pipeline [MEGAHIT]
In this step the DNA of the microbacteria is assembled back together from the reads. By overlapping short reads, long contigs can be created. We used MEGAHIT as an assembler, since it allows time and memory efficient computation.
##### Scripts:  
- `samplesheet_mag.csv`
- `ASMB_nfcore_mag.sh`
- `ENVS_taxpasta.yaml`
##### Commands:
```
snakemake --configfile /mnt/archgen/users/richtermeier/bachelor_thesis/genome_assembly_test/ancient_metagenome_assembly/config/config.yaml \
	--use-conda --conda-prefix /mnt/archgen/users/huebner/ancient_metagenome_assembly/conda \
	--profile sge_archgenq -j 8 --cores 24 --latency-wait 60 \
	--restart-times 3
```
## functional annotation
#### Snakemake pipeline (antiSMASH)
This step involves annotating genes in the coding regions to find BGCs. This is done by antiSMASH. The commands to run antiSMASH itself is defined in the Snakefile. The environments for the correct packages can be found in the ENVS files.
##### Scripts:  
- `ANNO_antiSMASH.Snakefile`
- `ENVS_comBGC.yaml`
- `ENVS_bakta.yaml`
##### Commands: 
```
snakemake -s 02-scripts/ANNO_antiSMASH.Snakefile \
    --use-conda --conda-prefix /mnt/archgen/users/huebner/conda \
    --use-singularity --singularity-prefix /mnt/archgen/users/huebner/containers \
    --singularity-args '-B /mnt/archgen/microbiome_paleobiotech/calcBGCecoevo/03-data/refdbs/bakta/db' \
    --latency-wait 60 --profile sge_archgenq -j 10 \
    --rerun-triggers mtime
```

## contig annotation 
This script identifies all contigs for which antiSMASH detected a BGC and taxonomically classifies them against the GTDB classification. 
##### Scripts:
- `QUAL_MMseqs2_pyDamage.Snakefile`
##### Commands:
```
snakemake -s 02-scripts/QUAL_MMseqs2_pyDamage.Snakefile \
    --use-conda --conda-prefix /mnt/archgen/users/huebner/conda \
    --use-singularity --singularity-prefix /mnt/archgen/users/huebner/containers \
    --singularity-args '-B /mnt/archgen/users/huebner/refdbs/mmseqs2_gtdb_r207_db' \
    --profile sge_archgenq --latency-wait 30 -j 10
```

## Data visualization
The `plots` folder includes all scripts for visualizing the data in plots using python

- `preprocessing_results.py`: plotting the distribution of the read length among all samples and percentage of reads filtered out during preprocessing
- `metaphlan_profile.py`: bar plot for the most abundant clades as identified by MetaPhlan4
- `assembly_results.py`: different bar plots showing the assembly sucess by looking at the number of contigs and their length
- `MAG_quality.py`: plots the number of high and low quality MAGs and their assigned genus
