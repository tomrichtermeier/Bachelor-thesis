# This directory contains the scripts used in the bachelore_thesis project
## Data preparation
### downloading data
This is a simple bash script for downloading the data from ENA using cURL.  
##### Scripts:  
- `curl_download_script.sh`
##### Commands:  
```
bash curl_download_script.sh
```
### nextflow/eager
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
### Removal of host DNA with snakemake
##### Scripts:
- `PREP_remove_hostDNA.Snakefile`
- `ENVS_bioawk.yaml`
- `ENVS_samtools.yaml`
##### Commands:  
```
snakemake -s 02-scripts/PREP_remove_hostDNA.Snakefile --use-conda --conda-prefix conda --profile sge_archgenq --cores 10 -j 8
```
## Composition analysis
### MetaPhlAn4
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
## de novo assembly
### nf-core mag with MEGAHIT and metaSPAdes
##### Scripts:  
- `samplesheet_mag.csv`
- `ASMB_nfcore_mag.sh`
- `ENVS_taxpasta.yaml`
##### Commands:  
```
bash 02-scripts/ASMB_nfcore_mag.sh
```