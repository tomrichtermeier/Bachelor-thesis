# This directory contains the scripts used in the bachelore_thesis project
## Data preparation
### downloading data
- `curl_download_script.sh`
```
bash curl_download_script.sh
```
### nextflow/eager
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
- `PREP_remove_hostDNA.Snakefile`
- `ENVS_bioawk.yaml`
- `ENVS_samtools.yaml`
```
snakemake -s 02-scripts/PREP_remove_hostDNA.Snakefile --use-conda --conda-prefix conda --profile sge_archgenq --cores 10 -j 8
```
## Composition analysis
### MetaPhlAn4
- `COMP_MetaPhlAn4.Snakefile`
- `ENVS_MetaPhlAn4.yaml`