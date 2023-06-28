# De novo assembly with contig-correction using MEGAHIT

nextflow run nf-core/mag \
    -r 2.3.2 \
    -profile eva,archgen \
    --input samplesheet.csv \
    --outdir 04-analysis/assembly/megahit \
    --skip_clipping \
    --skip_spades \
    --skip_prodigal \
    --skip_binning \
    --skip_prokka \
    --skip_binqc \
    --ancient_dna \
    -c nfcore_mag.conf

# De novo assembly without contig-correction using metaSPAdes

nextflow run nf-core/mag \
    -r 2.3.2 \
    -profile eva,archgen \
    --input samplesheet.csv \
    --outdir 04-analysis/assembly/metaspades \
    --skip_clipping \
    --spades_options "-k 21,33,55,77" \
    --skip_megahit \
    --skip_prodigal \
    --skip_binning \
    --skip_prokka \
    --skip_binqc \
    --ancient_dna --skip_ancient_damagecorrection \
    -c nfcore_mag.conf
