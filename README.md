# EPITOME: Enhanced Phylogenetic Inference Through Optimized Mapping Efficiency

EPITOME condenses a diverse set of DNA sequences (species-level or closer) into discrete, composite sequences that represent the overarching diversity of the dataset. In other words, EPITOME creates sequences that are the _epitome_ of the dataset diversity. This is accomplished by clustering the input based on pairwise genetic distances and then averaging the nucleotides at each genomic position (ties selected at random). When the genetic distance is based on read mapping efficiency, EPITOME can be used to create a set of reference genomes for consensus-based assembly pipelines, like [VAPER](https://github.com/DOH-JDJ0303/VAPER) or [viralrecon](https://github.com/nf-core/viralrecon).

See the [wiki](https://github.com/DOH-JDJ0303/epitome/wiki) for more information.

# Quick Start
## Step 1. Create your samplesheet
> Note: Nextflow requires absolute paths in samplesheets
Create a samplesheet containing the taxa name, genome segment, path to a multi-fasta file of sequences for the taxa, and the expected sequence length (default within 15%).
`samplesheet.csv`:
```
taxa,segment,assembly,length
Influenza_A,HA,flu-a_HA_NCBI_2024-4-1.fasta,1710
Influenza_A,NA,flu-a_NA_NCBI_2024-4-1.fasta,1400
Measles,wg,measles_NCBI_2024-4-1.fasta,16000
```
## Step 2. Run EPITOME
Run EPITOME using the command below.
> Note: See the [wiki](https://github.com/DOH-JDJ0303/epitome/wiki) for how to assign references with existing subtype classifications (e.g., H1-H9) using the `--seeds` parameter.
```
nextflow run DOH-JDJ0303/epitome \
    -r main \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results
```
