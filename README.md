# EPITOME: Enhanced Phylogenetic Inference Through Optimized Mapping Efficiency

EPITOME condenses a diverse set of DNA sequences into discrete, composite sequences that represent the overarching diversity of the dataset. In other words, EPITOME creates sequences that are the _epitome_ of the dataset diversity. This is accomplished by clustering the input based on pairwise genetic distances and then selecting the most common nucleotide at each genomic position (ties selected at random) and / or by selecting the centroid of each cluster. When the genetic distance is based on read mapping efficiency, EPITOME creates a set of reference genomes for consensus-based assembly pipelines, like [VAPER](https://github.com/DOH-JDJ0303/VAPER).

See the [wiki](https://github.com/DOH-JDJ0303/epitome/wiki) for more information.

# Pipeline Overview
```mermaid
---
config:
  theme: 'base'
  logLevel: 'debug'
  themeVariables:
    gitBranchLabel0: '#fff'
    git0: "#264653"
    git1: "#2A9D8F"
    git2: "#E9C46A"
    git3: "#F4A261"
    git4: "#E76F51"
    git5: "#8D6E63"
    git6: "#6D6875"
    git7: "#B5838D"
  gitGraph:
    mainBranchName: "Entry [Main]"
---

gitGraph:
  commit id: "Samplesheet 🗒️"
  commit id: " "
  branch "NCBI [Subworkflow]" order: 1
  checkout "NCBI [Subworkflow]"
  commit id: "Datasets"
  commit id: "EDirect"
  checkout "Entry [Main]"
  commit id: "Samplesheet Data 🗂️"
  checkout "NCBI [Subworkflow]"
  commit id: "epitome_ncbi.py"
  commit id: "NCBI Data 🗂️"
  checkout "Entry [Main]"
  merge "NCBI [Subworkflow]"
  commit id: "epitome_metadata.py" tag: "Format Inputs"

  branch "Create [Subworkflow]" order: 3
  checkout "Create [Subworkflow]"
  commit id: "epitome_qc.py"
  commit id: "epitome_cluster.py" tag: "Centroid [Method]"

  branch "Consensus [Method]" order: 4
  checkout "Consensus [Method]"
  commit id: "mafft"
  commit id: "epitome_consensus.py"
  commit id: "epitome_condense.py"

  checkout "Create [Subworkflow]"
  merge "Consensus [Method]"

  checkout "Entry [Main]"
  commit id: "{outdir}/**/inputs/"

  checkout "Create [Subworkflow]"
  commit id: "epitome_summary.py"
  commit id: "{outdir}/**/outputs/"
```

# Quick Start
## Step 1. Create your samplesheet
> [!NOTE]
> The example below creates a reference set for each taxon using data automatically downloaded from NCBI. It is also possible to supply your own sequence data with or without the NCBI data - learn more [here](https://doh-jdj0303.github.io/epitome-docs/).

`samplesheet.csv`:
```
taxon,segmented
Hantavirus,true
Norovirus,false
```
## Step 2. Run EPITOME
Run EPITOME using the command below.
```
nextflow run DOH-JDJ0303/epitome \
    -r main \
    -profile docker \
    --input samplesheet.csv \
    --outdir results
```


