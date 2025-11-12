# tropical-transcriptomes
Bioinformatics pipeline and configs for de novo transcriptome assembly of Alzatea verticillata, Cedrela montana, Graffenrieda emarginata, and Handroanthus chrysanthus.

The same code runs for all species using a ${spp} placeholder.

## Data overview

## Data availability

Raw RNA-Seq reads are available at the NCBI SRA. Example: PRJNA1359388 (H. chrysanthus). Añade aquí el resto de accesiones cuando estén disponibles.

## Table of Contents





Explicación de dónde se corrieron estos scripts (xanadú) y usando SLURM.

Folder structure:

```text
${spp}/
├── 01_Raw_Reads
├── 02_Quality_Control
├── 03_Assembly
├── 04_Coding_Regions
├── 05_Clustering
├── 06_Assembly_Control
├── 07_Annotation
├── 08_Expression
├── 09_DEA
├── 10_Functional_analysis
```

## Quality Control

Quality control of Illumina reads using FastQC/MultiQC (see [`scripts/01_fastqc.sh`](scripts/01_fastqc.sh) and [`scripts/01_multiqc.sh`](scripts/01_multiqc.sh))

## Trimming

## Assembling Transcriptomes


## Identifying the Coding Regions

## Determining and Removing Redundant Transcripts


## Evaluating the Assembly

## Quantifying gene expression

### Creating an index

### Counting reads mapping to transcripts

## Functional Annotation

## Differential expression analysis

