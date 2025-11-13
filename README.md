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

Quality control of Illumina reads using FastQC/MultiQC, following the structure:

```bash
module load fastqc/0.12.1

SAM=CEDRELA-3-LEAF
fastqc --threads 4 --outdir ./${DIR}_fastqc/ ../01_Raw_Reads/${SAM}_R1_001.fastq.gz ../01_Raw_Reads/${SAM}_R2_001.fastq.gz
```

## Trimming

Trimmomatic was used to remove low quality and residual adapter sequence. 

```bash


```

## Assembling Transcriptomes

To create a single reference transcriptome for all samples, each sample was separately assembled, pool all the resulting transcripts, cluster them into groups that hopefully represent single genes, and finally select a single "best" representative transcript for each gene.


```bash
SAM=K21

Trinity --seqType fq \
  --left ../02_Quality_Control/trim_${SAM}_R1.fastq.gz \
  --right ../02_Quality_Control/trim_${SAM}_R2.fastq.gz \
  --min_contig_length 300 \
  --CPU 36 \
  --max_memory 100G \
  --output trinity_${SAM} \
  --full_cleanup
```

## Identifying the Coding Regions



## Determining and Removing Redundant Transcripts


## Evaluating the Assembly

## Quantifying gene expression

### Creating an index

### Counting reads mapping to transcripts

## Functional Annotation

## Differential expression analysis

