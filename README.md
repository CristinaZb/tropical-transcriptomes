# Tropical-transcriptomes
Bioinformatics pipeline and configs for *de novo* transcriptome assembly of *Alzatea verticillata*, *Cedrela montana*, *Graffenrieda emarginata*, and *Handroanthus chrysanthus*.

## Data overview

We sampled four tissues —cambium, flower buds, leaf buds, and leaves— from four Neotropical tree species (*Alzatea verticillata*, *Cedrela montana*, *Graffenrieda emarginata*, *Handroanthus chrysanthus*). Libraries were sequenced on Illumina NovaSeq 6000 (2×150 bp). Leaf buds and leaves include three biological replicates from the same trees; for other tissues we pooled equal RNA amounts from 2–3 individuals per species. In total, 39 tissue-specific libraries were generated.
Analyses were executed on the Xanadú HPC (Institute for Systems Genomics Computational Biology Core - UConn) using SLURM. The command blocks in this README reproduce the exact parameters we used, following the structure of the runs. We followed the RNAseq_nonmodel tutorial (https://github.com/CBC-UCONN/RNAseq_nonmodel) as the operational reference; for detailed system setup and folder scaffolding, please refer to that tutorial.

## Data availability

Raw RNA-Seq reads are available at the NCBI:

| Species                     | Tissue       | Library         | SRA accession | BioSample accession | BioProject   |
|----------------------------|--------------|-----------------|---------------|---------------------|--------------|
| Handroanthus chrysanthus   | Cambium      | CAM-POOL        | SRR35996128   | SAMN53163668        | PRJNA1359388 |
| Handroanthus chrysanthus   | Flower buds  | FB-POOL         | SRR35996127   | SAMN53163669        | PRJNA1359388 |
| Handroanthus chrysanthus   | Leaf buds    | H04-B1          | SRR35996126   | SAMN53163670        | PRJNA1359388 |
| Handroanthus chrysanthus   | Leaf buds    | H06-B1          | SRR35996124   | SAMN53163672        | PRJNA1359388 |
| Handroanthus chrysanthus   | Leaf buds    | H10-B1          | SRR35996122   | SAMN53163674        | PRJNA1359388 |
| Handroanthus chrysanthus   | Leaf         | H04-L1          | SRR35996125   | SAMN53163671        | PRJNA1359388 |
| Handroanthus chrysanthus   | Leaf         | H06-L1          | SRR35996123   | SAMN53163673        | PRJNA1359388 |
| Handroanthus chrysanthus   | Leaf         | H10-L1          | SRR35996121   | SAMN53163675        | PRJNA1359388 |
| Alzatea verticillata       | Cambium      | ALZATEA-CAM     | SRR36004175   | SAMN53199287        | PRJNA1362998 |
| Alzatea verticillata       | Flower buds  | ALZATEA-BF      | SRR36004174   | SAMN53199288        | PRJNA1362998 |
| Alzatea verticillata       | Flowers      | ALZATEA-FL      | SRR36004173   | SAMN53199289        | PRJNA1362998 |
| Alzatea verticillata       | Leaf buds    | ALZATEA-2-YH    | SRR36004172   | SAMN53199290        | PRJNA1362998 |
| Alzatea verticillata       | Leaf buds    | ALZATEA-3-YH    | SRR36004171   | SAMN53199291        | PRJNA1362998 |
| Alzatea verticillata       | Leaf buds    | ALZATEA-6-YH    | SRR36004170   | SAMN53199292        | PRJNA1362998 |
| Alzatea verticillata       | Leaf         | ALZATEA-6-LEAF  | SRR36004169   | SAMN53199293        | PRJNA1362998 |
| Alzatea verticillata       | Leaf         | ALZATEA-3-LEAF  | SRR36004168   | SAMN53199294        | PRJNA1362998 |
| Alzatea verticillata       | Leaf         | ALZATEA-2-LEAF  | SRR36004167   | SAMN53199295        | PRJNA1362998 |
| Cedrela montana            | Cambium      | CEDRELA-CAM     | SRR36016291   | SAMN53202879        | PRJNA1363130 |
| Cedrela montana            | Flower buds  | CEDRELA-BF      | SRR36016290   | SAMN53202880        | PRJNA1363130 |
| Cedrela montana            | Flowers      | CEDRELA-FL      | SRR36016289   | SAMN53202881        | PRJNA1363130 |
| Cedrela montana            | Leaf buds    | CEDRELA-3-YH    | SRR36016288   | SAMN53202882        | PRJNA1363130 |
| Cedrela montana            | Leaf buds    | CEDRELA-5-YH    | SRR36016287   | SAMN53202883        | PRJNA1363130 |
| Cedrela montana            | Leaf buds    | CEDRELA-6-YH    | SRR36016286   | SAMN53202884        | PRJNA1363130 |
| Cedrela montana            | Leaf         | CEDRELA-3-L     | SRR36016285   | SAMN53202885        | PRJNA1363130 |
| Cedrela montana            | Leaf         | CEDRELA-6-L     | SRR36016284   | SAMN53202886        | PRJNA1363130 |
| Cedrela montana            | Leaf         | CEDRELA-5-L     | SRR36016283   | SAMN53202887        | PRJNA1363130 |
| Graffenrieda emarginata    | Flowers      | GRAFF-FL        | SRR36016315   | SAMN53197954        | PRJNA1362976 |
| Graffenrieda emarginata    | Leaf buds    | GRAFF-YH        | SRR36016314   | SAMN53197964        | PRJNA1362976 |
| Graffenrieda emarginata    | Leaf         | GRAFF-LEAF      | SRR36016313   | SAMN53198155        | PRJNA1362976 |


## Table of Contents



## Quality Control

Quality control of Illumina reads with ``FastQC`` and aggregation with ``MultiQC``.

From the ``02_Quality_Control/`` directory:

```bash
module load fastqc/0.12.1

# Single sample
SAMPLE="${SPS}_${ID}"

# FastQC → raw_fastqc/
fastqc --outdir ./raw_fastqc/ ../01_Raw_Reads/${SAMPLE}_R1_001.fastq.gz ../01_Raw_Reads/${SAMPLE}_R2_001.fastq.gz

# MultiQC → raw_multiqc/
multiqc --outdir raw_multiqc ./raw_fastqc/
```

### Taxonomic filtering (Kraken2)

Classify reads and split into unclassified (kept for plant assembly) and classified (non-plant) fractions.

**Inputs:**

``${SAMPLE}_R1_001.fastq.gz``
``${SAMPLE}_R2_001.fastq.gz``

**Outputs:**

``${SAMPLE}_unclassified#.fastq``
``${SAMPLE}_classified#.fastq``
``${SAMPLE}_report.txt``

Script structure for Kraken2:

```bash
module load kraken/2.1.2

SAMPLE="${file}"
THREADS=16          
KRAKEN2_DB="${KRAKEN2_DB}"  # b+a+v+f

mkdir -p kraken_output/kraken_unclassified kraken_output/kraken_classified

kraken2 --db "${KRAKEN2_DB}" \
  --paired "${SAMPLE}_R1_001.fastq.gz" "${SAMPLE}_R2_001.fastq.gz" \
  --use-names \
  --threads "${THREADS}" \
  --unclassified-out kraken_output/kraken_unclassified/${SAMPLE}_unclassified#.fastq \
  --classified-out   kraken_output/kraken_classified/${SAMPLE}_classified#.fastq \
  --report           kraken_output/${SAMPLE}_report.txt
```

## Trimming

Adapter and quality trimming on Kraken2-unclassified reads with Trimmomatic 0.39.

**Inputs:**

``${SAMPLE}_unclassified_1.fastq``
``${SAMPLE}_unclassified_2.fastq``

**Outputs:**

Paired: ``trim_${SAMPLE}_R1.fastq.gz``
``trim_${SAMPLE}_R2.fastq.gz``

Singletons: ``singles_trim_${SAMPLE}_R1.fastq.gz`` 
``singles_trim_${SAMPLE}_R2.fastq.gz``

Script structure for Trimmomatic:

```bash
module load Trimmomatic/0.39

THREADS=16
ADAPTERS="TruSeq3-PE-2.fa"

java -jar "$Trimmomatic" PE -threads "${THREADS}" \
  kraken_output/kraken_unclassified/${SAMPLE}_unclassified_1.fastq \
  kraken_output/kraken_unclassified/${SAMPLE}_unclassified_2.fastq \
  trim_${SAMPLE}_R1.fastq.gz singles_trim_${SAMPLE}_R1.fastq.gz \
  trim_${SAMPLE}_R2.fastq.gz singles_trim_${SAMPLE}_R2.fastq.gz \
  ILLUMINACLIP:${ADAPTERS}:2:30:10 \
  HEADCROP:10 SLIDINGWINDOW:4:25 MINLEN:45
```

Parameter notes:

``ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10`` — removes Illumina adapters using the TruSeq PE list; allows 2 seed mismatches, requires palindrome clip score ≥30 (for overlapping PE reads) and simple clip score ≥10 for standard matches.

``HEADCROP:10`` — drops the first 10 bases from each read to eliminate low-quality/biased starts.

``SLIDINGWINDOW:4:25`` — scans with a 4-bp window and trims once the average Phred quality <25 (~0.3% error per base).

``MINLEN:45`` — discards reads shorter than 45 bp after trimming to avoid spurious mappings/assemblies.

## Assembling Transcriptomes

To create a single reference transcriptome for all samples, each sample was separately assembled, resulting transcripts were pooled and clustered into groups that represent single genes. Finally, select a single "best" representative transcript for each gene. *De novo* transcriptome assembly with Trinity 2.15.0 on trimmed paired reads.

**Inputs:**

``trim_${SAMPLE}_R1.fastq.gz`` 
``trim_${SAMPLE}_R2.fastq.gz``
(singletons not used here)

**Outputs:**

``trinity_out_dir/Trinity.fasta``

Trinity logs and intermediate files in ``trinity_out_dir/``

```bash
module load Trinity/2.15.0

SAMPLE="${SAMPLE}"         
CPU=16
MEM="100G" # high-mem node

Trinity \
  --seqType fq \
  --left  trim_${SAMPLE}_R1.fastq.gz \
  --right trim_${SAMPLE}_R2.fastq.gz \
  --CPU ${CPU} \
  --max_memory ${MEM} \
  --min_contig_length 300 \
  --output trinity_${SAMPLE} \
  --full_cleanup
```

Parameter notes:

``--min_contig_length 300`` — filters very short contigs to reduce noise before downstream steps.

``--max_memory 100G`` — ensures sufficient RAM for in-memory k-mer graphs on large plant datasets.

## Identifying the Coding Regions

We detect long ORFs, annotate PFAM domains on the candidate peptides, and then predict coding regions while retaining PFAM-supported hits.

**Inputs:**

``trinity_combine_${spp}.fasta`` (assembled transcripts)

**Key outputs:**

``${ASSEMBLY}.transdecoder.gff3`` 
``${ASSEMBLY}.transdecoder.bed``

Predicted peptides ``*.transdecoder.pep`` and CDS ``*.transdecoder.cds``

PFAM domain table: ``pfam.domtblout``

**1) Long ORFs scan (TransDecoder.LongOrfs)**

```bash
module load TransDecoder/5.3.0

SPS="${SPS}"             
ASSEMBLY="trinity_combine_${SPS}.fasta"

TransDecoder.LongOrfs -t "${ASSEMBLY}"
```

**2) PFAM domain search on candidate peptides (HMMER hmmscan)**

```bash
module load hmmer/3.2.1

# PFAM HMM library
PFAM_DB="${PFAM_DB}"            
CPU="${SLURM_NTASKS:-16}"

hmmscan --cpu "${CPU}" \
  --domtblout pfam.domtblout \
  "${PFAM_DB}" \
  "${ASSEMBLY}.transdecoder_dir/longest_orfs.pep"
```

**3) Coding region prediction with PFAM retention (TransDecoder.Predict)**

```bash
module load TransDecoder/5.3.0

CPU="${SLURM_NTASKS:-16}"

TransDecoder.Predict -t "${ASSEMBLY}" \
  --retain_pfam_hits pfam.domtblout \
  --cpu "${CPU}"
```

Notes:

``TransDecoder.LongOrfs`` identifies candidate ORFs and writes longest_orfs.pep inside ``${ASSEMBLY}.transdecoder_dir/``

``hmmscan`` populates ``pfam.domtblout``, which is passed to ``TransDecoder.Predict`` via ``--retain_pfam_hits`` to prioritize ORFs with conserved domains.

## Determining and Removing Redundant Transcripts

Cluster predicted coding sequences (CDS) to collapse near-duplicate transcripts and retain one representative (“centroid”) per cluster.

**Inputs:**

``trinity_combine_${SPS}.fasta.transdecoder.cds`` (CDS from TransDecoder)

**Outputs:**

``centroids_${SPS}.fasta`` (non-redundant CDS / “unigenes”)

``clusters.uc`` (cluster membership table)

``LOGFile`` (vsearch run log)

```bash
module load vsearch/2.4.3

SPS="${SPS}"
CPU="${SLURM_NTASKS:-16}"

vsearch \
  --threads "${CPU}" \
  --log LOGFile \
  --cluster_fast trinity_combine_${SPS}.fasta.transdecoder.cds \
  --id 0.90 \
  --centroids centroids_${SPS}.fasta \
  --uc clusters.uc
```
Parameter notes:

``--cluster_fast`` — greedy clustering of nucleotide CDS; sequences are compared in their given orientation.

``--id 0.90`` — merges sequences with ≥90% global identity.

``--centroids`` — writes the representative sequence for each cluster (use as the non-redundant CDS set downstream).

``--uc`` — maps every input sequence to its cluster ID and centroid (traceability).

## Evaluating the Assembly

Assess contiguity/ORF-level metrics with rnaQUAST 1.5.2 (enabling GeneMarkS-T 5.1) and estimate completeness with BUSCO 5.4.5 using the Viridiplantae lineage.

**Inputs:**

``centroids_${SPS}.fasta`` (non-redundant CDS from clustering)

**Outputs:**

``Genemark_${SPS}/`` (rnaQUAST reports, tables, plots; GeneMarkS-T predictions)

``busco_out/`` (BUSCO summary, full tables, single-copy/duplicated/fragmented/missing counts)

rnaQUAST (with GeneMarkS-T)

```bash
module load rnaQUAST/1.5.2
module load GeneMarkS-T/5.1

SPS="${SPS}"
CPU="${SLURM_NTASKS:-16}"

rnaQUAST.py \
  --transcripts centroids_${SPS}.fasta \
  --gene_mark \
  --threads "${CPU}" \
  --output_dir Genemark_${SPS}
```

Notes:

``--gene_mark`` enables GeneMarkS-T on transcripts to report ORF-level metrics.

rnaQUAST summarizes N50/L50, lengths, %GC, and ORF statistics in ``Genemark_${SPS}``.

BUSCO (transcriptome mode; Viridiplantae odb10)

```bash
module load busco/5.4.5

SPS="${SPS}"
CPU="${SLURM_NTASKS:-16}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:viridiplantae_odb10}"

busco \
  -i centroids_${SPS}.fasta \
  -o busco_out \
  -l "${BUSCO_LINEAGE}" \
  -m transcriptome \
  -c "${CPU}" \
  -f
```

Notes

``-m transcriptome`` expects potentially fragmented gene models.

## Quantifying gene expression

Estimate transcript abundances against the non-redundant CDS (“centroids”) using kallisto 0.46.1. This involves (1) building an index from the clustered CDS FASTA and (2) quantifying each paired library against that index. Singletons were not used.

### 1) Creating an index

**Input:**

``centroids_${SPS}.fasta`` (from clustering)

**Output:**

``centroids_${SPS}.fasta.index`` (kallisto index)

```bash
module load kallisto/0.46.1

SPS="${SPS}"

kallisto index -i centroids_${SPS}.fasta.index centroids_${SPS}.fasta
```

### 2) Counting reads mapping to transcripts

**Inputs:**

Trimmed reads: ``trim_${SAMPLE}_R1.fastq.gz``
``trim_${SAMPLE}_R2.fastq.gz``

Index: ``centroids_${SPS}.fasta.index``

**Outputs:**

``abundance.tsv`` 
``abundance.h5``
``run_info.json``

```bash
module load kallisto/0.46.1

SPS="${SPS}"
SAMPLE="${SAMPLE}"
CPU="${SLURM_NTASKS:-16}"

kallisto quant \
  -i centroids_${SPS}.fasta.index \
  -o "${SAMPLE}" \
  -t "${CPU}" \
 trim_${SAMPLE}_R1.fastq.gz \
 trim_${SAMPLE}_R2.fastq.gz
```

## Functional Annotation

Annotate the non-redundant peptides with EnTAP, which orchestrates DIAMOND homology searches against curated protein databases, applies taxonomic contaminant filtering, and transfers functional terms (e.g., GO/orthology) from the best-supported hits according to ``entap_config.ini``.

**Inputs:**

``centroids_${SPS}.pep`` (non-redundant peptide sequences)

``entap_config.ini``

DIAMOND databases: RefSeq (complete proteins), Swiss-Prot, and nr (*.dmnd)

**Outputs:**

EnTAP run directory (e.g., entap_out/) containing summary tables, best-hit assignments, GO/orthology mappings, and logs.

```bash
module load singularity

SPS="${SPS}"
CPU="${SLURM_NTASKS:-16}"

# DIAMOND databases
REFSEQ_DMND="refseq_complete_proteins.dmnd"
SWISSPROT_DMND="uniprot_sprot.dmnd"
NR_DMND="nr_protein.dmnd"

singularity exec entap.sif EnTAP \
  --runP \
  --entap-ini entap_config.ini \
  --threads "${CPU}" \
  --input "centroids_${SPS}.pep" \
  -d "${REFSEQ_DMND}" \
  -d "${SWISSPROT_DMND}" \
  -d "${NR_DMND}"
```

Parameter notes:

``--runP`` runs the full pipeline (filtering → DIAMOND searches → functional assignment).

``--entap-ini`` controls coverage/e-value thresholds, taxonomy, and optional expression filters (e.g., FPKM ≥ 0.5).

Multiple ``-d`` flags define the search order; EnTAP selects the best supported hit considering alignment quality and taxonomic relevance.
