---
title: "Processing Stool Sample Whole Genome Sequencing Data"
author: "Ankan Choudhury"
date: "2025-03-05"
output:
  html_document:
    df_print: paged
header-includes:
- \usepackage{times}
- \usepackage{listings}
- \lstset{breaklines=true}
- \usepackage{fvextra}
- \DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\{\}}
---

```{r setup, echo=FALSE}
knitr::opts_chunk$set(warning = FALSE, message=FALSE, tidy.opts=list(width.cutoff=50),tidy=TRUE)
```

```{r packages_deafult,echo=FALSE}

library(biomformat)
library(tidyverse)
library(data.table)
library(fs)
library(phyloseq)
library(writexl)
library(phylostratr)
library(ape)
library(taxize)
```

# **Introduction**

This code includes the bash code for performing Whole-genome Sequence analysis on the Stool samples using Trimmomatic, Bowtie, Kraken, Bracken, Kaiju, Metaphlan, and Humann.The paired sequences ("raw") were first trimmed using Trimmomatic v0.39, which trimmed the sequences and separated out into paired ("trimmed") and unpaired, sequences who do not have a matching read pairing with it. The trimmed sequences were then filtered for host reads (*Homo sapiens*) using Bowtie2 resultant sequences ("trimmed-filtered") were used for taxonomic classification. Taxonomic classification was done at all three levels of sequence refinement were done using Kraken2, Kaiju, and Metaphlan. Kraken2 classification used standard Kraken database and was normalized by Bracken using standard Bracken database. Kaiju classification was done using the Kaiju RefSeq database. Metaphlan 4.1.1 used its own ChocoPhlan database (updated October 2024). The trimmed-filtered sequences were also assembled into contigs using Megahit, and the assembled sequences ("trimmed-filtered-contigs") were also classified using Kraken2-Bracken algorithm. The raw, trimmed, trimmed-filtered sequences were also functionally annotated using Humann 3.9. The resultant taxonomic classification data and functional annotation data were transformed and compiled into phyloseq objects in R.

## **Prerequisites**

 1. Install qiime2-metagenome-2024.10
 2. Install Trimmomatic 0.39
 3. Install Metaphlan 4.4.1
 4. Install Humann 3.9
 5. Install R 4.3.1 or higher

## **1. Getting the sequences**

All the code chunks given below should be made into a .sh files and executed by batch submission in kodiak server as they usually take long time to run.

```{bash,eval=FALSE,echo=TRUE}
ZIP_FILE="061824KGmetagenome.zip"  # The zip file containing the sequences
EXTRACT_DIR="extracted_files"
DEST_DIR="demultiplexed/raw-paired"
mkdir -p "$EXTRACT_DIR" "$DEST_DIR"

unzip -o "$ZIP_FILE" -d "$EXTRACT_DIR"

find "$EXTRACT_DIR" -type f -name "*.fastq.gz" -exec mv {} "$DEST_DIR" \;
```

For subsequent operations, please use a custom temporary folder when using kodiak server as if the temporary files generated during any operation exceeds the given space, it will shut the operation down. Add this piece of code to your bash script every time

```{bash,eval=FALSE,echo=TRUE}
TEMP_DIR="tmp"
mkdir -p "$TEMP_DIR"
export TMPDIR="$TEMP_DIR"
export TMP="$TEMP_DIR"
export TEMP="$TEMP_DIR"
```

## **2. Trimming the sequences using Trimmomatic v0.39**

Install Trimmomatic v0.39 prior to this. Include the path where trimmomatic-0.39.jar executable file is stored.

```{bash,eval=FALSE,echo=TRUE}

cd demultiplexed 

mkdir -p trimmed-paired # to store trimmed and pair matched sequences
mkdir -p trimmed-unpaired # to store unmatched sequences

for f in $(find raw-paired -maxdepth 1 -name "*R1_001.fastq.gz" -printf "%f\n" | sed 's/1_001.fastq.gz//' | sort -u)
do 
java -jar /path/to/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 ${f}1_001.fastq.gz ${f}2_001.fastq.gz trimmed-paired/${f}1_001.fastq.gz trimmed-unpaired/${f}1_001.fastq.gz trimmed-paired/${f}2_001.fastq.gz trimmed-unpaired/${f}2_001.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35;
echo Processed: ${f}1_001.fastq.gz & ${f}2_001.fastq.gz
done

```

Storing the trimmed sequences as a qiime2 artifact (.qza)

```{bash,eval=FALSE,echo=TRUE}

source activate qiime2-metagenome-2024.10 #activating qiime2 metagenome distribution

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path demultiplexed/trimmed-paired \   
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired/trimmed.qza
  
qiime demux summarize \
  --i-data demux-paired/trimmed.qza \
  --o-visualization demux-paired/trimmed.qzv  #summary visualization of the sequences
```

## **3. Filtering the sequences with Bowtie2**

Bowtie2 requires a Bowtie database built from genomic data provided by user. Currently, the database used here was built from Human Genome Data gathered from NCBI. If you are using Bowtie2 standalone installation then you can download it from the Bowtie indices for H. sapiens, GRCh38 no-alt analysis set from <a href="https://bowtie-bio.sourceforge.net/bowtie2/index.shtml">here</a>. Here Bowtie2 was used from the qiime2 metagenome distribution platform, so the database has to be in a qimme2 artifact (.qza) format. For that, we can download the FASTA version of H. sapiens genome GRCh38 in NCBI from <a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/">here</a>, and build it into Bowtie2 datbase format using the following bash script.

```{bash,eval=FALSE,echo=TRUE}
cd /path/to/where/it/is/stored/GRCh38_FASTA

source activate qiime2-metagenome-2024.10

qiime tools import \                #Importing the genome stored in .fna to .qza format
  --input-path GRCh38.fna \
  --output-path GRCh38.qza \
  --type 'FeatureData[Sequence]'

qiime quality-control bowtie2-build \   #Building Bowtie2 Index from it
	--i-sequences GRCh38.qza \
	--o-database Hsapiens_BowTie2Index.qza

```

After BowTie2 Indices are built, we can use that to filter human reads from the sequences using the following bash script.

```{bash,eval=FALSE,echo=TRUE}
source activate qiime2-metagenome-2024.10

qiime quality-control filter-reads \
	--i-demultiplexed-sequences demux-paired/trimmed.qza \           #Input trimmed sequences 
	--i-database /path/to/where/it/is/stored/Hsapiens_BowTie2Index.qza \        #Bowtie2 Index built from Human genome
	--o-filtered-sequences demux-paired/trimmed-filtered.qza        #Output trimmed & filtered sequences
```

## **4. Assembling the sequences into contigs using Megahit**

The assemble contigs can be used for taxonomic classification using Kraken2 and functional annotation using eggNOG.
```{bash,eval=FALSE,echo=TRUE}
mkdir -p contigs 

source activate qiime2-metagenome-2024.10

qiime assembly assemble-megahit \
    --i-seqs demux-paired/trimmed-filtered.qza \        #Input trimmed and filtered sequences
    --p-presets "meta-sensitive" \        #Preset which looks for more rarer species in the metagenome
    --p-num-cpu-threads 24 \        #number of cpus to be used
    --p-num-partitions 1 \        
    --o-contigs contigs/trimmed-filtered-contigs.qza \
    --verbose
```

## **5. Taxonomic classification using Kraken2**

We can use Kraken2 to classify both sequences and contigs using the same basic script. Kraken2 is a taxonomic classification system using exact k-mer matches to achieve high accuracy and fast classification speeds. This classifier matches each k-mer within a query sequence to the lowest common ancestor (LCA) of all genomes containing the given k-mer. You can install Kraken2 as a stand-alone, but here we have used kraken2 through the qiime2 metagenome platform's moshpit distribution package. Before that we will have to build databses for Kraken2 and Bracken using this script. We will be storing the database in the qiime cache. A cache will store the "data" of a dataset in a separate data folder with a random generated name, which is stored in a separate "key" file named after the dataset in  the keys folder. In a qiime artifact (.qza), they are all bundled up together into one thing that has to zipped and unzipped everytime you use it.

```{bash,eval=FALSE,echo=TRUE}
source activate qiime2-metagenome-2024.10

qiime tools cache-create --cache cache        #creating the cache

qiime moshpit build-kraken-db \
    --p-collection standard \       #Using the standard Kraken2 and Bracken database
    --o-kraken2-database cache:kraken_standard \
    --o-bracken-database cache:bracken_standard \
    --verbose

```

After the databases are built we can go ahead and classify the data using Kraken2.

```{bash,eval=FALSE,echo=TRUE}
source activate qiime2-metagenome-2024.10

qiime moshpit classify-kraken2 \
	--i-seqs contigs/trimmed-filtered.qza \       #Input sequences/contigs
	--i-kraken2-db cache:kraken_standard \
	--p-threads 40 \        #Number of threads to use in computation
	--p-confidence 0.1 \       #Confidence level (0 to 1) of the matches that will be made, higher the value more stringent the algorithm
	--p-minimum-base-quality 0 \        #Minimum base quality at which Kraken2 classifies a read
	--o-hits cache:kraken_hits-trimmed-filtered \       #Kraken2 hits stored in cache format
	--o-reports cache:kraken_reports-trimmed-filtered \       #Kraken2 reports stored in cache format
	--use-cache cache \
	--verbose \
  --p-memory-mapping False

```

Kraken outputs conatin a feature table (.biom), taxnomy file (.txt) and kraken reports (.txt). Kraken reports contain taxon information on each line: Percentage of fragments covered by the clade root, number of fragments covered by clade root, Number of fragments assigned directly to this taxon, a rank code: indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies, NCBI taxonomic ID number, and taxonomic annotation.

```{r Kraken2Report, echo=FALSE}
kraken<-readLines("Kraken Outputs/New/kraken_reports-trimmed-filtered/1-MDcol002-5-2.report.txt")
writeLines(kraken[1:25])

```
### **Renormalization with Bracken**

Bracken (Bayesian Reestimation of Abundance with KrakEN) is a highly accurate statistical method that computes the abundance of species in DNA sequences from a metagenomics sample. Bracken is a companion program to Kraken 2. While Kraken classifies reads to multiple levels in the taxonomic tree, Bracken allows estimation of abundance at a single level using those classifications (e.g. Bracken can estimate abundance of species within a sample). Bracken uses a Bracken database, the length of your reads and the kraken reports to give you a feature frequency table that renormalizes the data by dropping the unclassified reads and according to the Bayesian probability of the kraken hits being correct.

```{bash,eval=FALSE,echo=TRUE}
source activate qiime2-metagenome-2024.10

mkdir kraken-outputs

name="trimmed-filtered"

qiime moshpit estimate-bracken \
    --i-bracken-db cache:bracken_standard \
    --p-read-len 150 \
    --i-kraken-reports cache:kraken_reports-$name \
    --o-reports kraken-outputs/bracken-reports-$name.qza \
    --o-taxonomy kraken-outputs/taxonomy-bracken-$name.qza \
    --o-table kraken-outputs/table-bracken-$name.qza


```

The Bracken reports are very similar to Kraken2 reports except for having no unclassified reads
```{r BrackenTaxonomy, echo=FALSE}
bracken<-readLines("BrackenOutputs/bracken-reports-trimmed-filtered/1-MDcol002-5-2.report.txt")
writeLines(bracken[1:20])

```

The taxonomy file contains the NCBI taxa ID number and the string of taxonomic names leading to the species or the final taxonomic ranks for that ID.

```{r BrackenReport, echo=FALSE}
brackentaxonomy<-readLines("BrackenOutputs/taxonomy-bracken-trimmed-filtered/taxonomy.tsv")
writeLines(brackentaxonomy[1:20])

```

The feature table (.biom) is a BIOM format matrix table with Feature IDs as the row number and the numbe rof reads for that feature in each samples as column data.

```{r BrackenTable, echo=FALSE}
brackentable<-as.data.frame(as.matrix(biomformat::biom_data(read_biom("BrackenOutputs/table-bracken-trimmed-filtered/feature-table.biom"))))
head(brackentable)

```
### **Creating phylogenetic trees from the bracken data**

We can generate phylogenetic trees from the brcaken data using a package called gracken, which can be installed by ```pip install gracken``` or from this link (https://github.com/jonasoh/gracken). But gracken uses scientific names as phylogenetic tree tip labels whereas here the feature table contains the NCBI taxa IDs, so you can replace the executable python file for gracken in the miniconda path "miniconda/lib/python3.9/site-packages/gracken/" with this file <a href="data:application/octet-stream;base64,aW1wb3J0IG9zCmltcG9ydCBzeXMKaW1wb3J0IGdsb2IKaW1wb3J0IGFyZ3BhcnNlCmltcG9ydCBwYW5kYXMgYXMgcGQKZnJvbSAuIGltcG9ydCBfX3ZlcnNpb25fXwpmcm9tIGV0ZTMgaW1wb3J0IE5DQklUYXhhCmZyb20gLnZhcnMgaW1wb3J0IHRheF9jb2xzCmZyb20gZXRlMyBpbXBvcnQgVHJlZSBhcyBUcmVlMwpmcm9tIC4gaW1wb3J0IG5jYmlfdXRpbHMgYXMgbmMKZnJvbSAuIGltcG9ydCBndGRiX3V0aWxzIGFzIGd0CmZyb20gLmRhdGFfbG9hZGVycyBpbXBvcnQgcmVhZF9icmFja2VuX3JlcG9ydCwgcmVhZF9rcmFrZW4yX3JlcG9ydAoKcGQub3B0aW9ucy5tb2RlLmNvcHlfb25fd3JpdGUgPSBUcnVlCgoKZGVmIHBhcnNlX2FyZ3MoKToKICAgIHBhcnNlciA9IGFyZ3BhcnNlLkFyZ3VtZW50UGFyc2VyKAogICAgICAgIGRlc2NyaXB0aW9uPSJDcmVhdGVzIGEgcGh5bG9nZW5ldGljIHRyZWUgYW5kIE9UVSB0YWJsZSBmcm9tIEJyYWNrZW4vS3Jha2VuMiByZXBvcnRzIGJ5IHBydW5pbmcgR1REQi9OQ0JJIHRyZWVzIgogICAgKQogICAgcGFyc2VyLmFkZF9hcmd1bWVudCgKICAgICAgICAiLS12ZXJzaW9uIiwKICAgICAgICAiLXYiLAogICAgICAgIGFjdGlvbj0idmVyc2lvbiIsCiAgICAgICAgdmVyc2lvbj1fX3ZlcnNpb25fXywKICAgICAgICBoZWxwPSJzaG93IHByb2dyYW0gdmVyc2lvbiBhbmQgZXhpdCIsCiAgICApCiAgICBwYXJzZXIuYWRkX2FyZ3VtZW50KAogICAgICAgICItLWlucHV0X2RpciIsCiAgICAgICAgIi1pIiwKICAgICAgICByZXF1aXJlZD1UcnVlLAogICAgICAgIGhlbHA9ImRpcmVjdG9yeSBjb250YWluaW5nIEJyYWNrZW4vS3Jha2VuMiByZXBvcnQgZmlsZXMiLAogICAgKQogICAgcGFyc2VyLmFkZF9hcmd1bWVudCgKICAgICAgICAiLS1iYWNfdGF4b25vbXkiLAogICAgICAgIGhlbHA9InBhdGggdG8gR1REQiBiYWN0ZXJpYWwgdGF4b25vbXkgZmlsZSIsCiAgICApCiAgICBwYXJzZXIuYWRkX2FyZ3VtZW50KAogICAgICAgICItLWFyX3RheG9ub215IiwKICAgICAgICBoZWxwPSJwYXRoIHRvIEdUREIgYXJjaGFlYWwgdGF4b25vbXkgZmlsZSIsCiAgICApCiAgICBwYXJzZXIuYWRkX2FyZ3VtZW50KAogICAgICAgICItLWJhY190cmVlIiwKICAgICAgICBoZWxwPSJwYXRoIHRvIEdUREIgYmFjdGVyaWFsIHRyZWUgZmlsZSIsCiAgICApCiAgICBwYXJzZXIuYWRkX2FyZ3VtZW50KAogICAgICAgICItLWFyX3RyZWUiLAogICAgICAgIGhlbHA9InBhdGggdG8gR1REQiBhcmNoYWVhbCB0cmVlIGZpbGUiLAogICAgKQogICAgcGFyc2VyLmFkZF9hcmd1bWVudCgKICAgICAgICAiLS1vdXRfcHJlZml4IiwKICAgICAgICAiLW8iLAogICAgICAgIGRlZmF1bHQ9Im91dHB1dCIsCiAgICAgICAgaGVscD0icHJlZml4IGZvciBvdXRwdXQgZmlsZXMgKGRlZmF1bHQ6IG91dHB1dCkuIENyZWF0ZXMgPHByZWZpeD4udHJlZSBhbmQgPHByZWZpeD4ub3R1LmNzdiIsCiAgICApCiAgICBwYXJzZXIuYWRkX2FyZ3VtZW50KAogICAgICAgICItLW1vZGUiLAogICAgICAgICItbSIsCiAgICAgICAgY2hvaWNlcz1bImJyYWNrZW4iLCAia3Jha2VuMiJdLAogICAgICAgIGRlZmF1bHQ9ImJyYWNrZW4iLAogICAgICAgIGhlbHA9ImlucHV0IGZpbGUgZm9ybWF0IChkZWZhdWx0OiBicmFja2VuKSIsCiAgICApCiAgICBwYXJzZXIuYWRkX2FyZ3VtZW50KAogICAgICAgICItLWtlZXAtc3BhY2VzIiwKICAgICAgICBhY3Rpb249InN0b3JlX3RydWUiLAogICAgICAgIGRlZmF1bHQ9RmFsc2UsCiAgICAgICAgaGVscD0ia2VlcCBzcGFjZXMgaW4gc3BlY2llcyBuYW1lcyAoZGVmYXVsdDogRmFsc2UpIiwKICAgICkKICAgIHBhcnNlci5hZGRfYXJndW1lbnQoCiAgICAgICAgIi0tdGF4b25vbXkiLAogICAgICAgICItdCIsCiAgICAgICAgY2hvaWNlcz1bImd0ZGIiLCAibmNiaSJdLAogICAgICAgIGRlZmF1bHQ9Im5jYmkiLAogICAgICAgIGhlbHA9InRheG9ub215IHNvdXJjZSB0byB1c2UgKGRlZmF1bHQ6IG5jYmkpLiIsCiAgICApCiAgICBwYXJzZXIuYWRkX2FyZ3VtZW50KAogICAgICAgICItLWZ1bGwtdGF4b25vbXkiLAogICAgICAgICItZiIsCiAgICAgICAgYWN0aW9uPSJzdG9yZV90cnVlIiwKICAgICAgICBkZWZhdWx0PUZhbHNlLAogICAgICAgIGhlbHA9ImluY2x1ZGUgZnVsbCB0YXhvbm9teSBpbmZvIGluIE9UVSB0YWJsZSAoZGVmYXVsdDogRmFsc2UpIiwKICAgICkKICAgIGFyZ3MgPSBwYXJzZXIucGFyc2VfYXJncygpCiAgICBpZiBhcmdzLnRheG9ub215ID09ICJndGRiIjoKICAgICAgICBtaXNzaW5nX2FyZ3MgPSBbXQogICAgICAgIGlmIG5vdCBhcmdzLmJhY190YXhvbm9teToKICAgICAgICAgICAgbWlzc2luZ19hcmdzLmFwcGVuZCgiLS1iYWNfdGF4b25vbXkiKQogICAgICAgIGlmIG5vdCBhcmdzLmFyX3RheG9ub215OgogICAgICAgICAgICBtaXNzaW5nX2FyZ3MuYXBwZW5kKCItLWFyX3RheG9ub215IikKICAgICAgICBpZiBub3QgYXJncy5iYWNfdHJlZToKICAgICAgICAgICAgbWlzc2luZ19hcmdzLmFwcGVuZCgiLS1iYWNfdHJlZSIpCiAgICAgICAgaWYgbm90IGFyZ3MuYXJfdHJlZToKICAgICAgICAgICAgbWlzc2luZ19hcmdzLmFwcGVuZCgiLS1hcl90cmVlIikKICAgICAgICBpZiBtaXNzaW5nX2FyZ3M6CiAgICAgICAgICAgIHBhcnNlci5lcnJvcigKICAgICAgICAgICAgICAgICJGb3IgR1REQiB0YXhvbm9teSwgdGhlIGZvbGxvd2luZyBhcmd1bWVudHMgYXJlIHJlcXVpcmVkOiAiCiAgICAgICAgICAgICAgICArICIsICIuam9pbihtaXNzaW5nX2FyZ3MpCiAgICAgICAgICAgICkKICAgIHJldHVybiBhcmdzCgoKZGVmIG1haW4oKToKICAgIGFyZ3MgPSBwYXJzZV9hcmdzKCkKCiAgICAjIGxvYWQgZGF0YSBhbmQgY3JlYXRlIE9UVSB0YWJsZQogICAgb3R1X3RhYmxlID0gcGQuRGF0YUZyYW1lKCkKICAgIGZpbGVfcGF0dGVybiA9ICIqLnJlcG9ydC50eHQiCgogICAgbWF0Y2hpbmdfZmlsZXMgPSBnbG9iLmdsb2Iob3MucGF0aC5qb2luKGFyZ3MuaW5wdXRfZGlyLCBmaWxlX3BhdHRlcm4pKQogICAgaWYgbm90IG1hdGNoaW5nX2ZpbGVzOgogICAgICAgIHByaW50KGYiRXJyb3I6IE5vIHthcmdzLm1vZGV9IGZpbGVzIGZvdW5kIGluIHthcmdzLmlucHV0X2Rpcn0iLCBmaWxlPXN5cy5zdGRlcnIpCiAgICAgICAgc3lzLmV4aXQoMSkKCiAgICBmb3IgZiBpbiBtYXRjaGluZ19maWxlczoKICAgICAgICBuYW1lID0gb3MucGF0aC5zcGxpdGV4dChvcy5wYXRoLmJhc2VuYW1lKGYpKVswXQogICAgICAgIGlmIGFyZ3MubW9kZSA9PSAiYnJhY2tlbiI6CiAgICAgICAgICAgIHNhbXBsZV9vdHUgPSByZWFkX2JyYWNrZW5fcmVwb3J0KGYpCiAgICAgICAgZWxzZToKICAgICAgICAgICAgc2FtcGxlX290dSA9IHJlYWRfa3Jha2VuMl9yZXBvcnQoZikKCiAgICAgICAgIyBDaGVjayBpZiBzYW1wbGVfb3R1IGlzIGVtcHR5IG9yIGRvZXNuJ3QgaGF2ZSB0aGUgcmVxdWlyZWQgY29sdW1ucwogICAgICAgIGlmIHNhbXBsZV9vdHUuZW1wdHkgb3Igbm90IGFsbCgKICAgICAgICAgICAgY29sIGluIHNhbXBsZV9vdHUuY29sdW1ucwogICAgICAgICAgICBmb3IgY29sIGluICgKICAgICAgICAgICAgICAgIFsidGF4aWQiLCAiYWJ1bmRhbmNlIl0KICAgICAgICAgICAgICAgIGlmIGFyZ3MudGF4b25vbXkgPT0gIm5jYmkiCiAgICAgICAgICAgICAgICBlbHNlIFsic3BlY2llcyIsICJhYnVuZGFuY2UiXQogICAgICAgICAgICApCiAgICAgICAgKToKICAgICAgICAgICAgcHJpbnQoCiAgICAgICAgICAgICAgICBmIldhcm5pbmc6IFNraXBwaW5nIHtmfSBiZWNhdXNlIGl0J3MgZW1wdHkgb3IgbWFsZm9ybWVkLiIsCiAgICAgICAgICAgICAgICBmaWxlPXN5cy5zdGRlcnIsCiAgICAgICAgICAgICkKICAgICAgICAgICAgY29udGludWUKCiAgICAgICAgc2FtcGxlX290dVsic2FtcGxlIl0gPSBuYW1lCiAgICAgICAgb3R1X3RhYmxlID0gcGQuY29uY2F0KFtvdHVfdGFibGUsIHNhbXBsZV9vdHVdLCBpZ25vcmVfaW5kZXg9VHJ1ZSkKCiAgICBkZWYgZ2V0X3NhbXBsZV9jb2xzKGNvbHMpOgogICAgICAgIHJldHVybiBbYyBmb3IgYyBpbiBjb2xzIGlmIGMgbm90IGluIHRheF9jb2xzXQoKICAgIGlmIGFyZ3MudGF4b25vbXkgPT0gIm5jYmkiOgogICAgICAgIG5jYmkgPSBOQ0JJVGF4YSgpCgogICAgICAgICMgR2V0IHVuaXF1ZSB0YXhpZHMgYW5kIHRyYW5zbGF0ZSB0byBzcGVjaWVzIG5hbWVzCiAgICAgICAgdGF4aWRfbGlzdCA9IGxpc3Qob3R1X3RhYmxlWyJ0YXhpZCJdLnVuaXF1ZSgpKQogICAgICAgIHRyYW5zbGF0b3IgPSBuY2JpLmdldF90YXhpZF90cmFuc2xhdG9yKHRheGlkX2xpc3QpCiAgICAgICAgb3R1X3RhYmxlWyJzcGVjaWVzIl0gPSBvdHVfdGFibGVbInRheGlkIl0ubWFwKAogICAgICAgICAgICBsYW1iZGEgdGlkOiB0cmFuc2xhdG9yLmdldChpbnQodGlkKSwgc3RyKHRpZCkpCiAgICAgICAgKQoKICAgICAgICBpZiBub3QgYXJncy5rZWVwX3NwYWNlczoKICAgICAgICAgICAgb3R1X3RhYmxlWyJzcGVjaWVzIl0gPSBvdHVfdGFibGVbInNwZWNpZXMiXS5zdHIucmVwbGFjZSgiICIsICJfIikKCiAgICAgICAgd2lkZV9vdHVfdGFibGUgPSBvdHVfdGFibGUucGl2b3RfdGFibGUoCiAgICAgICAgICAgIGluZGV4PSJzcGVjaWVzIiwgY29sdW1ucz0ic2FtcGxlIiwgdmFsdWVzPSJhYnVuZGFuY2UiCiAgICAgICAgKS5yZXNldF9pbmRleCgpCgogICAgICAgIGlmIGFyZ3MuZnVsbF90YXhvbm9teToKICAgICAgICAgICAgbWFwcGluZyA9IG5jLmJ1aWxkX2xpbmVhZ2VfbWFwcGluZyhuY2JpLCBvdHVfdGFibGUpCiAgICAgICAgICAgIGZvciBjb2wgaW4gdGF4X2NvbHNbOi0xXTogICMgZXhjbHVkZSBzcGVjaWVzCiAgICAgICAgICAgICAgICB3aWRlX290dV90YWJsZVtjb2xdID0gd2lkZV9vdHVfdGFibGVbInNwZWNpZXMiXS5tYXAoCiAgICAgICAgICAgICAgICAgICAgbGFtYmRhIHNwOiBtYXBwaW5nLmdldChzcCwge30pLmdldChjb2wsICIiKQogICAgICAgICAgICAgICAgKQogICAgICAgICAgICBzYW1wbGVfY29scyA9IGdldF9zYW1wbGVfY29scyh3aWRlX290dV90YWJsZS5jb2x1bW5zKQogICAgICAgICAgICB3aWRlX290dV90YWJsZSA9IHdpZGVfb3R1X3RhYmxlW3RheF9jb2xzICsgc2FtcGxlX2NvbHNdCgogICAgICAgICMgQnVpbGQgdHJlZSB1c2luZyBOQ0JJIHRheG9ub215IGFuZCB1cGRhdGUgbGVhZiBuYW1lcwogICAgICAgIHRyZWUgPSBuY2JpLmdldF90b3BvbG9neSh0YXhpZF9saXN0KQogICAgICAgICNmb3IgbGVhZiBpbiB0cmVlLmdldF9sZWF2ZXMoKToKICAgICAgICAgICAgI2xlYWYubmFtZSA9IHRyYW5zbGF0b3IuZ2V0KGludChsZWFmLm5hbWUpLCBsZWFmLm5hbWUpCiAgICAgICAgICAgICNpZiBub3QgYXJncy5rZWVwX3NwYWNlczoKICAgICAgICAgICAgICAgICNsZWFmLm5hbWUgPSBsZWFmLm5hbWUucmVwbGFjZSgiICIsICJfIikKCiAgICAjIEdUREIgdGF4b25vbXkKICAgIGVsc2U6CiAgICAgICAgd2lkZV9vdHVfdGFibGUgPSBvdHVfdGFibGUucGl2b3RfdGFibGUoCiAgICAgICAgICAgIGluZGV4PSJzcGVjaWVzIiwgY29sdW1ucz0ic2FtcGxlIiwgdmFsdWVzPSJhYnVuZGFuY2UiCiAgICAgICAgKS5yZXNldF9pbmRleCgpCgogICAgICAgICMgd2UgZHJvcCB0aGUgc3BlY2llcyBwcmVmaXggaGVyZSBzaW5jZSBzb21lIGRhdGFiYXNlcyBoYXZlIGl0IGFuZCBvdGhlcnMgZG9uJ3QKICAgICAgICB3aWRlX290dV90YWJsZVsic3BlY2llcyJdID0gd2lkZV9vdHVfdGFibGVbInNwZWNpZXMiXS5zdHIucmVwbGFjZSgKICAgICAgICAgICAgciJec19fIiwgIiIsIHJlZ2V4PVRydWUKICAgICAgICApCgogICAgICAgIGlmIG5vdCBhcmdzLmtlZXBfc3BhY2VzOgogICAgICAgICAgICB3aWRlX290dV90YWJsZVsic3BlY2llcyJdID0gd2lkZV9vdHVfdGFibGVbInNwZWNpZXMiXS5zdHIucmVwbGFjZSgiICIsICJfIikKCiAgICAgICAgaWYgYXJncy5mdWxsX3RheG9ub215OgogICAgICAgICAgICBtYXBwaW5nID0gZ3QuYnVpbGRfdGF4b25vbXlfbWFwcGluZygKICAgICAgICAgICAgICAgIFthcmdzLmJhY190YXhvbm9teSwgYXJncy5hcl90YXhvbm9teV0sIGFyZ3Mua2VlcF9zcGFjZXMKICAgICAgICAgICAgKQogICAgICAgICAgICBmb3IgY29sIGluIHRheF9jb2xzWzotMV06ICAjIGV4Y2x1ZGUgc3BlY2llcwogICAgICAgICAgICAgICAgd2lkZV9vdHVfdGFibGVbY29sXSA9IHdpZGVfb3R1X3RhYmxlWyJzcGVjaWVzIl0ubWFwKAogICAgICAgICAgICAgICAgICAgIGxhbWJkYSBzcDogbWFwcGluZy5nZXQoc3AsIHt9KS5nZXQoY29sLCAiIikKICAgICAgICAgICAgICAgICkKICAgICAgICAgICAgc2FtcGxlX2NvbHMgPSBnZXRfc2FtcGxlX2NvbHMod2lkZV9vdHVfdGFibGUuY29sdW1ucykKCiAgICAgICAgICAgIHdpZGVfb3R1X3RhYmxlID0gd2lkZV9vdHVfdGFibGVbdGF4X2NvbHMgKyBzYW1wbGVfY29sc10KCiAgICAgICAgIyB1bmlxdWUgc3BlY2llcwogICAgICAgIHNwZWNpZXNfc2V0ID0gc2V0KG90dV90YWJsZVsic3BlY2llcyJdKQoKICAgICAgICAjIGxvYWQgdGF4b25vbWllcwogICAgICAgIGJhY19zcGVjaWVzX3RvX2dlbm9tZXMsIGJhY19nZW5vbWVfdG9fc3BlY2llcyA9IGd0LnByb2Nlc3NfdGF4b25vbXkoCiAgICAgICAgICAgIGFyZ3MuYmFjX3RheG9ub215CiAgICAgICAgKQogICAgICAgIGFyX3NwZWNpZXNfdG9fZ2Vub21lcywgYXJfZ2Vub21lX3RvX3NwZWNpZXMgPSBndC5wcm9jZXNzX3RheG9ub215KAogICAgICAgICAgICBhcmdzLmFyX3RheG9ub215CiAgICAgICAgKQoKICAgICAgICAjIGFkZCBiYWNrIHNwZWNpZXMgcHJlZml4IGZvciBzZWFyY2hpbmcgaW4gR1REQiB0YXhvbm9teSBmaWxlCiAgICAgICAgbXlfc3BlY2llcyA9IHsic19fIiArIHggaWYgbm90IHguc3RhcnRzd2l0aCgic19fIikgZWxzZSB4IGZvciB4IGluIHNwZWNpZXNfc2V0fQoKICAgICAgICBiYWNfZm91bmQsIF8gPSBndC5maW5kX21hdGNoaW5nX3NwZWNpZXMobXlfc3BlY2llcywgYmFjX3NwZWNpZXNfdG9fZ2Vub21lcykKICAgICAgICBhcl9mb3VuZCwgXyA9IGd0LmZpbmRfbWF0Y2hpbmdfc3BlY2llcyhteV9zcGVjaWVzLCBhcl9zcGVjaWVzX3RvX2dlbm9tZXMpCgogICAgICAgIGJhY190cmVlID0gZ3QucHJvY2Vzc190cmVlKGFyZ3MuYmFjX3RyZWUsIGJhY19nZW5vbWVfdG9fc3BlY2llcywgYmFjX2ZvdW5kKQogICAgICAgIGFyX3RyZWUgPSBndC5wcm9jZXNzX3RyZWUoYXJncy5hcl90cmVlLCBhcl9nZW5vbWVfdG9fc3BlY2llcywgYXJfZm91bmQpCgogICAgICAgIHRyZWUgPSBUcmVlMygpCiAgICAgICAgdHJlZS5uYW1lID0gInJvb3QiCiAgICAgICAgdHJlZS5hZGRfY2hpbGQoYmFjX3RyZWUpCiAgICAgICAgdHJlZS5hZGRfY2hpbGQoYXJfdHJlZSkKCiAgICB3aWRlX290dV90YWJsZS50b19jc3YoZiJ7YXJncy5vdXRfcHJlZml4fS5vdHUuY3N2Iiwgc2VwPSIsIiwgaW5kZXg9RmFsc2UpCiAgICB0cmVlLndyaXRlKG91dGZpbGU9ZiJ7YXJncy5vdXRfcHJlZml4fS50cmVlIiwgZm9ybWF0PTEsIHF1b3RlZF9ub2RlX25hbWVzPUZhbHNlKQoKCmlmIF9fbmFtZV9fID09ICJfX21haW5fXyI6CiAgICBtYWluKCkK" download="gracken.py">gracken.py</a>, which contains an updated code. Once you have it installed, this code will create the phylogentic tree for each dataset ("raw", "trimmed","filtered"...etc) and an OTU table , but we can ignore that.

```{bash,eval=FALSE,echo=TRUE}
mkdir -p bracken-outputs

source activate qiime2-metagenome-2024.10

for f in $(find kraken-outputs -maxdepth 2 -name "table-bracken*")        #Converting the bracken feature tables from a .qza to folder format
do 
filename=$(basename "$f" .qza)
qiime tools export \
	--input-path $f \
	--output-path bracken-outputs/$filename
echo $filename
done

for f in $(find kraken-outputs -maxdepth 2 -name "taxonomy-bracken*")       #Converting the bracken taxonomy from a .qza to folder format
do 
filename=$(basename "$f" .qza)
qiime tools export \
	--input-path $f \
	--output-path bracken-outputs/$filename
echo $filename
done

for f in $(find kraken-outputs -maxdepth 2 -name "bracken-reports*")        #Converting the bracken reports from a .qza to folder format
do 
filename=$(basename "$f" .qza)
qiime tools export \
	--input-path $f \
	--output-path bracken-outputs/$filename
echo $filename
done

find bracken-outputs -type d -name "bracken-reports-*" | while read dir; do #Create a phylogenetic tree for each dataset using the bracken reports
     gracken --input_dir $dir --out_prefix $dir/output --mode bracken
done


```

## **6. Taxonomic Classification using Kaiju**

Kaiju is another taxonomic classification tool that assigns each sequencing read to a taxon in the NCBI taxonomy by comparing it to a reference protein database containing microbial and viral protein sequences. By using protein-level classification, Kaiju achieves a higher sensitivity compared with methods based on nucleotide comparison. You can use any of the several reference protein databases, such as complete genomes from NCBI RefSeq or the microbial subset of the NCBI BLAST non-redundant protein database nr, optionally also including fungi and microbial eukaryotes. We have used the NCBI RefSeq and have run the kaiju classifation through the qiime2 metagenome platform's moshpit distribution package. The databses can be built in the following manner and stored as a cache artifact.

```{bash,eval=FALSE,echo=TRUE}
source activate qiime2-metagenome-2024.10
qiime moshpit fetch-kaiju-db \
	 --p-database-type 'refseq' \
	 --o-database cache:KaijuDB_refseq \
	 --verbose
```

After the databasehas been buit we can classify the datasets by Kaiju to any taxonomic level of choice. Here we have classiifed it to the species level, hence any reads that could not be assigned to any species level with enough confidence but can be assigned to other taxonomic ranks will be put in the "cannot be assigned to a (non-viral) species" category, which is different than the "unclassified" category, which contains reads that cannot be assigned to any taxonomic category at all. To classify using Kaiju, use the following code.

```{bash,eval=FALSE,echo=TRUE}
source activate qiime2-metagenome-2024.10
qiime moshpit classify-kaiju \
	--i-seqs /data/choudhurya/COLON_MD/contigs/filtered-contigs.qza \
	--i-db cache:KaijuDB_refseq \
	--p-z 40 \
	--p-a 'greedy' \
	--p-c 0 \
	--o-abundances cache:kaiju_abun-trimmed-filtered \        #Stores the abundance table
	--o-taxonomy cache:kaiju_taxonomy-trimmed-filtered \        #Stores the taxonomy file
	--use-cache cache \
	--verbose
```
 
The output from Kaiju is a BIOM formated feature table (.biom) and a taxonmy file.

Feature Table from Kaiju:

```{r KaijuTable, echo=FALSE}
kaijutable<-as.data.frame(as.matrix(biomformat::biom_data(read_biom("Kaiju Taxonomy Abundance/kaiju_feature-table-trimmed-filtered.biom"))))
head(kaijutable)

```
Taxonomy from Kaiju:

```{r KaijuTaxonomy, echo=FALSE}
kaijutaxonomy<-read.table("Kaiju Taxonomy Abundance/kaiju_taxonomy-trimmed-filtered.tsv",sep='\t',header = TRUE)
head(kaijutaxonomy)

```

## **7. Taxonomic Classification using Metaphlan**

MetaPhlAn 4 relies on unique clade-specific marker genes identified from ~1M microbial genomes for taxonomic assignment of metgenome sequence reads. Metaphlan 4.1.1 was used as a standalone installation as a miniconda environment, whose installation information can be gathered from here (https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4.1). Metaphlan automatically downloads its latest database "ChocoPhlan" when you use it for the first time, but you can also do it manually with this ```$ metaphlan --install --index <latest database index> --bowtie2db <database folder>```. The latest database indices can be found here (http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/). Once everything is installed, the following code was used to perform classification using Metaphlan for our datasets.

```{bash,eval=FALSE,echo=TRUE}

name="trimmed-filtered"        #dataset name

source activate qiime2-metagenome-2024.10

qiime tools export \        #Fetching the datset from .qza format to a folder format
	--input-path demux-paired/$name \
	--output-path demultiplexed/$name
	
conda activate mpa        #Activating the miniconda environment, here named as "mpa"	

mkdir -p MetaphlanAnalysis/BowTieOutput/$name MetaphlanAnalysis/Output/$name        #folders for storing Metaphlan output and the intermediated Bowtie alignment ouput

folder="demultiplexed/$name"        #folder storing our sequences

# Loop over sample names
for f in $(find ${folder} -maxdepth 1 -name "*R1_001.fastq.gz" | sed 's|.*/||' | sed 's/R1_001.fastq.gz//' | sort -u)
do 
    # Define output filename without _S[number]_L[number]_R
    clean_f=$(echo ${f} | sed -E 's/_S[0-9]+_L[0-9]+_//')

    # Run MetaPhlAn
     metaphlan ${folder}/${f}R1_001.fastq.gz,${folder}/${f}R2_001.fastq.gz --bowtie2db <database folder> --bowtie2out MetaphlanAnalysis/BowTieOutput/$name/${clean_f}.bowtie2.bz2 --input_type fastq --nproc 28 -o MetaphlanAnalysis/Output/$name/${clean_f}.output.txt	
     echo ${clean_f} is done
done

```

A typical Metaphlan output contains taxonomic relative abundances

```{r MetaphlanReport, echo=FALSE}
content<-readLines("MetaphlanOutput/trimmed-filtered/1-MDcol002-5-2.output.txt")
numbered_lines <- paste("[",seq_along(content), "]  ",content, sep = "")
writeLines(numbered_lines[1:20])
```
The file has a five line header (lines beginning with #). We've added [i]s (line numbers) to the header text above to help call out their roles more clearly.

  1. Indicates the version of the marker gene database that MetaPhlAn used in this run. There are currently ~1.1M unique clade-specific marker genes identified from ~100k reference genomes (~99,500 bacterial and archaeal and ~500 eukaryotic) included in the database.
  2. Provides a copy of the command that was run to produce this profile. This includes the path to the MetaPhlAn executable, the input file that was analyzed, and any custom parameters that were configured.
  3. Indicates the number of sample reads processed.
  4. Lists your sample name.
  5. Provides column headers from the profile data that follows.
More specifically, [5] includes four headers to be familiar with:

**clade_name**: The taxonomic lineage of the taxon reported on this line, ranging from Kingdom (e.g. Bacteria/Archaea) to species genome bin (SGB). Taxon names are prefixed to help indicate their rank: Kingdom: k__, Phylum: p__, Class: c__, Order: o__, Family: f__, Genus: g__, Species: s__, and SGB: t__.

**NCBI_tax_id**: The NCBI-equivalent taxon IDs of the named taxa from clade_name.

**relative_abundance**: The taxon's relative abundance in %. Since typical shotgun-sequencing-based taxonomic profile is relative (i.e. it does not provide absolute cell counts), clades are hierarchically summed. Each taxonomic level will sum to 100%. That is, the sum of all kingdom-level clades is 100%, the sum of all phylum-level clades is 100%, and so forth.

**additional_species**: Additional species names for cases where the metagenome profile contains clades that represent multiple species. The species listed in column 1 is the representative species in such cases.

## **8. Functional annotation using Humann 3.9**

HUMAnN is a pipeline for efficiently and accurately profiling the presence/absence and abundance of microbial pathways in a community from metagenomic or metatranscriptomic sequencing data (typically millions of short DNA/RNA reads). This process, referred to as functional profiling, aims to describe the metabolic potential of a microbial community and its members. You can install Humann 3.9 as a stand-alone through here (https://github.com/biobakery/humann). Once installed, it automatically downloads the ChocoPhlan database (if not already installed) to do taxonomic analysis using Metaphlan on its first run, but you will have to download the translated search database for anotation. You can do it using this code ```humann_databases --download uniref uniref90_diamond $INSTALL_LOCATION```, which will download the full UniRef90 database used in our analysis. After that, run this code to annotate the sequences using HUMAnN. (P.S: Unlike MetaPhlan, HUMAnN cannpot parse paired read sequences, so we will have to first concatenate both reads from each sample into a single .fastq file and then send it off for functional annotation)

```{bash,eval=FALSE,echo=TRUE}
name="trimmed-filtered"       #Dataset name

mkdir -p HumannAnalysis/HumannOutput/$name       #where the output is to be stored

folder="demultiplexed/$name"        #Dataset folder

cd $folder
forward=$(echo *R1*.fastq.gz)
reverse=$(echo *R2*.fastq.gz)

mkdir -p concatenated

for f in $forward
do
	for r in $reverse
	do
		clean_f=$(echo ${f} | sed -E 's/_S[0-9]+_L[0-9]+_R[0-9]_[0-9]+.fastq.gz//')
		clean_r=$(echo ${r} | sed -E 's/_S[0-9]+_L[0-9]+_R[0-9]_[0-9]+.fastq.gz//')
		
		if [ $clean_f = $clean_r ]; then
			cat $f $r > concatenated/$clean_f.fastq.gz        #concatenating the forward and reverse reads of the same sample
			humann --input concatenated/$clean_f.fastq.gz --threads 16 --output HumannAnalysis/HumannOutput/$name
			fi
	done
done
rm -rf concatenated       #deleting the folder storing the conactenated reads
```

HUMAnN output consists of three files: Gene families file (```$SAMPLENAME_genefamilies.tsv```), Path Coverage file (```$SAMPLENAME_pathcoverage.tsv```), and Path Abundance file (```$SAMPLENAME_pathabundance.tsv```).

The **gene families file** contains the the abundance of each gene family in the community. Gene families are groups of evolutionarily-related protein-coding sequences that often perform similar functions. Gene family abundance at the community level is stratified to show the contributions from known and unknown species. Individual species' abundance contributions sum to the community total abundance. Gene family abundance is reported in RPK (reads per kilobase) units to normalize for gene length; RPK units reflect relative gene (or transcript) copy number in the community.

```{r Humanngene, echo=FALSE}
humanngenfamily<-read.delim("HumannOutput/trimmed-filtered/1-MDcol002-5-2_genefamilies.tsv",sep='\t',header = TRUE,comment.char = "")
df<-colnames(humanngenfamily)%>%
    str_replace_all("X\\.\\.","")%>%
    str_replace_all("X","")%>%
    str_replace_all("\\.","-")
colnames(humanngenfamily)<-df
head(humanngenfamily)
```


The **path coverage file** provides information on the completeness of metabolic pathways in a given sample. It contains two main columns: Pathway, which lists metabolic pathway identifiers, and Coverage, which represents the proportion of the pathway detected in the sample as a value between 0 and 1. A coverage value of 0.00 indicates that the pathway is absent, 0.50 suggests that half of its reactions are present, and 1.00 means the pathway is fully reconstructed. These pathway identifiers are derived from the MetaCyc database.

```{r Humancoverage, echo=FALSE}
humanncoverage<-read.delim("HumannOutput/trimmed-filtered/1-MDcol002-5-2_pathcoverage.tsv",sep='\t',header = TRUE,comment.char = "")
df<-colnames(humanncoverage)%>%
    str_replace_all("X\\.\\.","")%>%
    str_replace_all("X","")%>%
    str_replace_all("\\.","-")
colnames(humanncoverage)<-df
head(humanncoverage)
```

The **path abundance file** details the abundance of each pathway in the community as a function of the abundances of the pathway's component reactions, with each reaction's abundance computed as the sum over abundances of genes catalyzing the reaction. Pathway abundance is proportional to the number of complete "copies" of the pathway in the community but unlike gene abundance, a pathway's community-level abundance is not necessarily the sum of its stratified abundance values. It consists of two main columns: Pathway, which contains metabolic pathway identifiers from MetaCyc, and Abundance, which represents the cumulative abundance of reactions contributing to the pathway

```{r Humanabundance, echo=FALSE}
humannabundance<-read.delim("HumannOutput/trimmed-filtered/1-MDcol002-5-2_pathabundance.tsv",sep='\t',header = TRUE,comment.char = "")
df<-colnames(humannabundance)%>%
    str_replace_all("X\\.\\.","")%>%
    str_replace_all("X","")%>%
    str_replace_all("\\.","-")
colnames(humannabundance)<-df
head(humannabundance)
```

After gathering all the data, we will move to performing furtehr data transfromation in R in steps 9 onwards.

### **Required packages for the following R code**

```{r packages, message=FALSE, warning=FALSE}

library(biomformat)
library(tidyverse)
library(data.table)
library(fs)
library(phyloseq)
library(writexl)
library(phylostratr)
library(ape)
library(taxize)
```

## **9. Transforming Kaiju Data**

The following code is to read the Kaiju feature table (.biom) and Kaiju taxonomy, merge them together to create an abundance table that has both the taxonomy and the reads in each sample for each taxa, provided all the biom files and taxonomy files are kept in a folder named "Kaiju Taxonomy Abundance".

```{r Kaiju data transformation}

Kaijubioms<-dir("Kaiju Taxonomy Abundance",pattern="*biom")    #Listing all biom files from the kaiju folder 
Kaijutaxa<-dir("Kaiju Taxonomy Abundance",pattern="kaiju_taxonomy*")  #Listing all taxonomy files from the kaiju folder
for (biom in Kaijubioms){
  for(taxa in Kaijutaxa){
    if(substring(biom,21,nchar(biom)-5)==substring(taxa,16,nchar(taxa)-4))
      name=substring(biom,21,nchar(biom)-5)
    OTUtable<-as.data.frame(as.matrix(biom_data(read_biom(paste("Kaiju Taxonomy Abundance",biom,sep="/")))))%>% #Reading biom files
      rownames_to_column("Feature.ID")
    taxatable<-read.table(paste("Kaiju Taxonomy Abundance",taxa,sep="/"),sep='\t',header = TRUE)      #Reading the taxonomy file
    Abundance<-inner_join(OTUtable,taxatable)%>%      #merging them together and separating the taxonomy string into separate columns for each taxonomic ranks
      select(-"Feature.ID")%>%
      relocate(Taxon)%>%
      mutate(Taxon=str_replace_all(Taxon,".__",""))%>%
      separate(Taxon,c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")%>%
      mutate(Phylum=if_else(is.na(Phylum),Domain,Phylum),
             Class=if_else(is.na(Class),Phylum,Class),
             Order=if_else(is.na(Order),Class,Order),
             Family=if_else(is.na(Family),Order,Family),
             Genus=if_else(is.na(Genus),Family,Genus),
             Species=if_else(is.na(Species),Genus,Species))
    write_xlsx(Abundance,paste("Kaiju Taxonomy Abundance/kaiju_abun-",name,".xlsx",sep=""))   #Writing the table as an excel file
  }
}
```

Sample Kaiju feature table
```{r, echo=FALSE}
head(OTUtable)
```

Sample Kaiju taxonomy file
```{r, echo=FALSE}
head(taxatable)
```

Sample Kaiju combined abundance table
```{r, echo=FALSE}
head(Abundance)
```

### **Creating Phyloseq Objects for the Kaiju data**

The following code packages the feature table and the taxonomy file into a single phyloseq object for the datasets ("raw", "trimmed-filtered" etc). A Phyloseq object usually contains,

  - a **sample metadata**, which contains the sample IDs and metadata pertaining to it.
  - a **feature table**, which is a matrix where the rows are the taxa/feature and columns are the sample IDs with read count for each taxa/feature.
  - a **taxonomy file**, which is a matrix where the rows are the taxa/feature and columns are the taxonomy path for each taxa/feature.
  - a **phylogenetic tree**, a  phylogenetic tree with the taxa/feature as the tip labels.

The phyloseq object can be made with one or more of these data file.

```{r Creating Kaiju phyloseq object}

Kaijubioms<-dir("Kaiju Taxonomy Abundance",pattern="*biom")   #Listing all biom files from the kaiju folder
Kaijutaxa<-dir("Kaiju Taxonomy Abundance",pattern="kaiju_taxonomy*")  #Listing all taxonomy files from the kaiju folder
dir.create("PhyloseqOutputs/KaijuPhyloseq")
for (biom in Kaijubioms){
  for(taxa in Kaijutaxa){
    if(substring(biom,21,nchar(biom)-5)==substring(taxa,16,nchar(taxa)-4))
      name=substring(biom,21,nchar(biom)-5)
    OTUtable<-as.matrix(biom_data(read_biom(paste("Kaiju Taxonomy Abundance",biom,sep="/")))) #Reading biom files
    taxatable<-read.table(paste("Kaiju Taxonomy Abundance",taxa,sep="/"),sep='\t',header = TRUE)%>%      #Reading the taxonomy file
      mutate(Taxon=str_replace_all(Taxon,".__",""))%>%
      separate(Taxon,c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")%>%
      mutate(Phylum=if_else(is.na(Phylum),Domain,Phylum),
             Class=if_else(is.na(Class),Phylum,Class),
             Order=if_else(is.na(Order),Class,Order),
             Family=if_else(is.na(Family),Order,Family),
             Genus=if_else(is.na(Genus),Family,Genus),
             Species=if_else(is.na(Species),Genus,Species))%>%
      column_to_rownames("Feature.ID")%>%
      as.matrix()
    physeq <- phyloseq(otu_table(OTUtable, taxa_are_rows = TRUE), tax_table(taxatable))   #Compiling the phyloseq object
    saveRDS(physeq,paste("PhyloseqOutputs/KaijuPhyloseq/phyloseq_kaiju-",name,".rds",sep="")) #Saving it as an .RDS file
  }
}
```

Summary of the sample phyloseq object
```{r, echo=FALSE}
summary(physeq)
physeq
```
## **10. Transforming Kraken2-Bracken Data**

The following code is to read the Bracken feature table (.biom) and Bracken taxonomy, merge them together to create an abundance table that has both the taxonomy and the reads in each sample for each taxa.

```{r Creating Bracken Taxonomic abundance Table}

brackentables<-dir(path="BrackenOutputs",pattern="table-bracken*")   #Listing all the folders containing the bracken feature table
brackentaxonomy<-dir(path="BrackenOutputs",pattern="taxonomy-bracken*")  #Listing all the folders containing the bracken taxonomy

dir.create("Kraken-Bracken Taxonomy Abundance")
for(each in brackentables){
  name<-substring(each,15)
  for(each2 in brackentaxonomy){
    name2<-substring(each2,18)
    if(name==name2)
    {table<-as.data.frame(as.matrix(biom_data(read_biom(paste("BrackenOutputs",each,"feature-table.biom",sep='/')))))%>% #Reading biom file
      rownames_to_column("Feature.ID")
    taxonomy<-read.table(paste("BrackenOutputs",each2,"taxonomy.tsv",sep='/'),sep="\t",header = TRUE)%>%  #Reading taxonomy file
      mutate(Feature.ID=as.character(Feature.ID))
    taxatable<-inner_join(table,taxonomy)%>%  #Joining into an abundance table
      select(-"Feature.ID")%>%
      relocate(Taxon)%>%
      mutate(Taxon=str_replace_all(Taxon,".__",""))%>%
      separate(Taxon,c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep=";")
    write_xlsx(taxatable,paste("Kraken-Bracken Taxonomy Abundance/Kraken-Bracken_abun-",name,".xlsx",sep="")) #Writing the abundance table as excel files
    }
  }
}
```
Sample Bracken feature table, taxonomy and combined abundance table

```{r, echo=FALSE}
head(table)
head(taxonomy)
head(taxatable)
```
### **Creating Phyloseq Objects for the Kraken2-Bracken data**

The following code packages the feature table, taxonomy file, and the phylogeentic tree (obtained from gracken, see Step 5) into a single phyloseq object for the datasets ("raw", "trimmed-filtered" etc). 

```{r Creating Bracken Phyloseq}

brackentables<-dir(path="BrackenOutputs",pattern="table-bracken*")   #Listing all the folders containing the bracken feature table
brackentaxonomy<-dir(path="BrackenOutputs",pattern="taxonomy-bracken*")  #Listing all the folders containing the bracken taxonomy
brackentree<-dir(path="BrackenOutputs",pattern="bracken-reports*")  #Listing all the folders containing the bracken reports and tree

dir.create("PhyloseqOutputs/BrackenPhyloseq")
counter=0
for(each in brackentables){
  name<-substring(each,15)
  table<-as.matrix(biom_data(read_biom(paste("BrackenOutputs",each,"feature-table.biom",sep='/')))) #Reading biom file
  taxonomy<-read.table(paste("BrackenOutputs/taxonomy-bracken-",name,"/taxonomy.tsv",sep=''),sep="\t",header = TRUE)%>%  #Reading taxonomy file
    mutate(Taxon=str_replace_all(Taxon,".__",""))%>%
    separate(Taxon,c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")%>%
    mutate(Phylum=if_else(is.na(Phylum),Domain,Phylum),
           Class=if_else(is.na(Class),Phylum,Class),
           Order=if_else(is.na(Order),Class,Order),
           Family=if_else(is.na(Family),Order,Family),
           Genus=if_else(is.na(Genus),Family,Genus),
           Species=if_else(is.na(Species),Genus,Species))%>%
    column_to_rownames("Feature.ID")%>%
    as.matrix()
  tree<-read_tree(paste("BrackenOutputs/bracken-reports-",name,"/output.tree",sep=''))   #Reading the phylogenetic tree for that dataset
  counter=counter+1
  physeqBracken<-phyloseq(otu_table(table,taxa_are_rows = TRUE),tax_table(taxonomy))  #creating the phyloseq object
  physeqBracken<-merge_phyloseq(physeqBracken,tree)   #Adding the phylogenetic tree to the phyloseq object
  saveRDS(physeqBracken,paste("PhyloseqOutputs/BrackenPhyloseq/phyloseq_bracken-",name,".rds",sep="")) #Saving it as an .RDS file
  
}
```

Summary of the phyloseq object generated from the Kraken2-Bracken data

```{r, echo=FALSE}
summary(physeqBracken)
physeqBracken
```

## **10. Transforming Metaphlan Data and creating phyloesq objects from it**

Since, MetaPhlan data does not have the regular feature table and taxonomy format as Kraken-2 or Kaiju, we will have to tarnsform the data to make it fit that format. Also, we will use this tree generated from all the species level features in ChocoPhlan found <a href="https://www.dropbox.com/scl/fi/9ld98f1rjc0vtzygb9atj/20201013_mpav30_speciesMod.nwk.zip?rlkey=wcuppse0uw9z4pv1j3k19tbwg&e=1&dl=0">here</a>.

```{r Creating Metaphlan Phyloseq}

MetaphlanOutputs<-dir("MetaphlanOutput")    #Folder where the MetaPhlan outputs are stored
dir.create("PhyloseqOutputs/MetaphlanPhyloseq")
MetaphlanTree<-read_tree("20201013_mpav30_speciesMod.nwk")    #The Tree made from ChocoPhlan database
new_names <- MetaphlanTree$tip.label
new_names<-new_names%>%
  str_replace_all("s__","")
MetaphlanTree$tip.label<-new_names
for(folder in MetaphlanOutputs){    #Compiling all the outputs from different sample into one feature abundance table and one taxonomy file
  Alltaxonomy<-NULL
  AllOTU<-NULL
  for(file in list.files(paste("MetaphlanOutput",folder,sep="/"))){
    content<-readLines(paste("MetaphlanOutput",folder,file,sep="/"))[-c(1:4)]   #Parsing the MetaPhlan files
    id<-substring(basename(file),1,nchar(basename(file))-11)    #Stores the Sample ID
    totalreads<-readLines(paste("MetaphlanOutput",folder,file,sep="/"))[3]    #Parsing the total reads for that sample
    totalreads<-as.numeric(gsub("^#(\\d+).*", "\\1", totalreads))
    table<-read.table(text=content,header=FALSE,sep="\t")%>%    #Transforming the Metaphlan file into feature table with taxonomy
      select(-4)%>%
      mutate(V1=str_replace_all(V1,".__",""),
             V1 = str_replace_all(V1, "\\|", ";"),
             reads=round(V3/100*totalreads)) %>%
      separate(V1,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep=";")%>%
      unique()%>%
      filter(!(is.na(Species)))
    taxonomy<-table[,c(1:7)]    #Storing the taxonomy
    OTU<-table[,c(7,10)]%>%     #Storing the feature table
      unique()%>%
      mutate(Sample=id)
    Alltaxonomy<-rbind(Alltaxonomy,taxonomy)    #Adding taxonomy extracted from that sample with the rest in that dataset
    AllOTU<-rbind(AllOTU,OTU)    #Adding feature table extracted from that sample with the rest in that dataset
}
  
  taxa_table<-Alltaxonomy%>%    #Prepping taxatable for phyloseq
    arrange("Species")%>%
    unique()
  rownames(taxa_table)<-taxa_table$Species
  otu_table<-AllOTU%>%      #Prepping OTU/feature table for phyloseq
    group_by(Sample,Species)%>%
    summarize(reads=sum(reads))%>%
    pivot_wider(id_cols = "Species",values_from = "reads",names_from = "Sample",values_fill = 0)%>%
    column_to_rownames("Species")
  physeqMetaphlan<-phyloseq(otu_table(as.matrix(otu_table),taxa_are_rows = TRUE),tax_table(as.matrix(taxa_table)))    #Creating the phyloseq object
  physeqMetaphlan<-merge_phyloseq(physeqMetaphlan,phy_tree(MetaphlanTree))    #Adding the tree to the phyloseq object
  saveRDS(physeqMetaphlan,paste("PhyloseqOutputs/MetaphlanPhyloseq/phyloseq_Metaphlan-",substring(folder,1,nchar(folder)),".rds",sep=""))    #Saving the phyloseq object as an .RDS file
}
```

Sample Metaphlan output for a sample from our datasets

```{r MetaphlanReport2, echo=FALSE}
content<-readLines("MetaphlanOutput/trimmed-filtered/1-MDcol002-5-2.output.txt")
numbered_lines <- paste("[",seq_along(content), "]  ",content, sep = "")
writeLines(numbered_lines[1:20])
```
Sample taxonomy and feature table extracted from the MetaPhlan outputs

```{r, echo=FALSE}
head(taxa_table)
head(select(as.data.frame(otu_table),c(1:10)))
```

Summary of the phyloseq object generated from the MetaPhlan data

```{r, echo=FALSE}
summary(physeqMetaphlan)
physeqMetaphlan
```

## **11. Creating Phyloseq Objects with the HUMAnN data**

The following code parses the genefemilies.tv, pathcoverage.tsv, and pathabundance.tv files and combines them from individual sample files into three phyloseq objects for each dataset.

```{r Creating Humann Phyloseq}

HumannOutput<-dir("HumannOutput")     #Folder containing all the HUMAnN outputs
dir.create("PhyloseqOutputs/HumannPhyloseq")
for (folder in HumannOutput){
  AllCoverage<-data.frame()
  AllAbundance<-data.frame()
  counter=1
  for(file in list.files(paste("HumannOutput",folder,sep="/"),pattern="*_pathabundance.tsv")){
    id<-substring(basename(file),1,nchar(basename(file))-18)    #Parsing the Sample ID
    abundance<-read.delim(paste("HumannOutput",folder,file,sep="/"),sep='\t',header = TRUE,comment.char = "")   #Reading the pathway abundance file for a sample
    df<-colnames(abundance)%>%
      str_replace_all("X\\.\\.","")%>%
      str_replace_all("X","")%>%
      str_replace_all("\\.","_")%>%
      str_replace_all("_Abundance","")
    colnames(abundance)<-df
    coverage<-read.delim(paste("HumannOutput/",folder,"/",id,"_pathcoverage.tsv",sep=""),sep='\t',header = TRUE,comment.char = "")   #Reading the pathway coverage file for a sample
    df<-colnames(coverage)%>%
  str_replace_all("X\\.\\.","")%>%
  str_replace_all("X","")%>%
  str_replace_all("\\.","-")%>%
  str_replace_all("_Coverage","")
    colnames(coverage)<-df
    
    #Combining the pathway abundances and pathway coverage data with those from the previous samples
    if (counter==1){
      AllAbundance<-abundance}
    else{AllAbundance<-left_join(AllAbundance,abundance)}
    if(counter==1){
    AllCoverage<-coverage}
    else{AllCoverage<-left_join(AllCoverage,coverage)}
    counter=counter+1
  }

  AllAbundance<-AllAbundance%>%
    column_to_rownames("Pathway")%>%
    as.matrix()
  AllCoverage<-AllCoverage%>%
    column_to_rownames("Pathway")%>%
    as.matrix()
  physeqcoverage<-phyloseq(otu_table(AllCoverage,taxa_are_rows=TRUE))     #Creating phyloseq object for the pathway abundance data
  physeqabundance<-phyloseq(otu_table(AllAbundance,taxa_are_rows=TRUE))     #Creating phyloseq object for the pathway coverage data
  saveRDS(physeqabundance,paste("PhyloseqOutputs/HumannPhyloseq/phyloseq_humann-",folder,"_pathwayabundance.rds",sep=""))
  saveRDS(physeqcoverage,paste("PhyloseqOutputs/HumannPhyloseq/phyloseq_humann-",folder,"_pathwaycoverage.rds",sep=""))
}

for (folder in HumannOutput){
  AllRPK<-data.frame()
  counter=1
  for(file in list.files(paste("HumannOutput",folder,sep="/"),pattern="*_genefamilies.tsv")){
    id<-substring(basename(file),1,nchar(basename(file))-17)
    #Reading the gene family abundance for a sample
    RPK<-read.delim(paste("HumannOutput",folder,file,sep="/"),sep='\t',header = TRUE,comment.char = "")
    df<-colnames(RPK)%>%
      str_replace_all("X\\.\\.","")%>%
      str_replace_all("X","")%>%
      str_replace_all("\\.","-")%>%
      str_replace_all("_Abundance-RPKs","")
    colnames(RPK)<-df
    #Joining the gene family abundance data of the current sample to the rest of the samples
    if (counter==1){
      AllRPK<-RPK}
    else{AllRPK<-left_join(AllRPK,RPK)}
    counter=counter+1
    }

  AllRPK<-AllRPK%>%
    column_to_rownames("Gene-Family")%>%
    as.matrix()
  physeqRPK<-phyloseq(otu_table(AllRPK,taxa_are_rows=TRUE))     #Creating phyloseq object for the gene family abundance data
  saveRDS(physeqRPK,paste("PhyloseqOutputs/HumannPhyloseq/phyloseq_humann-",folder,"_genefamilies.rds",sep=""))
}
```

Sample pathway abundance file
```{r, echo=FALSE}
head(abundance)
```
Combined pathway abundances of all samples in a dataset
```{r, echo=FALSE}
head(as.data.frame(AllAbundance))
```
Sample pathway coverage file
```{r, echo=FALSE}
head(coverage)
```
Combined pathway coverage of all samples in a dataset
```{r, echo=FALSE}
head(as.data.frame(AllCoverage))
```
Sample gene family abundance file
```{r, echo=FALSE}
head(RPK)
```
Combined gene family abundances of all samples in a dataset
```{r, echo=FALSE}
head(as.data.frame(AllRPK))
rm(AllRPK)
```

