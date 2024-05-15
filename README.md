# Identification & Assembly of Unknown Bacterial Genomes
### Damien Mann, Disha Hegde, Tim Ralich


## Background
* Data for project was pulled from tmp @ /tmp/gen711_project_data/genome-assembly-fqs
* 3 genomes were picked out from about 81
* sequenced w/ next-gen sequencing through Illumina
* Raw reads already on Ron
* Unknowns Chosen for Analysis:
  - Sample 12
  - Sample 2
  - Sample CC4

## Workflow
![image](https://github.com/dishahegde13/final_project_ddt/assets/158323878/4b349b18-0d2e-4203-8c34-aef2e8ad5a84)

Figure 1: Flowchart for bacterial genome assembly/identification

## Methods

### Raw reads
Initial data gathered from Ron
Untrimmed reads

### Trimmomatic 
Used on raw reads to cut out adapters and poorly matched reads
* Trimmed reads used for later analysis

### FASTQC
Used on trimmed and untrimmed reads for quality score.
* Used once before trim, and again after trim

### SPADes
Used trimmed reads to create workable assembled genome to use for assesment steps
* very long, dependent on server activity
* Utilized nohup to run in background

### PROKKA
Used to annotate completed genome from SPADes
* Also utilized nohup
* 16S sequence from output needed for BLAST

### QUAST
Used SPADes output to test for contig quality/length
* Gives data for potential genome coverage
* Also measures genes w/GC content

### BUSCO
Took SPADes output and tested for genome completeness
* Checking for single-copy orthologs that span 90%+ in most genomes

### BLAST
Two Different Tests
* 16S BLAST: Testing for sequences that match 16S seuqences
* nucleotide BLAST: Testing for occurance of contaminants (Non-bacterial DNA)
Used to give possible ID's to bacterial genome
* Multiple results with varying % DNA matches

### BWA-mem
Aligning SPADes-constructed de novo reference to original raw reads
* Creates coverage table to support ID predictions

### BlobTools
Combines three separate inputs to create output graph that confidently ID's bacterial genome
Inputs:
* contigs.fasta file from SPADes
* "hits" file from BLAST
* BAM file from BWA-mem

### Genome Filtering
Took BlobTools assessment and filtered through outliers to find reliable contigs to make conclusions on
* Similar contigs part of larger organism ID
* Required intuitive thinking as no direct coding pathway was given

## Results

### ID's Based on Results
* Sample 12: Streptomyces
  - Largest genome
  - Coverage full of "hypothetical proteins", no specific strain identified
* Sample 2: Streptomyces
  - Smaller version of sample 12
  - more contigs than 12, indicated worse overall coverage
* Sample CC4: Mycobacterium shigaense
  - Another small genome w/many contigs, not great coverage

### FASTQC Comparison: Mycobacterium shigaense

![image](https://github.com/dishahegde13/final_project_ddt/assets/158323878/8d78ddb6-ecd2-426a-a3da-3819f3b78c34)

Figure 2: Untrimmed forward read for Mycobacterium shigaense (CC4)

![image](https://github.com/dishahegde13/final_project_ddt/assets/158323878/1a5d0fbf-31a2-402b-a959-2ba3b924ffd6)

Figure 3: Trimmed forward read for Mycobacterium shigaense (CC4)

### QUAST Results

![image](https://github.com/dishahegde13/final_project_ddt/assets/158323878/32cdf872-2962-4f9d-9884-9acc922d8e82)

Figure 4: QUAST for Mycobacterium shigaense outlining contig amount/lengths

![image](https://github.com/dishahegde13/final_project_ddt/assets/158323878/2cb80c23-43d0-4af5-b255-c99126e08811)

Figure 5: QUAST for Streptomyces 12 outlining contig amount/lengths

### BUSCO Results

![image](https://github.com/dishahegde13/final_project_ddt/assets/158323878/2f61a586-bfe8-4183-a86c-c28a2fe9b2ef)

Figure 6: BUSCO report for Mycobacterium shigaense 

![image](https://github.com/dishahegde13/final_project_ddt/assets/158323878/f5f67ac5-5265-47c1-9498-1404318dfd96)

Figure 7: BUSCO report for Streptomyces 12

### Proksee Genetic Maps

![image](https://github.com/dishahegde13/final_project_ddt/assets/158323878/dee9546f-5760-4518-8e2f-3b1a756e0bae)

Figure 8: Genetic map of Mycobacterium shigaense

![image](https://github.com/dishahegde13/final_project_ddt/assets/158323878/f792c82f-735c-4431-b5eb-ebfadccbe521)

Figure 9: Genetic map of Streptomyces 12

## Script 
``` bash
#! /bin/bash

# Getting Raw reads
mkdir -p ~/final_project_exe/raw_reads
cp /tmp/gen711_project_data/genome-assembly-fqs/12* ~/final_project_exe/raw_reads
cp /tmp/gen711_project_data/genome-assembly-fqs/CC4* ~/final_project_exe/raw_reads
cp /tmp/gen711_project_data/genome-assembly-fqs/2_*  ~/final_project_exe/raw_reads

# Checking Read Quality
conda activate genomics
mkdir  ~/final_project/fastqc_raw_reads
fastqc ~/final_project/raw_reads/* -o ~/final_project/fastqc_raw_reads
## used filezilla to move .html files to desktop and open fastqc.html files on chrome

# Genome 2
# Trimming reads using trimmomatic 
cd ~/final_project/raw_reads
trim_scriptV2.sh 2_S26_L001_R1_001.fastq.gz 2_S26_L001_R2_001.fastq.gz
cd trimmed-reads

# Checking Read Quality of Trimmed reads
fastqc 2* -o ~/final_project/fastqc_trimmed_raw_reads
## used filezilla to move .html files to desktop and open fastqc.html files 

# Genome assembly
spades.py -1 ~/final_project/raw_reads/trimmed-reads/2_S26_L001_R1_001.fastq.gz -2 ~/final_project/raw_reads/trimmed-reads/2_S26_L001_R2_001.fastq.gz -s ~/final_project/raw_reads/trimmed-reads/unpaired-2_S26_L001_R1_001.fastq.gz -s ~/final_project/raw_reads/trimmed-reads/unpaired-2_S26_L001_R2_001.fastq.gz -o ~/final_project/genome_assembly_2 -t 24 
## used nohup to run in the background
cd ~/final_project/genome_assembly_2/
rm -r tmp K* corrected 

# Quast for genome structure assessment 
quast.py ~/final_project/genome_assembly_2/contigs.fasta -o  ~/final_project/genome_assembly_2/quast_results
## used filezilla to move .html files to desktop and open .html files

# BUSCO for genome completeness assessment 
conda activate busco
busco -i ~/final_project/genome_assembly_2/contigs.fasta -m genome -o busco-results -l bacteria
conda activate genomics

# Prokka for genome annotation 
prokka ~/final_project/genome_assembly_2/contigs.fasta --outdir prokka_output --cpus 24 --mincontiglen 200 
grep -o "product=.*" prokka_output/PROKKA_*.gff | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_abundances.txt
extract_sequences "16S ribosomal RNA" prokka_output/PROKKA_*.ffn > 16S_sequence.fasta

# BLAST
makeblastdb -in ~/final_project/genome_assembly_2/contigs.fasta -dbtype nucl -out contigs_db 
blastn -query 16S_sequence.fasta -db contigs_db -out 16S_vs_contigs_6.tsv -outfmt 6
blob_blast.sh ~/final_project/genome_assembly_2/contigs.fasta

# Read Mapping
bwa index ~/final_project/genome_assembly_2/contigs.fasta
bwa mem -t 24 ~/final_project/genome_assembly_2/contigs.fasta ~/final_project/raw_reads/trimmed-reads/2_S26_L001_R1_001.fastq.gz ~/final_project/raw_reads/trimmed-reads/2_S26_L001_R2_001.fastq.gz > raw_mapped.sam
samtools view -@ 24 -Sb  raw_mapped.sam  | samtools sort -@ 24 - sorted_mapped
samtools index sorted_mapped.bam
bedtools genomecov -ibam sorted_mapped.bam > coverage.out
gen_input_table.py  --isbedfiles ~/final_project/genome_assembly_2/contigs.fasta coverage.out > coverage_table.tsv

# Non-target contig removal
blobtools create -i ~/final_project/genome_assembly_2/contigs.fasta -b sorted_mapped.bam -t contigs.fasta.vs.nt.cul5.1e5.megablast.out -o blob_out
blobtools view -i blob_out.blobDB.json -r all -o blob_taxonomy
blobtools plot -i blob_out.blobDB.json -r genus

# Filtering the genome 
grep -v '#' blob_taxonomy.blob_out.blobDB.table.txt | awk -F'\t' '$5 > 50' | awk -F'\t' '$3 > 0.6' | awk -F'\t' '{print $1}' > list_of_contigs_to_keep_gc60%_cov50.txt
filter_contigs_by_list.py ~/final_project/genome_assembly_2/contigs.fasta list_of_contigs_to_keep_gc60%_cov50.txt bacteria2_filtered.fasta
## Average coverage from Blobtools table
grep -f list_to_keep.txt blob_taxonomy.blob_out.blobDB.table.txt | awk '{w = w + $2; e = e + $5 * $2;} END {print e/w}'

# BLAST the final contigs against UniVec
wget "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec"
blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -query bacteria2_filtered.fasta -subject UniVec  -outfmt 6 -out genome_vs_univec.6

# Genome 12
# Trimming reads using trimmomatic 
cd ~/final_project/raw_reads
trim_scriptV2.sh 12_S36_L001_R1_001.fastq.gz 12_S36_L001_R2_001.fastq.gz
cd trimmed-reads

# Checking Read Quality of Trimmed reads
fastqc 12* -o ~/final_project/fastqc_trimmed_raw_reads
## used filezilla to move .html files to desktop and open fastqc.html files 

# Genome assembly
spades.py -1 ~/final_project/raw_reads/trimmed-reads/12_S36_L001_R1_001.fastq.gz -2 ~/final_project/raw_reads/trimmed-reads/12_S36_L001_R2_001.fastq.gz -s ~/final_project/raw_reads/trimmed-reads/unpaired-12_S36_L001_R1_001.fastq.gz -s ~/final_project/raw_reads/trimmed-reads/unpaired-12_S36_L001_R2_001.fastq.gz -o ~/final_project/genome_assembly_12 -t 24 
## used nohup to run in the background
cd ~/final_project/genome_assembly_12/
rm -r tmp K* corrected 

# Quast for genome structure assessment 
quast.py ~/final_project/genome_assembly_12/contigs.fasta -o  ~/final_project/genome_assembly_12/quast_results
## used filezilla to move .html files to desktop and open .html files

# BUSCO for genome completeness assessment 
conda activate busco
busco -i ~/final_project/genome_assembly_12/contigs.fasta -m genome -o busco-results -l bacteria
conda activate genomics

# Prokka for genome annotation 
prokka ~/final_project/genome_assembly_12/contigs.fasta --outdir prokka_output --cpus 24 --mincontiglen 200 
grep -o "product=.*" prokka_output/PROKKA_*.gff | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_abundances.txt
extract_sequences "16S ribosomal RNA" prokka_output/PROKKA_*.ffn > 16S_sequence.fasta

# BLAST
makeblastdb -in ~/final_project/genome_assembly_12/contigs.fasta -dbtype nucl -out contigs_db 
blastn -query 16S_sequence.fasta -db contigs_db -out 16S_vs_contigs_6.tsv -outfmt 6
blob_blast.sh ~/final_project/genome_assembly_12/contigs.fasta

# Read Mapping
bwa index ~/final_project/genome_assembly_12/contigs.fasta
bwa mem -t 24 ~/final_project/genome_assembly_12/contigs.fasta ~/final_project/raw_reads/trimmed-reads/12_S36_L001_R1_001.fastq.gz ~/final_project/raw_reads/trimmed-reads/12_S36_L001_R2_001.fastq.gz > raw_mapped.sam
samtools view -@ 24 -Sb  raw_mapped.sam  | samtools sort -@ 24 - sorted_mapped
samtools index sorted_mapped.bam
bedtools genomecov -ibam sorted_mapped.bam > coverage.out
gen_input_table.py  --isbedfiles ~/final_project/genome_assembly_12/contigs.fasta coverage.out > coverage_table.tsv

# Non-target contig removal
blobtools create -i ~/final_project/genome_assembly_12/contigs.fasta -b sorted_mapped.bam -t contigs.fasta.vs.nt.cul5.1e5.megablast.out -o blob_out
blobtools view -i blob_out.blobDB.json -r all -o blob_taxonomy
blobtools plot -i blob_out.blobDB.json -r genus

# Filtering the genome 
grep -v '#' blob_taxonomy.blob_out.blobDB.table.txt | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 20' | awk -F'\t' '{print $1}' > list_of_contigs_to_keep_gc60%_cov50.txt
filter_contigs_by_list.py ~/final_project/genome_assembly_12/contigs.fasta list_of_contigs_to_keep_len500_cov20.txt bacteria12_filtered.fasta
## Average coverage from Blobtools table
grep -f list_to_keep.txt blob_taxonomy.blob_out.blobDB.table.txt | awk '{w = w + $2; e = e + $5 * $2;} END {print e/w}'

# BLAST the final contigs against UniVec
wget "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec"
blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -query bacteria12_filtered.fasta -subject UniVec  -outfmt 6 -out genome_vs_univec.6

# Genome CC4
# Trimming reads using trimmomatic 
cd ~/final_project/raw_reads
trim_scriptV2.sh CC4-5-1_S40_L001_R1_001.fastq.gz CC4-5-1_S40_L001_R2_001.fastq.gz
cd trimmed-reads

# Checking Read Quality of Trimmed reads
fastqc 2* -o ~/final_project/fastqc_trimmed_raw_reads
## used filezilla to move .html files to desktop and open fastqc.html files 

# Genome assembly
spades.py -1 ~/final_project/raw_reads/trimmed-reads/CC4-5-1_S40_L001_R1_001.fastq.gz -2 ~/final_project/raw_reads/trimmed-reads/CC4-5-1_S40_L001_R2_001.fastq.gz -s ~/final_project/raw_reads/trimmed-reads/unpaired-CC4-5-1_S40_L001_R1_001.fastq.gz -s ~/final_project/raw_reads/trimmed-reads/unpaired-CC4-5-1_S40_L001_R2_001.fastq.gz -o ~/final_project/genome_assembly_CC4 -t 24 
## used nohup to run in the background
cd ~/final_project/genome_assembly_CC4/
rm -r tmp K* corrected 

# Quast for genome structure assessment 
quast.py ~/final_project/genome_assembly_CC4/contigs.fasta -o  ~/final_project/genome_assembly_CC4/quast_results
## used filezilla to move .html files to desktop and open .html files

# BUSCO for genome completeness assessment 
conda activate busco
busco -i ~/final_project/genome_assembly_CC4/contigs.fasta -m genome -o busco-results -l bacteria
conda activate genomics

# Prokka for genome annotation 
prokka ~/final_project/genome_assembly_CC4/contigs.fasta --outdir prokka_output --cpus 24 --mincontiglen 200 
grep -o "product=.*" prokka_output/PROKKA_*.gff | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_abundances.txt
extract_sequences "16S ribosomal RNA" prokka_output/PROKKA_*.ffn > 16S_sequence.fasta

# BLAST
makeblastdb -in ~/final_project/genome_assembly_CC4/contigs.fasta -dbtype nucl -out contigs_db 
blastn -query 16S_sequence.fasta -db contigs_db -out 16S_vs_contigs_6.tsv -outfmt 6
blob_blast.sh ~/final_project/genome_assembly_CC4/contigs.fasta

# Read Mapping
bwa index ~/final_project/genome_assembly_CC4/contigs.fasta
bwa mem -t 24 ~/final_project/genome_assembly_2/contigs.fasta ~/final_project/raw_reads/trimmed-reads/2_S26_L001_R1_001.fastq.gz ~/final_project/raw_reads/trimmed-reads/2_S26_L001_R2_001.fastq.gz > raw_mapped.sam
samtools view -@ 24 -Sb  raw_mapped.sam  | samtools sort -@ 24 - sorted_mapped
samtools index sorted_mapped.bam
bedtools genomecov -ibam sorted_mapped.bam > coverage.out
gen_input_table.py  --isbedfiles ~/final_project/genome_assembly_CC4/contigs.fasta coverage.out > coverage_table.tsv

# Non-target contig removal
blobtools create -i ~/final_project/genome_assembly_CC4/contigs.fasta -b sorted_mapped.bam -t contigs.fasta.vs.nt.cul5.1e5.megablast.out -o blob_out
blobtools view -i blob_out.blobDB.json -r all -o blob_taxonomy
blobtools plot -i blob_out.blobDB.json -r genus

# Filtering the genome 
grep -v '#' blob_taxonomy.blob_out.blobDB.table.txt | awk -F'\t' '$5 > 50' | awk -F'\t' '$2 > 1000' | awk -F'\t' '{print $1}' > list_of_contigs_to_keep_len1000_cov50.txt
filter_contigs_by_list.py ~/final_project/genome_assembly_CC4/contigs.fasta list_of_contigs_to_keep_len1000_cov50.txt bacteriaCC4_filtered.fasta
## Average coverage from Blobtools table
grep -f list_to_keep.txt blob_taxonomy.blob_out.blobDB.table.txt | awk '{w = w + $2; e = e + $5 * $2;} END {print e/w}'

# BLAST the final contigs against UniVec
wget "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec"
blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -query bacteriaCC4_filtered.fasta -subject UniVec  -outfmt 6 -out genome_vs_univec.6
```

## Conclusions
* Guide very helpful for completing project
  - Straightforward and concise
* Proksee great additional tool for helping with visualization
* Some graphs harder to understand than others
  - Certain analysis metrics hard to understand without prior understanding of program

### Helpful Links
* Proksee site: https://proksee.ca/
* Whole-Genome Assembly and Assessment Tutorial: https://github.com/Joseph7e/MDIBL-T3-WGS-Tutorial
