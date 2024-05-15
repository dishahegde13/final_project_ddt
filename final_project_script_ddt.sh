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