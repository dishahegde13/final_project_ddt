# Identification & Assembly of Unknown Bacterial Genomes
### Damien Mann, Disha Hegde, Tim Ralich


## Background
* Data for project was pulled from tmp @ /tmp/gen711_project_data/genome-assembly-fqs
* 3 genomes were picked out from about 81
* sequenced w/ next-gen sequencing through Illumina
* Raw reads already on Ron

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
