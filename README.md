# Identification & Assembly of Unknown Bacterial Genomes
### Damien Mann, Disha Hegde, Tim Ralich


## Background
* Data for project was pulled from tmp @ /tmp/gen711_project_data/genome-assembly-fqs
* 3 genomes were picked out from about 81
* sequenced w/ next-gen sequencing through Illumina
* Raw reads already on Ron
Unknowns Chosen for Analysis:
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

## Conclusions
* Guide very helpful for completing project
  - Straightforward and concise
* Proksee great additional tool for helping with visualization
* Some graphs harder to understand than others
  - Certain analysis metrics hard to understand without prior understanding of program

### Helpful Links
Proksee site: https://proksee.ca/
Whole-Genome Assembly and Assessment Tutorial: https://github.com/Joseph7e/MDIBL-T3-WGS-Tutorial
