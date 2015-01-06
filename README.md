TFBS_footprinting
=================
Pipeline: Identification of cis-regulatory elements by matrix scoring and conservation.

## Purpose
The goal of these scripts is to provide a pipeline for the transcription factor binding site (TFBS) footprinting method.

## Input
- Dictionary of: {transcript_id : Ensembl_transcript_id}
- Di-nucleotide TFBS position frequency matrix dictionary as JSON file (PFM) (output from 'sites2dfm.py')
- Mono-nucleotide TFBS PFM dictionary as JSON file (output from 'convert_jaspar_matrix.py')
- Dictionary of unique TFBSs for each TF as JSON file (output from 'sites2dfm.py')
- Ensembl species group designation (e.g. mammals)
- Phylogenetic footprinting predictions from Bigfoot program (optional, run separately)
- Statistical alignment of orthologous sequences from Bigfoot program (optional, run separately)
- Ist.medisapiens.com tissue coexpression tables (optional)
- Various score thresholds


## Process
Iterate through each transcript in the dictionary and:
 1. Retrieve EPO (Enredo-Pecan-Ortheus) aligned orthologous sequences from Ensembl database (determined by species group designation) corresponding to 1000 nt upstream of transcript
 2. Edit retrieved alignment:
- Remove sequences that are less than 75% length of human sequence
- Replace characters not corresponding to nucleotides (ACGT), with gaps
- Remove gap-only columns from alignment
- Starting from the beginning of alignment, remove all columns until the first position where at least 50% of sequences are non-gap
 3. Generate di-nucleotide weight matrices from Jaspar TFBSs for each ortholog sequence
 4. Score all experimentally determined TFBSs, and set lowest score as cut-off for denovo TFBS prediction
 5. Score each sequence using all weight matrices while identifying matches to sites determined experimentally previously
 6. Identify predicted/matched human TFBSs conserved across multiple species
 7. For each conserved (human and 1+ additional species) TFBS, generate 'combined affinity score' as a sum of position weight scores across species

## Output
- Figure showing 10 highest scoring (combined affinity score) TFBSs mapped onto human promoter
- CSV table file for each transcript with all predicted binding sites for all species
- CSV table file for each transcript with conserved binding sites (human and 1+ additional species)
- FASTA file of 1,000 nucleotide (nt) EPO alignment
- JSON file of dictionary of predicted TFBSs

## Dependencies
- Currently only tested on Linux
- Python 2.7
- Biopython (http://biopython.org/wiki/Download#Installation_Instructions)
- clustalo `sudo apt-get install clustalo`
- matplotlib `sudo apt-get install python-matplotlib`

## Usage
 1. Download TFBS_footprinting git project files to same directory
 2. Download 'sites' folder from 'http://jaspar.genereg.net/html/DOWNLOAD/' to current directory
 3. Execute sites2dfm.py to create di-nucleotide frequency matrices and corresponding .json output files
 4. Execute convert_jaspar_matrix.py to create mono-nucleotide frequency matrices and corresponding .json output file from 'pfm_vertebrates.txt'
 5. Edit 'Variables/Thresholds/Input' section of TFBS_analyzer_dinuc.py' and execute
