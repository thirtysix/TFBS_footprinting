TFBS_footprinting
=================
Pipeline: Identification of cis-regulatory elements by matrix scoring and conservation across groups of species (mammals, primates, sauropsids, fish) catalogued in the Ensembl database.

## Purpose
Pipeline for the transcription factor binding site (TFBS) footprinting method.  Predict TFBSs in a target species (e.g. homo sapiens), and identify if it is conserved in a group of related species (e.g. mammals or primates) using alignments extracted from the Ensembl database.  Uses the Jaspar 2016 vertebrate binding site data of 519 TFBSs.

## Dependencies
- `git clone https://github.com/thirtysix/TFBS_footprinting.git`
- Python 2.7
- Numpy `sudo apt-get install python-numpy`
- Biopython `sudo apt-get install python-biopython`
- matplotlib `sudo apt-get install python-matplotlib`
- Currently only tested on Linux

## User Input
<code>
Example Usage:
simplest:
TFBS_analyzer2.py PATH_TO/sample_ids.txt

all arguments:
TFBS_analyzer2.py PATH_TO/sample_ids.txt -s homo_sapiens -g mammals -pb 900 -pa 100 -l 5 -c 2 -tx 10 -o PATH_TO/Results/

positional arguments:
Required: Location of a file containing Ensembl mammal
transcript ids (see sample file: sample_ids.txt)")

- target_species , -s 
[default: "homo_sapiens"] - Target species (string),options are located at (https://rest.ensembl.org/info/compara/species_sets/EPO_LOW_COVERAGE?content-type=application/json). Conservation of TFs acrossother species will be based on identifying them inthis species first.
- species_group , -g 
[default: "mammals"] - Group of species (string) toidentify conservation of TFs within. Your targetspecies should be a member of this species group (e.g."homo_sapiens" and "mammals" or "primates". Options:"mammals", "primates", "sauropsids", "fish". Groupsand members are listed at (https://rest.ensembl.org/info/compara/species_sets/EPO_LOW_COVERAGE?content-type=application/json)
- promoter_before_tss , -pb 
[default: 900] - Number (integer) of nucleotidesupstream of TSS to include in analysis.
- promoter_after_tss , -pa 
[default: 100] - Number (integer) of nucleotidesdownstream of TSS to include in analysis.
- locality_threshold , -l 
[default: 5] - Nucleotide distance (integer)upstream/downstream in which TF predictions in otherspecies will be included to support a hit in the target species.
- conservation_min , -c 
[default: 2] - Minimum number (integer) of species apredicted TF is found in, in alignment, to beconsidered conserved.
- top_x_tfs , -tx [default: 10] - Number (integer) of unique TFs toinclude in output .svg figure.
- output_dir , -o [default: /home/harlan/Dropbox/manuscripts/tfbs_footprinting/8.somewhere/scripts/testing/Results ] - Fullpath of directory where result directories will beoutput.
</code>

## Process
Iterate through each user provided Ensembl transcript id:
 1. Retrieve EPO aligned orthologous sequences from Ensembl database for user-chosen species group (mammals, primates, fish, sauropsids).
 2. Edit retrieved alignment:
- Remove sequences that are less than 75% length of target_species sequence
- Replace characters not corresponding to nucleotides (ACGT), with gaps
- Remove gap-only columns from alignment
 3. Generate position weight matrices (PWMs) from Jaspar position frequency matrices (PFMs).
 4. Score each species sequence in the user-chosen species group using all weight matrices.
 5. Keep predictions with a score greater than score threshold corresponding to p-value of 0.001.
 6. Identify predicted TFBSs in target_species which are conserved in non-target_species species of the the species_group within the locality_threshold and totaling at least conservation_min.
 7. For each conserved TFBS, compute 'combined affinity score' as a sum of position weight scores of species possessing a prediction.
 8. Sort target_species predictions by combined affinity score, generate a vector graphics figure showing top_x_tfs unique TFs mapped onto the promoter of the target transcript.


## Output
- Original alignment as retrieved from Ensembl (alignment_uncleaned.fasta).
- Cleaned alignment (alignment_cleaned.fasta).
- Regulatory information for the target transcripts user-defined promoter region (regulatory_decoded.json).
- Transcript properties for target transcript (transcript_dict.json).
- All predicted TFBSs for all species which satisfy p-value threshold (TFBSs_found.all.json).
- All predicted TFBSs for target species which are supported by at least conservation_min predictions in other species, and those supporting species, grouped into clusters (TFBSs_found.clusters.csv).
- All predicted TFBSs for target species which are supported by at least conservation_min predictions in other species, sorted by combined affinity score (TFBSs_found.sortedclusters.csv).
- Figure showing top_x_tfs highest scoring (combined affinity score) TFBSs mapped onto target_species promoter (ENST00000285379_mammals.Promoterhisto.svg). 
