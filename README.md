TFBS_footprinting
=================
Pipeline: Identification of cis-regulatory elements by matrix scoring and conservation.

Purpose:
The goal of this script is to provide a pipeline for the transcription factor binding site (TFBS) footprinting method.

Input:
- Dictionary of: {transcript_id : Ensembl_transcript_id}
- Gene symbol (e.g. 'CA1')
- Dinucleotide TFBS position frequency matrix dictionary as JSON file (PFM) (output from x.py)
- Mononucleotide TFBS PFM dictionary as JSON file (output from x.py)
- Dictionary of unique TFBSs for each TF as JSON file (output from x.py)
- Ensembl species group designation (e.g. mammals)
- Phylogenetic footprinting predictions from Bigfoot program (optional, run separately)
- Statistical alignment of orthologous sequences from Bigfoot program (optional, run separately)
- Various score thresholds

Process:
Iterate through each transcript in the dictionary and:
- Retrieve EPO (Enredo-Pecan-Ortheus) aligned orthologous sequences from Ensembl database (determined by species group designation) corresponding to 1000 nt upstream of transcript
- Generate dinucleotide weight matrices from Jaspar TFBSs for each ortholog sequence
- Score all experimentally determined TFBSs, and set lowest score as cut-off for denovo TFBS prediction
- Score each sequence using all weight matrices while identifying matches to sites determined experimentally previously
- Identify predicted/matched human TFBSs conserved across multiple species
- For each conserved (human and 1+ additional species) TFBS, generate 'combined affinity score' as a sum of position weight scores across species

Output:

-hello
