##TFBS_footprinting
=================
Pipeline: Identification of cis-regulatory elements by matrix scoring and conservation across groups of species (mammals, primates, sauropsids, fish) catalogued in the Ensembl database.

## 1 Background
The TFBS footprinting method computationally predicts transcription factor binding sites (TFBSs) in a target species (e.g. homo sapiens), and identifies if they are conserved in a group of related species (e.g. mammals or primates) using alignments extracted from, and established within, the Ensembl database.  Our analyses have shown this method is superior in prediction of functional TFBSs compared to methods (PWMs and DeepBind) analyzing single species alone.  At the same time, TFBS footprinting is an extension of these methods.  This implementation uses the Jaspar 2016 vertebrate binding site data of 519 TFBSs.  This tool allows the user to analyze the user-defined promoter regions near the transcription start site of many vertebrate transcripts defined in the Ensembl database.

## 2 Usage 
Predict TFBSs in the promoter of from 1-80,000 human protein coding transcripts, and identify if they are conserved in mammals or primates.
Model organisms such as mouse and zebrafish.
39 Mammals
ee the Ensembl species groups to plan your analysis: https://rest.ensembl.org/info/compara/species_sets/EPO_LOW_COVERAGE?content-type=application/json

### 2.1 Installation
- `pip install TFBS_footprinting`

### 2.2 Dependencies
- Currently only tested on Linux
- Python 2.7
- Installed by *pip* automatically: Http2, Numpy, Biopython, Matplotlib


### 2.3 User Input Examples - File of Ensembl Transcript IDs
![sample_ids](https://raw.githubusercontent.com/thirtysix/TFBS_footprinting/master/sample_ids.txt.png)


### 2.4 User Input Examples - Command Line
```
$ TFBS_footprinter PATH_TO/sample_ids.txt
$ TFBS_footprinter PATH_TO/sample_ids.txt -s homo_sapiens -g mammals -pb 900 -pa 100 -l 5 -c 2 -tx 10 -o PATH_TO/Results/
```

### 2.5 Arguments
```
usage: TFBS_analyzer2.py [-h] [--tf_ids_file] [--target_species]
                         [--species_group] [--coverage]
                         [--promoter_before_tss] [--promoter_after_tss]
                         [--locality_threshold] [--conservation_min]
                         [--top_x_tfs] [--output_dir]

TFBS Footprinting - Identification of conserved vertebrate transcription factor binding sites (TFBSs)')

------------------------------------------------------------------------------------------------------
Example Usage:
    simplest:
    TFBS_analyzer2.py PATH_TO/sample_ids.txt

    all arguments:
    TFBS_analyzer2.py PATH_TO/sample_ids.txt -s homo_sapiens -g mammals -e low -pb 900 -pa 100 -l 5 -c 2 -tx 10 -o PATH_TO/Results/
------------------------------------------------------------------------------------------------------

positional arguments:
                        Required: Location of a file containing Ensembl
                        target_species transcript ids (see sample file
                        sample_ids.txt at
                        https://github.com/thirtysix/TFBS_footprinting)")

optional arguments:
  -h, --help            show this help message and exit
  --tf_ids_file , -tfs 
                        Optional: Location of a file containing a limited list
                        of Jaspar TFs to use in scoring alignment (see sample
                        file tf_ids.txt at
                        https://github.com/thirtysix/TFBS_footprinting)
                        [default: all Jaspar TFs]
  --target_species , -s 
                        [default: "homo_sapiens"] - Target species (string),
                        options are located at (https://github.com/thirtysix/T
                        FBS_footprinting/blob/master/README.md#species).
                        Conservation of TFs across other species will be based
                        on identifying them in this species first.
  --species_group , -g 
                        ("mammals", "primates", "sauropsids", or "fish")
                        [default: "mammals"] - Group of species (string) to
                        identify conservation of TFs within. Your target
                        species should be a member of this species group (e.g.
                        "homo_sapiens" and "mammals" or "primates"). The
                        "primates" group does not have a low-coverage version.
                        Groups and members are listed at (https://github.com/t
                        hirtysix/TFBS_footprinting/blob/master/README.md#speci
                        es)
  --coverage , -e       ("low" or "high") [default: "low"] - Which Ensembl EPO
                        alignment of species to use. The low coverage contains
                        significantly more species and is recommended. The
                        primate group does not have a low-coverage version.
  --promoter_before_tss , -pb 
                        (0-100,000) [default: 900] - Number (integer) of
                        nucleotides upstream of TSS to include in analysis
                        (0-100,000).
  --promoter_after_tss , -pa 
                        (0-100,000) [default: 100] - Number (integer) of
                        nucleotides downstream of TSS to include in analysis.
  --locality_threshold , -l 
                        (0-100) [default: 5] - Nucleotide distance (integer)
                        upstream/downstream within which TF predictions in
                        other species will be included to support a hit in the
                        target species.
  --conservation_min , -c 
                        (1-20)[default: 2] - Minimum number (integer) of
                        species a predicted TF is found in, in alignment, to
                        be considered conserved .
  --top_x_tfs , -tx     (1-20) [default: 10] - Number (integer) of unique TFs
                        to include in output .svg figure.
  --output_dir , -o     [default: /home/harlan/Dropbox/github/TFBS_footprintin
                        g/tfbs_results ] - Full path of directory where result
                        directories will be output.

```


## 3 Process
Iterate through each user provided Ensembl transcript id:
 1. Retrieve EPO aligned orthologous sequences from Ensembl database for user-defined species group (mammals, primates, fish, sauropsids) at for promoter of transcript id.
 2. Edit retrieved alignment:
- Remove species sequences that are less than 75% length of target_species sequence.
- Replace characters not corresponding to nucleotides (ACGT), with gaps characters "-".
- Remove gap-only columns from alignment.
 3. Generate position weight matrices (PWMs) from 519 TF Jaspar position frequency matrices (PFMs) for each species sequence based on nt background frequencies of target sequence.
 4. Score each species sequence in the alignment using all PWMs.
 5. Keep predictions with a score greater than score threshold corresponding to p-value of 0.001.
 6. Identify predicted TFBSs in target_species which are conserved in non-target_species species of the the species_group within the locality_threshold and totaling at least conservation_min.
 7. For each conserved TFBS, compute 'combined affinity score' as a sum of position weight scores of species possessing a prediction.
 8. Sort target_species predictions by combined affinity score, generate a vector graphics figure showing top_x_tfs unique TFs mapped onto the promoter of the target transcript, and additional output as described below.


## 4 Output
- Original alignment as retrieved from Ensembl (alignment_uncleaned.fasta).
- Cleaned alignment (alignment_cleaned.fasta).
- Regulatory information for the target transcripts user-defined promoter region (regulatory_decoded.json).
- Transcript properties for target transcript (transcript_dict.json).
- All predicted TFBSs for all species which satisfy p-value threshold (TFBSs_found.all.json).
- All predicted TFBSs for target species which are supported by at least conservation_min predictions in other species, and those supporting species, grouped into clusters (TFBSs_found.clusters.csv).
- All predicted TFBSs for target species which are supported by at least conservation_min predictions in other species, sorted by combined affinity score (TFBSs_found.sortedclusters.csv).
- Figure showing top_x_tfs highest scoring (combined affinity score) TFBSs mapped onto target_species promoter (ENSxxxxxxxxxxxx_mammals.Promoterhisto.svg). 

## 5 Species
The promoter region of any Ensembl transcript of any species within any column can be compared against the other members of the same column in order to identify a conserved binding site of the 519 transcription factors described in the Jaspar database.  The Enredo-Pecan-Ortheus pipeline was used to create whole genome alignments between the species in each column.  'EPO_LOW' indicates this column also contains genomes for which the sequencing of the current version is still considered low-coverage.  The TFBS footprinting pipeline partially accounts for this by removing sequences from alignments which appear to be missing segments.  Due to the significantly greater number of species, we recommend using the low coverage versions except for primate comparisons which do not have a low coverage version.

|EPO_LOW mammals           |EPO_LOW fish          |EPO_LOW sauropsids |EPO mammals          |EPO primates       |EPO fish              |EPO sauropsids     |
|--------------------------|----------------------|-------------------|---------------------|-------------------|----------------------|-------------------|
|ailuropoda_melanoleuca    |astyanax_mexicanus    |anas_platyrhynchos |bos_taurus           |callithrix_jacchus |danio_rerio           |anolis_carolinensis|
|bos_taurus                |danio_rerio           |anolis_carolinensis|callithrix_jacchus   |chlorocebus_sabaeus|gasterosteus_aculeatus|gallus_gallus      |
|callithrix_jacchus        |gadus_morhua          |ficedula_albicollis|canis_familiaris     |gorilla_gorilla    |lepisosteus_oculatus  |meleagris_gallopavo|
|canis_familiaris          |gasterosteus_aculeatus|gallus_gallus      |chlorocebus_sabaeus  |homo_sapiens       |oryzias_latipes       |taeniopygia_guttata|
|cavia_porcellus           |lepisosteus_oculatus  |meleagris_gallopavo|equus_caballus       |macaca_mulatta     |tetraodon_nigroviridis|                   |
|chlorocebus_sabaeus       |oreochromis_niloticus |pelodiscus_sinensis|felis_catus          |pan_troglodytes    |                      |                   |
|choloepus_hoffmanni       |oryzias_latipes       |taeniopygia_guttata|gorilla_gorilla      |papio_anubis       |                      |                   |
|dasypus_novemcinctus      |poecilia_formosa      |                   |homo_sapiens         |pongo_abelii       |                      |                   |
|dipodomys_ordii           |takifugu_rubripes     |                   |macaca_mulatta       |                   |                      |                   |
|echinops_telfairi         |tetraodon_nigroviridis|                   |mus_musculus         |                   |                      |                   |
|equus_caballus            |xiphophorus_maculatus |                   |oryctolagus_cuniculus|                   |                      |                   |
|erinaceus_europaeus       |                      |                   |ovis_aries           |                   |                      |                   |
|felis_catus               |                      |                   |pan_troglodytes      |                   |                      |                   |
|gorilla_gorilla           |                      |                   |papio_anubis         |                   |                      |                   |
|homo_sapiens              |                      |                   |pongo_abelii         |                   |                      |                   |
|ictidomys_tridecemlineatus|                      |                   |rattus_norvegicus    |                   |                      |                   |
|loxodonta_africana        |                      |                   |sus_scrofa           |                   |                      |                   |
|macaca_mulatta            |                      |                   |                     |                   |                      |                   |
|microcebus_murinus        |                      |                   |                     |                   |                      |                   |
|mus_musculus              |                      |                   |                     |                   |                      |                   |
|mustela_putorius_furo     |                      |                   |                     |                   |                      |                   |
|myotis_lucifugus          |                      |                   |                     |                   |                      |                   |
|nomascus_leucogenys       |                      |                   |                     |                   |                      |                   |
|ochotona_princeps         |                      |                   |                     |                   |                      |                   |
|oryctolagus_cuniculus     |                      |                   |                     |                   |                      |                   |
|otolemur_garnettii        |                      |                   |                     |                   |                      |                   |
|ovis_aries                |                      |                   |                     |                   |                      |                   |
|pan_troglodytes           |                      |                   |                     |                   |                      |                   |
|papio_anubis              |                      |                   |                     |                   |                      |                   |
|pongo_abelii              |                      |                   |                     |                   |                      |                   |
|procavia_capensis         |                      |                   |                     |                   |                      |                   |
|pteropus_vampyrus         |                      |                   |                     |                   |                      |                   |
|rattus_norvegicus         |                      |                   |                     |                   |                      |                   |
|sorex_araneus             |                      |                   |                     |                   |                      |                   |
|sus_scrofa                |                      |                   |                     |                   |                      |                   |
|tarsius_syrichta          |                      |                   |                     |                   |                      |                   |
|tupaia_belangeri          |                      |                   |                     |                   |                      |                   |
|tursiops_truncatus        |                      |                   |                     |                   |                      |                   |
|vicugna_pacos             |                      |                   |                     |                   |                      |                   |

