import pprint
import os
import json
from os import listdir
from os.path import isfile, join
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


################################################################################
# Functions ####################################################################
################################################################################

def dump_json(file_name, json_data):
    open_outfile =  open(file_name, 'w')
    json_data = json.dump(json_data, open_outfile)
    open_outfile.close()

def clean2Jaspar_defined(unique_sites):
    """Extract the portion of the TFBsite data that corresponds to the Jaspar Motif (uppercase letters).
    Jaspar sites contain additional sequence up/downstream of the defined TFBS, these nt are lowercase.
    Ensure that non-nucleotide characters are not included.
    
    """

    jaspar_unique_site_set = []
    error_chars = ['N', 'X', '*']
    for unique_site in unique_sites:
        cleaned_unique_set = ""
        for char in unique_site:
            if char.isupper():
                cleaned_unique_set+=char

        # make sure there are no non-nucleotide characters in binding site motif
        if not any(x in cleaned_unique_set for x in error_chars) and cleaned_unique_set != "":
            jaspar_unique_site_set.append(cleaned_unique_set)
    return jaspar_unique_site_set

def sites2dfm(sites_path, vertebrate_id2vertebrate_TF_name_dict, TFBS_likelihood_threshold):
    """Convert experimentally determined transcription factor binding sites (TFBSs) from Jaspar database,
    in 'sites' folder of current directory to dinucleotide frequency matrices DFMs
    """

    vertebrate_ids = list(vertebrate_id2vertebrate_TF_name_dict.iterkeys())
    used_ids = []
    site_files = [join(sites_path,f) for f in listdir(sites_path) if isfile(join(sites_path,f)) and ".sites" in f]
    vertebrate_file_name_paths = []
    TFBS_DiPFM_dict = {}
    TFBS_cleaned_sites_dict = {}


    for file_name in site_files:
        ID, file_extension = os.path.splitext(os.path.basename(file_name))
        if file_extension == '.sites' and ID in vertebrate_ids:
            vertebrate_file_name_paths.append(file_name)
            used_ids.append(ID)

            # access sites data, reduce to non-redundant
            unique_sites = list(set(SeqIO.parse(file_name, 'fasta')))

            # retrieve Jaspar defined TFBS (e.g. capitalized chars in 'tgcGCCCCGCCCCTcggc')
            TF_name = vertebrate_id2vertebrate_TF_name_dict[ID]
            jaspar_unique_site_set = clean2Jaspar_defined(unique_sites)


            # pre-populate dinucleotide frequency table with zeroes
            dinuc_motif_len = len(jaspar_unique_site_set[0]) - 1
            dinuc_freq = [
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len,
            [0] * dinuc_motif_len]

            # row locations of each dinuc
            dinuc_locs = {
                'AA' : 0,
                'AC' : 1,
                'AG' : 2,
                'AT' : 3,
                'CA' : 4,
                'CC' : 5,
                'CG' : 6,
                'CT' : 7,
                'GA' : 8,
                'GC' : 9,
                'GG' : 10,
                'GT' : 11,
                'TA' : 12,
                'TC' : 13,
                'TG' : 14,
                'TT' : 15}


            # iterate through each dinuc position in all cleaned experimental binding sites
            # populate dinucleotide freq table
            # i = position in list of all sites (e.g. 'ACGCGAGCCAATGGG'); j = position within each site (e.g. 'A'/'C'/'G'/'T')
            for i in range(0, len(jaspar_unique_site_set)):
                for j in range(0, dinuc_motif_len):
                    dinuc = jaspar_unique_site_set[i][j:j + 2]
                    dinuc_freq[dinuc_locs[dinuc]][j] += 1


            # test if the likelihood of finding one of the sites in this set of sites in a 1000 nt promoter region
            # if less than threshold, then add to the pfm dict and site dict
            TFBS_likelihood = ((1000.0 - len(jaspar_unique_site_set[0])) * len(set(jaspar_unique_site_set)))/(4**len(jaspar_unique_site_set[0]))
            if  TFBS_likelihood < TFBS_likelihood_threshold:

                # add dinucleotide freq table to dict
                TFBS_DiPFM_dict[TF_name] = dinuc_freq
                

                # separate entries into 'ACGT' subdicts and add to dict for use in TFBS finder script
                TFBS_cleaned_sites_dict[TF_name] = {}
                TFBS_cleaned_sites_dict[TF_name]['A'] = [x for x in set(jaspar_unique_site_set) if x[0] == 'A']
                TFBS_cleaned_sites_dict[TF_name]['C'] = [x for x in set(jaspar_unique_site_set) if x[0] == 'C']
                TFBS_cleaned_sites_dict[TF_name]['G'] = [x for x in set(jaspar_unique_site_set) if x[0] == 'G']
                TFBS_cleaned_sites_dict[TF_name]['T'] = [x for x in set(jaspar_unique_site_set) if x[0] == 'T']
                
    ##            # add the list of cleaned unique sites to dict and output .json file
    ##            TFBS_cleaned_sites_dict[TF_name] = list(set(jaspar_unique_site_set))

            else:
                print TF_name, "not included"
                print TF_name, "likelihood in 1000 nt:", TFBS_likelihood
                print TF_name, "len", len(jaspar_unique_site_set[0])
                         
    # create 1 entry dict of mammal TFBS names that cannot be made into dinuc due to lack of experimental site data
    non_dinuc_mammal_ids = list(set(vertebrate_ids) - set(used_ids))
    non_dinuc_mammal_TFBSs = {}
    non_dinuc_mammal_TFBSs['non_dinuc_mammal_TFBSs'] = [vertebrate_id2vertebrate_TF_name_dict[x] for x in non_dinuc_mammal_ids]

    # output .json file for each dict needed in the TFBS prediction script  
    # {key = TFBS name : value = duplicate-free list of experimentally determined binding sites}
    dump_json('./cleaned_Jaspar_sites.json', TFBS_cleaned_sites_dict)
    # {key = TFBS name : value = mono-nucleotide frequency matrix}
    dump_json('./non_dinuc_mammal_TFBSs.json', non_dinuc_mammal_TFBSs)
    # {key = TFBS name : value = di-nucleotide frequency matrix}
    dump_json('./dinucleotide_TFBS_PFM.json', TFBS_DiPFM_dict)

################################################################################
# Variables/Thresholds/Input ###################################################
################################################################################

# location of sites downloaded from Jaspar (http://jaspar.genereg.net/html/DOWNLOAD/sites/)
# assumed to be in 'sites' folder in current directory
sites_path = "./sites"

# threshold for likelihood of finding some TFBS in a 1000 nt sequence (ie promoter region)
# those TFBS with a likelihood lower than the threshold will be included in the output
# and used in subsequent TFBS analysis
TFBS_likelihood_threshold = 5.0

# Dictionary of Jaspar internal IDs to TFBS names
vertebrate_id2vertebrate_TF_name_dict = {
    'MA0004.1':'Arnt',
    'MA0006.1':'Arnt::Ahr',
    'MA0009.1':'T',
    'MA0017.1':'NR2F1',
    'MA0019.1':'Ddit3::Cebpa',
    'MA0025.1':'NFIL3',
    'MA0027.1':'En1',
    'MA0028.1':'ELK1',
    'MA0029.1':'Mecom',
    'MA0030.1':'FOXF2',
    'MA0031.1':'FOXD1',
    'MA0032.1':'FOXC1',
    'MA0033.1':'FOXL1',
    'MA0038.1':'Gfi1',
    'MA0040.1':'Foxq1',
    'MA0041.1':'Foxd3',
    'MA0042.1':'FOXI1',
    'MA0043.1':'HLF',
    'MA0046.1':'HNF1A',
    'MA0048.1':'NHLH1',
    'MA0051.1':'IRF2',
    'MA0056.1':'MZF1_1-4',
    'MA0057.1':'MZF1_5-13',
    'MA0059.1':'MYC::MAX',
    'MA0063.1':'Nkx2-5',
    'MA0066.1':'PPARG',
    'MA0067.1':'Pax2',
    'MA0068.1':'Pax4',
    'MA0069.1':'Pax6',
    'MA0070.1':'PBX1',
    'MA0071.1':'RORA_1',
    'MA0072.1':'RORA_2',
    'MA0073.1':'RREB1',
    'MA0074.1':'RXRA::VDR',
    'MA0075.1':'Prrx2',
    'MA0077.1':'SOX9',
    'MA0078.1':'Sox17',
    'MA0081.1':'SPIB',
    'MA0084.1':'SRY',
    'MA0087.1':'Sox5',
    'MA0088.1':'znf143',
    'MA0089.1':'NFE2L1::MafG',
    'MA0090.1':'TEAD1',
    'MA0091.1':'TAL1::TCF3',
    'MA0092.1':'Hand1::Tcfe2a',
    'MA0101.1':'REL',
    'MA0107.1':'RELA',
    'MA0108.2':'TBP',
    'MA0109.1':'Hltf',
    'MA0111.1':'Spz1',
    'MA0115.1':'NR1H2::RXRA',
    'MA0116.1':'Zfp423',
    'MA0117.1':'Mafb',
    'MA0119.1':'TLX1::NFIC',
    'MA0122.1':'Nkx3-2',
    'MA0124.1':'NKX3-1',
    'MA0125.1':'Nobox',
    'MA0130.1':'ZNF354C',
    'MA0131.1':'HINFP',
    'MA0132.1':'Pdx1',
    'MA0133.1':'BRCA1',
    'MA0135.1':'Lhx3',
    'MA0136.1':'ELF5',
    'MA0139.1':'CTCF',
    'MA0142.1':'Pou5f1::Sox2',
    'MA0149.1':'EWSR1-FLI1',
    'MA0062.2':'GABPA',
    'MA0039.2':'Klf4',
    'MA0138.2':'REST',
    'MA0002.2':'RUNX1',
    'MA0047.2':'Foxa2',
    'MA0112.2':'ESR1',
    'MA0065.2':'PPARG::RXRA',
    'MA0151.1':'ARID3A',
    'MA0152.1':'NFATC2',
    'MA0153.1':'HNF1B',
    'MA0155.1':'INSM1',
    'MA0156.1':'FEV',
    'MA0157.1':'FOXO3',
    'MA0158.1':'HOXA5',
    'MA0159.1':'RXR::RAR_DR5',
    'MA0160.1':'NR4A2',
    'MA0161.1':'NFIC',
    'MA0163.1':'PLAG1',
    'MA0164.1':'Nr2e3',
    'MA0018.2':'CREB1',
    'MA0099.2':'JUN::FOS',
    'MA0259.1':'HIF1A::ARNT',
    'MA0442.1':'SOX10',
    'MA0145.2':'Tcfcp2l1',
    'MA0146.2':'Zfx',
    'MA0141.2':'Esrrb',
    'MA0519.1':'Stat5a::Stat5b',
    'MA0518.1':'Stat4',
    'MA0517.1':'STAT2::STAT1',
    'MA0516.1':'SP2',
    'MA0515.1':'Sox6',
    'MA0514.1':'Sox3',
    'MA0513.1':'SMAD2::SMAD3::SMAD4',
    'MA0512.1':'Rxra',
    'MA0511.1':'RUNX2',
    'MA0510.1':'RFX5',
    'MA0509.1':'Rfx1',
    'MA0508.1':'PRDM1',
    'MA0507.1':'POU2F2',
    'MA0506.1':'NRF1',
    'MA0505.1':'Nr5a2',
    'MA0504.1':'NR2C2',
    'MA0503.1':'Nkx2-5 (var.2)',
    'MA0502.1':'NFYB',
    'MA0501.1':'NFE2::MAF',
    'MA0500.1':'Myog',
    'MA0499.1':'Myod1',
    'MA0498.1':'Meis1',
    'MA0497.1':'MEF2C',
    'MA0496.1':'MAFK',
    'MA0495.1':'MAFF',
    'MA0494.1':'Nr1h3::Rxra',
    'MA0493.1':'Klf1',
    'MA0492.1':'JUND (var.2)',
    'MA0491.1':'JUND',
    'MA0490.1':'JUNB',
    'MA0489.1':'JUN (var.2)',
    'MA0488.1':'JUN',
    'MA0486.1':'HSF1',
    'MA0485.1':'Hoxc9',
    'MA0484.1':'HNF4G',
    'MA0483.1':'Gfi1b',
    'MA0482.1':'Gata4',
    'MA0481.1':'FOXP1',
    'MA0480.1':'Foxo1',
    'MA0479.1':'FOXH1',
    'MA0478.1':'FOSL2',
    'MA0477.1':'FOSL1',
    'MA0476.1':'FOS',
    'MA0475.1':'FLI1',
    'MA0474.1':'Erg',
    'MA0473.1':'ELF1',
    'MA0472.1':'EGR2',
    'MA0471.1':'E2F6',
    'MA0470.1':'E2F4',
    'MA0469.1':'E2F3',
    'MA0468.1':'DUX4',
    'MA0467.1':'Crx',
    'MA0466.1':'CEBPB',
    'MA0465.1':'CDX2',
    'MA0464.1':'Bhlhe40',
    'MA0463.1':'Bcl6',
    'MA0462.1':'BATF::JUN',
    'MA0461.1':'Atoh1',
    'MA0520.1':'Stat6',
    'MA0521.1':'Tcf12',
    'MA0522.1':'Tcf3',
    'MA0523.1':'TCF7L2',
    'MA0524.1':'TFAP2C',
    'MA0525.1':'TP63',
    'MA0526.1':'USF2',
    'MA0527.1':'ZBTB33',
    'MA0528.1':'ZNF263',
    'MA0007.2':'AR',
    'MA0102.3':'CEBPA',
    'MA0024.2':'E2F1',
    'MA0154.2':'EBF1',
    'MA0162.2':'EGR1',
    'MA0076.2':'ELK4',
    'MA0258.2':'ESR2',
    'MA0098.2':'Ets1',
    'MA0148.3':'FOXA1',
    'MA0035.3':'Gata1',
    'MA0036.2':'GATA2',
    'MA0037.2':'GATA3',
    'MA0114.2':'HNF4A',
    'MA0050.2':'IRF1',
    'MA0058.2':'MAX',
    'MA0052.2':'MEF2A',
    'MA0100.2':'Myb',
    'MA0147.2':'Myc',
    'MA0104.3':'Mycn',
    'MA0150.2':'Nfe2l2',
    'MA0105.3':'NFKB1',
    'MA0060.2':'NFYA',
    'MA0014.2':'PAX5',
    'MA0080.3':'Spi1',
    'MA0143.3':'Sox2',
    'MA0079.3':'SP1',
    'MA0083.2':'SRF',
    'MA0137.3':'STAT1',
    'MA0144.2':'STAT3',
    'MA0140.2':'TAL1::GATA1',
    'MA0003.2':'TFAP2A',
    'MA0106.2':'TP53',
    'MA0093.2':'USF1',
    'MA0095.2':'YY1',
    'MA0103.2':'ZEB1',
    'MA0591.1':'Bach1::Mafk',
    'MA0592.1':'ESRRA',
    'MA0593.1':'FOXP2',
    'MA0594.1':'Hoxa9',
    'MA0595.1':'SREBF1',
    'MA0596.1':'SREBF2',
    'MA0597.1':'THAP1',
    'MA0598.1':'EHF',
    'MA0599.1':'KLF5',
    'MA0600.1':'RFX2',
    'MA0113.2':'NR3C1'}





################################################################################
# Process ######################################################################
################################################################################
sites2dfm(sites_path, vertebrate_id2vertebrate_TF_name_dict, TFBS_likelihood_threshold)







