import json
from time import sleep
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import math
import httplib2, sys
import os
import csv
import pprint
from operator import itemgetter
import matplotlib.pyplot as plt
from pylab import*
import random
from Bio import AlignIO
import collections
from bisect import bisect_left


################################################################################
# House-keeping ################################################################
################################################################################
def load_json(file_name):
    open_infile = open(file_name)
    json_data = json.load(open_infile)
    open_infile.close()
    return json_data

def dump_json(file_name, json_data):
    open_outfile =  open(file_name, 'w')
    json_data = json.dump(json_data, open_outfile)
    open_outfile.close()

def directoryCreator(directory_name):
    """Create directory if it does not already exist"""
    if not os.path.isdir('./'+directory_name):
        os.mkdir('./'+directory_name)

def occurrences(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count

################################################################################
# Sequence Retrieval ###########################################################
################################################################################

def ensembl_rest_rate(resp):
    """read ensembl REST headers and determine if rate-limit has been exceeded, sleep appropriately if necessary"""
    if int(resp['x-ratelimit-remaining']) == 0:       
        if 'Retry-After' in resp:
            sleep(resp['Retry-After'])
        else:
            sleep(40)

def ensemblrest(query_type, options, output_type, ensembl_id=None):
    """Retrieve REST data from Ensembl using provided ID, query type, and options"""
    http = httplib2.Http("ensembl_cache")
    http = httplib2.Http()
    server = "http://rest.ensembl.org"
    
    full_query = server + query_type + ensembl_id + options
    print full_query

    if output_type == 'json':
        resp, json_data = http.request(full_query, method="GET")
        pprint.pprint(resp)
        decoded_json = json.loads(json_data)
        ensembl_rest_rate(resp)
        return decoded_json

    if output_type == 'fasta':
        resp, fasta_content = http.request(server+ext, method="GET", headers={"Content-Type":"text/x-fasta"})
        pprint.pprint(resp)
        ensembl_rest_rate(resp)
        return fasta_content

def transfabulator(transcript):
    """Given a transcript ID, retrieve Ensembl data"""
    # Set parameters for retrieving Ensembl data via REST
    transcript_dict = {}

    query_type = '/overlap/id/'
    options = '?feature=transcript;content-type=application/json'

    # populate 'transcript_dict' dictionary with sub-dictionaries. key[transcript_id] = {chromosome, strand, start, end} for each ensembl transcript id
    decoded_json_description = ensemblrest(query_type, options, 'json', transcript)
    chromosome_dict = {}
    for described_transcript in decoded_json_description:
        described_transcript_id = described_transcript['id']
        if described_transcript_id == transcript:
            transcript_dict[transcript] = described_transcript

            # Extract position data
            chromosome = described_transcript['seq_region_name']
            chr_start = described_transcript['start']
            chr_end = described_transcript['end']
            strand = described_transcript['strand']
            if strand == -1:
                promoter_start = chr_end + 1
                promoter_end = chr_end + promoter_len

            if strand == 1:
                promoter_start = chr_start - promoter_len
                promoter_end = chr_start - 1


    return transcript_dict, chromosome, chr_start, chr_end, strand, promoter_start, promoter_end

def retrieve_genome_aligned(transcript_dict, promoter_len, species_group):
    """Takes as input human CCDS start position and size of promoter to be extracted.  Retrieves genome aligned, corresponding regions in all orthologs."""    

    # Retrieve alignment if alignment FASTA does not already exist
    if not os.path.isfile(aligned_file_name):
        query_type = "/alignment/block/region/"
        species = "homo_sapiens"
        pre_options = species + "/" + chromosome + ":" + str(promoter_start) + "-" + str(promoter_end) + ":" + str(strand)   
        
        options = pre_options + "?method=EPO;compact=1;content-type=application/json;species_set_group=" + species_group
        options = pre_options + "?method=EPO_LOW_COVERAGE;compact=1;content-type=application/json;species_set_group=" + species_group
        alignment_decoded = ensemblrest(query_type, options, 'json', "")

        # remove those entries which are computed ancestral species
        alignment = [x for x in alignment_decoded[0]['alignments'] if '[' not in x['species']]
    
    else:
        alignment_handle = open(aligned_file_name, 'r')
        alignment_list = list(SeqIO.parse(alignment_handle, 'fasta'))
        alignment_handle.close()
        alignment = [{'seq': str(entry.seq), 'species':entry.id} for entry in alignment_list if '[' not in entry.id]
        
    return alignment

def selective_alignment(alignment):
    """Remove sequences from the alignment if they have less then 75% of the nucleotides of the human sequence."""
    human_seq = alignment[0]['seq']
    human_seq_len = len(human_seq.replace("-","").replace(" ",""))
    cleaned_alignment = []
    for entry in alignment:
        entry_seq = entry['seq'].replace("-","").replace("N","").replace(" ","").replace(".","")
        
        entry_seq_len = len(entry_seq)
        if float(entry_seq_len)/human_seq_len >= 0.85:
            cleaned_alignment.append(entry)
        else:
            print entry['species'], "removed from alignment"

    return cleaned_alignment

def find_complexity_beginning(alignment):
    """Determine first position where more than 50% of sequences are non-gap"""
    new_alignment = []
    human_seq = alignment[0]['seq']
    alignment_length = len(human_seq)
    alignment_depth = len(alignment)
    high_complexity = False
    while high_complexity == False:
        for i in range(0, alignment_length):
            column = []
            for j in range(0, alignment_depth):
                column.append(alignment[j]['seq'][i])
            if column.count("-") <= len(alignment) * 0.5:
                high_complexity = True
                print "For analysis purposes, the start of the alignment is:", i
                break
    for entry in alignment:
        entry['seq'] = entry['seq'][i:]

    return alignment

def remove_non_ACGT(alignment):
    """remove non alignment characters and ambiguous nucleotides.  should consider changing to replacing any non ACGT char to '-' """
    for entry in alignment:
        entry['seq'] = entry['seq'].replace(' ', '-').replace('.', '-').replace('N', '-')

    return alignment

def remove_gap_only(alignment):
    """find columns in the alignment where the entire column is '-',
        replace the '-' with 'P', then remove the '*' """
    for entry in alignment:
        entry['seq'] = list(entry['seq'])

    for i in range(0,len(alignment[0]['seq'])):
        col = [x['seq'][i] for x in alignment]
        if col.count('-') == len(col):
            for entry in alignment:
                entry['seq'][i] = 'P'
    for entry in alignment:
        entry['seq'] = "".join(entry['seq']).replace(u'P',"")

    return alignment

def remove_missing_stretches(alignment):
    """ XXX Unused for now, need to account for places where an insert in one species causes all other species to be removed XXX"""
    """remove entries with more regions of '-' longer than 150 nt"""
    cleaned_alignment = []
    for entry in alignment:
        print "find - in entry", entry['seq'].count('-'*150)
        if entry['seq'].count('-'*150) == 0:
            cleaned_alignment.append(entry)

    return cleaned_alignment   
   
def fasta_writer(alignment, outfile):
    """write ensembl JSON alignment to fasta file"""

    if not os.path.isfile(outfile):
        aligned_file = open(outfile, "w")
        for entry in alignment:
            record = SeqRecord(Seq(entry['seq'], alphabet = IUPAC.ambiguous_dna), id = entry['species'], description = "")
            SeqIO.write(record, aligned_file, 'fasta')
        aligned_file.close()

def retrieve_regulatory(chromosome, strand, promoter_start, promoter_end):
    """retrieve ensembl JSON data for regulatory features within the coordinates provided"""
    query_type = "/overlap/region/"
    species = "homo_sapiens"
    pre_options = species + "/" + chromosome + ":" + str(promoter_start) + "-" + str(promoter_end) + ":" + str(strand)
    
    options = pre_options + "?feature=regulatory;content-type=application/json"
    regulatory_decoded = ensemblrest(query_type, options, 'json', "")
    print regulatory_decoded
    
    return regulatory_decoded

################################################################################
# TFBS Prediction ##############################################################
################################################################################

def PWM_scorer(seq, pwm, pwm_dict, mono):
    """generate score for current seq given a pwm"""
    # set relevant variables based on whether the pwm is mono or dinucleotide
    if mono:
        motif_dist = len(seq)
        span = 1

    if not mono:
        motif_dist = len(seq) - 1
        span = 2
    seq_score = 0.0

    # iterate through candidate sequence, and score each mono or dinucleotide
    for i in range(0, motif_dist):
        nuc = seq[i:i+span]
        row = pwm_dict[nuc]
        score = pwm[row][i]
        seq_score += score
    return seq_score

def TFBS_finder(alignment, TFBS_matrix_dict, dinuc_TFBS_matrix_dict,
                dinuc_TFBS_cleaned_sites_dict, non_dinuc_mammal_TFBSs,
                top_folder_name_input):
    """1b. Convert PFM to PWM for each TF in the Jaspar dictionary.
    Identify all possible TFBSs above a threshold of possible score using PWM.
    RULES:
    No rules found...party      
    """

    # Determine if the analysis has been done already
    TFBSs_found_dict_outfile_name = top_folder_name_input + "TFBSs_found_dict.json"
    if os.path.isfile(TFBSs_found_dict_outfile_name):
        TFBSs_found_dict = load_json(TFBSs_found_dict_outfile_name)

    # Build a PWM
    else:
        TFBSs_found_dict = {}

        # Remove alignment/ambiguous characters from the sequences
        for entry in alignment:
            align_chars = ['-', 'N', ' ', '.']
            cleaned_seq = ""
            for char in entry['seq']:
                if char not in align_chars:
                    cleaned_seq += char
            species = entry['species']
 
            print species
            start_time = time.time()
            entry_seqrecord = SeqRecord(Seq(cleaned_seq, alphabet=IUPAC.unambiguous_dna), id=species)
            forward_seq = entry_seqrecord.seq
            reverse_seq = entry_seqrecord.seq.reverse_complement()
            seq_dict = {"+1": str(forward_seq), "-1":str(reverse_seq)}

            # Di-nuc PWM 
            # Calculate ratios of each nucleotide in whole sequence       
            bg_dinuc_freq_dict = {}
            neg_bg_dinuc_freq_dict = {}
            
            dinuc_list = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
            dinuc_locs = {'AA':0, 'AC':1, 'AG':2, 'AT':3, 'CA':4, 'CC':5, 'CG':6, 'CT':7, 'GA':8, 'GC':9, 'GG':10, 'GT':11, 'TA':12, 'TC':13, 'TG':14, 'TT':15}

            for dinuc in dinuc_list:
                bg_dinuc_freq_dict[dinuc] = float(occurrences(forward_seq, dinuc) + 1)/(len(forward_seq) + 16)

            neg_bg_nuc_freq_dict = {}
            for dinuc in dinuc_list:
                neg_bg_dinuc_freq_dict[dinuc] = float(occurrences(reverse_seq, dinuc) + 1)/(len(reverse_seq) + 16)

            # iterate through the forward and reverse strand sequences
            for strand, seq in seq_dict.iteritems():

                # iterate through each of the promoter motifs, generate a best possible score using the equation.    
                # generate an empty matrix, for later population, each sublist corresponds to a dinuc.
                for TF_name, TF_motif in dinuc_TFBS_matrix_dict.iteritems():
                                        
                    pwm = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
                    ppm = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
                    best_possible_score = 0.0
                    best_possible_probability = 1
                    best_scores = []
                    motif_length = len(TF_motif[0]) + 1
                    dinuc_pfm_length = len(TF_motif[0])

                    # Create a (dinuc) PWM according to http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647310/
                    for i in range(0, dinuc_pfm_length):                   
                        col = [TF_motif[row][i] for row in range(0,16)]

                        # number of sequences
                        N = sum(col)

                        # for each position (col) in the PFM.
                        for j in range(0, len(TF_motif)):
                            
                            # nuc count at this position.
                            dinuc_count = TF_motif[j][i]

                            # pseudo-count = sqrt(total number of samples).
                            pseudo_count = 0.8

                            # background frequency for this nucleotide in the promoter.
                            if strand == "+1":
                                dinuc_bg = bg_dinuc_freq_dict[dinuc_list[j]]
                            if strand == "-1":
                                dinuc_bg = neg_bg_dinuc_freq_dict[dinuc_list[j]]

                            # probability of nuc
                            dinuc_prob = (dinuc_count + pseudo_count/16)/(N + pseudo_count)
                            dinuc_weight = math.log((dinuc_prob/dinuc_bg), 2)                       
                            pwm[j].append(dinuc_weight)
                            ppm[j].append(dinuc_prob)
                        
                    dinuc_pwm = pwm[:]

                    # Identify the greatest possible score in the PWM for this TF.
                    for i in range(0, dinuc_pfm_length):
                        col = [dinuc_pwm[row][i] for row in range(0,16)]
                        col.sort()
                        col_best_score = col[-1]
                        best_possible_score += col_best_score  

                    # score all of the experimentally validated sites, use the lowest sites score ratio as the cutoff
                    sites_scores = []
                    for char_sites in dinuc_TFBS_cleaned_sites_dict[TF_name].itervalues():
                        for site in char_sites:
                            site_score = PWM_scorer(site, dinuc_pwm, dinuc_locs, False)
                            sites_scores.append(site_score)
                    sites_scores.sort()
                    lowest_sites_score = sites_scores[0]
                    lowest_sites_score_ratio = lowest_sites_score/best_possible_score     

                    # score all possible motifs throughout the sequence
                    # iterate through the current sequence
                    seq_length = len(seq)
                    for i in range(0, seq_length - motif_length):
                        # set the current frame
                        current_frame = seq[i : i + motif_length]
                        first_nuc = current_frame[0]

                        # determine if current frame is in the dictionary of experimentally determined sites for this TF
                        if current_frame in set(dinuc_TFBS_cleaned_sites_dict[TF_name][first_nuc]):
                            current_frame_score = PWM_scorer(current_frame, dinuc_pwm, dinuc_locs, False)
                            jaspar_match = 'match'

                        # if not in experimental sites, score the frame using a pwm.
                        else:
                            current_frame_score = PWM_scorer(current_frame, dinuc_pwm, dinuc_locs, False)
                            jaspar_match = ''

                        # set the ratio of the current score to best possible score
                        current_frame_ratio = current_frame_score/best_possible_score    

                        # retain frames that are above some threshold, and add to table.
                        if current_frame_ratio >= lowest_sites_score_ratio:
                            current_frame_adjusted_score = current_frame_ratio * current_frame_score

                            if strand == "+1":
                                hit_loc_start = i
                                hit_loc_before_TSS_start = i - seq_length
                                hit_loc_end = i + motif_length
                                hit_loc_before_TSS_end = i - seq_length + motif_length
                            if strand == "-1":
                                hit_loc_start = seq_length - i - motif_length
                                hit_loc_before_TSS_start = (seq_length - i - motif_length) - seq_length
                                hit_loc_end = seq_length - i
                                hit_loc_before_TSS_end = (seq_length - i) - seq_length

                            if TF_name in TFBSs_found_dict:
                                TFBSs_found_dict[TF_name].append([species, current_frame, strand, hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end, current_frame_score, best_possible_score, current_frame_ratio, current_frame_adjusted_score, jaspar_match])
                            else:
                                TFBSs_found_dict[TF_name] = [[species, current_frame, strand, hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end, current_frame_score, best_possible_score, current_frame_ratio, current_frame_adjusted_score, jaspar_match]]

            # Mono-nuc PWM                    
            # generate background frequencies of each mono-nucleotide for forward and reverse strands
            mononucleotide_TFBS_prob = {}
            bg_nuc_freq_dict = {}
            nuc_list = ['A', 'C', 'G', 'T']
            mononuc_pwm_dict = {
                "A":0,
                "C":1,
                "G":2,
                "T":3}
            for nuc in nuc_list:
                bg_nuc_freq_dict[nuc] = forward_seq.count(nuc)/float(len(forward_seq))

            neg_bg_nuc_freq_dict = {}
            for nuc in nuc_list:
                neg_bg_nuc_freq_dict[nuc] = reverse_seq.count(nuc)/float(len(reverse_seq))

            # iterate through the forward and reverse strand sequences
            for strand, seq in seq_dict.iteritems():

                # iterate through each of the mono-nuc PFM, generate a best possible score using the equation.    
                # generate a 0 score matrix
                TF_name = "TBP"
                TF_motif = TFBS_matrix_dict[TF_name]
                cutoff_ratio = 0.5
                motif_length = len(TF_motif[0])
                if TF_name in non_dinuc_mammal_TFBSs and motif_length > 7:
                    pwm = [[],[],[],[]]
                    ppm = [[],[],[],[]]
                    best_possible_score = 0.0
                    best_possible_probability = 1
                    best_scores = []

                    # PWM according to http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647310/
                    for i in range(0, motif_length):                   
                        col = [TF_motif[0][i], TF_motif[1][i], TF_motif[2][i], TF_motif[3][i]]

                        # number of sequences
                        N = sum(col)

                        # for each position (col) in the PFM.
                        for j in range(0, len(TF_motif)):
                            
                            # nuc count at this position.
                            nuc_count = TF_motif[j][i]

                            # pseudo-count = sqrt(total number of samples).
                            pseudo_count = 0.8

                            # background frequency for this nucleotide in the promoter.
                            if strand == "+1":
                                nuc_bg = bg_nuc_freq_dict[nuc_list[j]]
                            if strand == "-1":
                                nuc_bg = neg_bg_nuc_freq_dict[nuc_list[j]]

                            # probability of nuc
                            nuc_probability = (nuc_count + pseudo_count/4)/(N + pseudo_count)

                            nuc_weight = math.log((nuc_probability/nuc_bg), 2)

                            pwm[j].append(nuc_weight)
                            ppm[j].append(nuc_probability)
                        
                    TF_motif = pwm[:]
                    
                    # Identify the greatest possible score in the PWM for this TF.
                    for i in range(0, motif_length):
                        col = [TF_motif[0][i], TF_motif[1][i], TF_motif[2][i], TF_motif[3][i]]
                        col.sort()
                        col_best_score = col[-1]
                        best_possible_score += col_best_score
                    
                    seq_length = len(seq)
                    for i in range(0, len(seq)- motif_length):
                        current_frame = seq[i:i+motif_length]
                        current_frame_score = PWM_scorer(current_frame, TF_motif, mononuc_pwm_dict, True)
                        current_frame_ratio = current_frame_score/best_possible_score

                        if current_frame_ratio >= cutoff_ratio:                            
                            current_frame_adjusted_score = current_frame_ratio * current_frame_score

                            if strand == "+1":
                                hit_loc_start = i
                                hit_loc_before_TSS_start = i - seq_length
                                hit_loc_end = i + motif_length
                                hit_loc_before_TSS_end = i - seq_length + motif_length
                            if strand == "-1":
                                hit_loc_start = seq_length - i - motif_length
                                hit_loc_before_TSS_start = (seq_length - i - motif_length) - seq_length
                                hit_loc_end = seq_length - i
                                hit_loc_before_TSS_end = (seq_length - i) - seq_length
                            if TF_name in TFBSs_found_dict:
                                TFBSs_found_dict[TF_name].append([species, current_frame, strand, hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end, current_frame_score, best_possible_score, current_frame_ratio, current_frame_adjusted_score, ""])
                            else:
                                TFBSs_found_dict[TF_name] = [[species, current_frame, strand, hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end, current_frame_score, best_possible_score, current_frame_ratio, current_frame_adjusted_score, ""]]                       

            end = time.time()
            print "total time:", end - start_time
        
        dump_json(TFBSs_found_dict_outfile_name, TFBSs_found_dict)
        
    return TFBSs_found_dict

def TFBSs_in_alignment(TFBSs_found_dict, aligned_file_name, current_strand_threshold):
    """3. Identify where the found TF binding sites are within the alignment
    (this allows for identification of supported TF binding sites among orthologs).
    Create 3 dictionaries [0. all hits sorted, 1. run-hits sorted into strands, 2. clustered hits sorted into strands]
    """

    cluster_dict = {}
    aligned_file_handle = open(aligned_file_name, 'rU')
    aligned_entries_dict = SeqIO.to_dict(SeqIO.parse(aligned_file_handle, 'fasta'))
    aligned_file_handle.close()
    TFBSs_found_dict_sorted = {}
    local_hitrun_threshold = 0
    for TF, value in TFBSs_found_dict.iteritems():
        # identify position in aligned sequence for each TFBS located in unaligned sequence
        for hit in value:
            hit_ID = hit[0]
            hit_loc = hit[3]
            seqRecord = aligned_entries_dict[hit_ID]
            hit_sequence = str(seqRecord.seq)
            aligned_position = 0
            unaligned_position = 0
            for i in range(0, len(hit_sequence)):
                if hit_sequence[i] != '-':
                    if unaligned_position == hit_loc:
                        hit.append(aligned_position)
                    unaligned_position += 1
                aligned_position += 1
                
    # sort the hits for each promoter list within the promoter dictionary by position in the alignment###
    for TF, value in TFBSs_found_dict.iteritems():
        try:
            TFBSs_found_dict_sorted[TF] = sorted(value, key = itemgetter(12))
        except:
            print "DANGER WILL ROBINSON, THAR BE SHORTNESS HERE", TF, value
    
    # newer method of clustering potential TFBSs (cluster_dict)
    # identify human hits
    for TF, values in TFBSs_found_dict_sorted.iteritems():
        cluster_dict[TF] = [] 
        clusters = []
        human_seeds = []
        for value in values:
            if 'sapiens' in value[0].lower():
                human_seeds.append(value)
        human_seeds = sorted(human_seeds, key = itemgetter(9), reverse=True)

        # for each human hit, find all values which are within the distance specified by threshold, and on the same strand
        used_human_seeds = []       
        for human_seed in human_seeds:
            if human_seed not in used_human_seeds:
                cluster = []
                cluster.append(human_seed)
                human_strand = human_seed[2]
                motif_len = len(human_seed[1])
                local_hitrun_threshold = 0
                
                for value in values:
                    value_strand = value[2]
                    if abs(human_seed[12] - value[12]) <= local_hitrun_threshold and value_strand == human_strand:
                        if human_seed != value:
                            cluster.append(value)
                        if value in human_seeds:
                            used_human_seeds.append(value)
                clusters.append(cluster)

        # determine if each cluster meets length (strand) threshold
        for cluster in clusters:
            if len(cluster) >= current_strand_threshold:
                cluster_dict[TF].append(cluster)

    return TFBSs_found_dict_sorted, cluster_dict, local_hitrun_threshold

def sort_greatest_hits(greatest_hits, num_TFs_to_include):
    """5. Sort the greatest hits by probability of occurrence, limit to the TFs_to_include"""
    # to keep track of how many TFs have been added
    top_x_TFs = []
    # to keep track of how many TFs have been added
    top_x_TFs_no_overlap = []
    # to store the x greatest TFs and their locations
    top_x_greatest_hits = {}
    # to store the x greatest TFs and their locations, keeping only the best TF in an overlap region
    top_x_greatest_hits_no_overlap = {}
    # to store a limited entry of all the entries in the greatest_hits dict
    greats = []

    # create a sorted list of the limited entries pulled from greatest_hits dict
    for TF, hits in greatest_hits.iteritems():
        if len(hits)>0:
            for hit in hits:
                temp = []
                temp = [TF] + [x for x in hit]
                greats.append(temp)
    greatest_hits_sorted = sorted(greats, key=itemgetter(6), reverse = True)

    # store the user-defined number of best TFs in a dictionary.
    # store the user-defined number of best TFs, keeping only the best scoring hit in an area of overlap, in a dictionary.
    # overlap loop
    for great in greatest_hits_sorted:
        if (len(top_x_TFs) < num_TFs_to_include) or (great[0] == 'TBP' and abs(great[1]) < 500):

            if great[0] not in top_x_TFs:
                top_x_TFs.append(great[0])

            if great[0] not in top_x_greatest_hits:
                top_x_greatest_hits[great[0]] = [great[1:]]
                
            elif great[0] in top_x_greatest_hits:
                top_x_greatest_hits[great[0]].append(great[1:])
 
    # non-overlap loop
    overlap_threshold = 0.5
    range_pool = []
    peaks_added_count = 0
    for great in greatest_hits_sorted:
        if peaks_added_count < num_TFs_to_include:
            great_start = great[1]
            great_end = great[2]
            great_range = range(great_start, great_end)
            pool_count = 0
            for x in great_range:
                if x in range_pool:
                    pool_count += 1.0
            range_pool = range_pool + great_range

            if (pool_count <= overlap_threshold * len(great_range) and peaks_added_count < num_TFs_to_include) or (great[0] == 'TBP' and abs(great[1]) < 500):
                if great[0] in top_x_greatest_hits_no_overlap:
                    top_x_greatest_hits_no_overlap[great[0]].append(great[1:])
                    peaks_added_count += 1

                if great[0] not in top_x_greatest_hits_no_overlap:
                    top_x_greatest_hits_no_overlap[great[0]] = [great[1:]]
                    peaks_added_count += 1

    return top_x_greatest_hits, greatest_hits_sorted, top_x_greatest_hits_no_overlap

def reg_position_translate(regulatory_decoded, promoter_start, promoter_end,
                           strand, promoter_len):
    """Convert positions of regulatory elements (json) to coordinates usable by the plot_promoter function.
    Requires supplying location relative to TSS (i.e. negative)"""

    converted_reg_dict = {}
    for reg in regulatory_decoded:
        reg_id = reg['ID']
        reg_start = reg['start']
        reg_end = reg['end']

        if strand == 1:
            converted_reg_start = (promoter_end - reg_start) * -1
            converted_reg_end = (promoter_end - reg_end) * -1
            if reg_start <= promoter_start:
                converted_reg_start = -1 * promoter_len
            if reg_end >= promoter_end:
                converted_reg_end = 0

        if strand == -1:
            converted_reg_start = (promoter_start - reg_start)
            converted_reg_end = (promoter_start - reg_end)
            if reg_start <= promoter_start:
                converted_reg_start = 0        
            if reg_end >= promoter_end:
                converted_reg_end = -1 * promoter_len
            
        converted_reg_dict[reg_id] = {'converted_start':converted_reg_start, 'converted_end':converted_reg_end}

    return converted_reg_dict

################################################################################
# Output data ##################################################################
################################################################################

def found_promoter_table_writer(TFBSs_found_dict_sorted, top_folder_name_input,
                                name, run_hit, human_percent_threshold,
                                local_hitrun_threshold):
    """4. find best hits and output that list.  take promoter results and output a .csv file containing all relevant information.  only print to file those promoter hits which have support at that position."""
    output_table_name = top_folder_name_input + "TFBSs-Table" + name  + "csv"
    output_table = open(output_table_name, 'wb')
    writerUS=csv.writer(output_table) 
    # 'ID' is position 0.
    writerUS.writerow(['binding_prot', 'ID', 'match', 'strand', 'start', 'end', 'pre-TSS start', 'pre-TSS end', 'frame score', 'best pos. score', 'score ratio', 'adjusted score', 'match', 'probability' , 'pos in align.', 'consensus'])
    greatest_hits = {}
    local_hitrun_threshold_backup = local_hitrun_threshold
    spacer_row = ["#############", "", "", "", "", "", "", "", "", "", "", "", "#############"]

    # Create output for cluster table
    if run_hit == 2:
        for TF, set_of_runs in TFBSs_found_dict_sorted.iteritems():
            greatest_hits[TF] = []
            for run in set_of_runs:
                run_score_total = 0
                human_hit = run[0]
                score = human_hit[7]
                score_ratio = human_hit[9]                    
                position_before_start = human_hit[5]
                position_before_start_end = human_hit[6]
                strand = human_hit[2]
                motif = human_hit[1]
                run_score_total = sum([hit[7] for hit in run])
                run_ratio_total = sum([hit[9] for hit in run])
                greatest_hits[TF].append((position_before_start, position_before_start_end, score, len(run), strand, run_score_total, score_ratio, motif))

                # write hits in cluster to cluster table
                for hit in run:
                    row = [TF] + hit
                    if 'sapiens' in hit[0]:
                         row = row + ['X']
                    writerUS.writerow(row)

                run_score_avg = run_score_total/len(run)
                writerUS.writerow(spacer_row + ["", str(len(run)), str(run_ratio_total/len(run))])

    # Create output for the Runhit table
    if run_hit == 1:
        for TF, set_of_runs in TFBSs_found_dict_sorted.iteritems():
            if TF == "Prrx2":
                local_hitrun_threshold = 1
            else:
                local_hitrun_threshold = local_hitrun_threshold_backup
            
            greatest_hits[TF] = []
            for run in set_of_runs:
                if len(run) > 0:
                    run_score_total = sum([hit[7] for hit in run])
                    run_ratio_total = sum([hit[9] for hit in run])

                    for i in range(0, len(run)):
                        hit = run[i]
                        hit_aligned_pos = hit[12]
                        species = hit[0]
                        position_before_start = hit[5]
                        position_before_start_end = hit[6]
                        strand = hit[2]
                        score_ratio = hit[9]
                        score =  hit[7]                       
                        row = [TF] + hit
                        
                        if 'homo_sapiens' in species and score_ratio >= human_percent_threshold:
                            row = row + ["X"]
                            greatest_hits[TF].append((position_before_start, position_before_start_end, score, len(run), strand, run_score_total, score_ratio))                           
                        writerUS.writerow(row)
                    run_score_avg = run_score_total/len(run)
                    writerUS.writerow(spacer_row + ["", str(len(run)), str(run_ratio_total/len(run))])

    # Create output for the larger non-runhit table
    if run_hit == 0:
        for TF, entry in TFBSs_found_dict_sorted.iteritems():
            for i in range(0, len(entry)):
                hit = entry[i]
                hit_aligned_pos = hit[12]
                if i > 0:
                    previous_hit = entry[i-1]
                    previous_hit_aligned_pos = previous_hit[12]
                    if abs(previous_hit_aligned_pos - hit_aligned_pos) > local_hitrun_threshold:
                        writerUS.writerow(spacer_row)
                
                row = [TF, hit[0], hit[1], hit[2], hit[3], hit[4], hit[5], hit[6], hit[7], hit[8], hit[9], hit[10], hit[11], hit[12]]
                writerUS.writerow(row)
    
    output_table.close() 
    return greatest_hits

def plot_promoter(identifier, transcript_id, greatest_hits, greatest_hits_sorted,
                  current_strand_threshold, promoter_len, top_folder_name_input,
                  name, num_TFs_to_include,converted_reg_dict, conservation, cpg_list):
    """6. Plot the predicted TFBSs, onto a 5000 nt promoter graph, which possess support above the current strand threshold."""
    fig = plt.figure(figsize=(12, 6))
    ax1 = plt.subplot2grid((10,1),(0,0), rowspan = 5)
    ax2 = plt.subplot2grid((10,1),(6,0), sharex=ax1, rowspan = 2)
    ax3 = plt.subplot2grid((10,1),(9,0), sharex=ax1, rowspan = 2)

    # Generate data for each of the greatest_hits and plot corresponding bar
    color_series = ['#800000', '#FF7F50', '#FF8C00', '#FFD700', '#B8860B', '#BDB76B', '#808000', '#00FF00', '#008080', '#4682B4', '#00BFFF', '#0000CD', '#8A2BE2', '#000000', '#696969', '#C0C0C0']
    color_dict = {'CTCF':'#FF0000', 'TBP':'#FF00FF'}
    y_range = []
    labels_used = []

    # Create a list sorted descending by number of supporting sequences, so that lower support hits that overlap can be seen.
    sorted_by_support_list = []
    for TF, great_hits in greatest_hits.iteritems():
        for great_hit in great_hits:
            sorted_by_support_list.append([TF] + great_hit)
    sorted_by_support_list = sorted(sorted_by_support_list, key=itemgetter(4), reverse=True)

    # Choose color and plot TFs
    for x in sorted_by_support_list:
        TF = x[0]
        great_hit = x[1:]
        if TF in labels_used:
            lab = ""
        else:
            lab = TF
            labels_used.append(TF)

        if TF not in color_dict:      
            pick = random.randint(0, len(color_series) - 1)
            picked_color = color_series[pick]
            color_series.remove(picked_color)
            color_dict[TF] = picked_color
        else:
            picked_color = color_dict[TF]
            
        print TF
        x_series = []
        y_series = []
        binding_site_start = great_hit[0]
        binding_site_end = great_hit[1]
        binding_support = great_hit[3]
        binding_strand = int(great_hit[4])

        for i in range(binding_site_start, binding_site_end):
            x_series.append(i)
            y_series.append(binding_support * binding_strand)
            y_range.append(binding_support)
                    
        ax1.bar(x_series, y_series, facecolor = picked_color, edgecolor = picked_color, align = 'center', width = 5, label = lab)

    # Plot bars for Ensembl regulatory information
    if len(converted_reg_dict) > 0:
        for reg_id, data in converted_reg_dict.iteritems():
            converted_start = data['converted_start']
            converted_end = data['converted_end']
            ax1.broken_barh([(converted_start, (converted_end - converted_start))], (-1.5,3), facecolors='red', edgecolors='', linewidths=0)
                
    # Conservation plot  
    ax2.plot(range(-1 * len(conservation), 0), conservation, color='0.55')
    
    # CPG plot
    # [1 C, 1 if G, 1 if CPG, CorG, num_cpg, obs2exp]
    obs2exp = [x[5] for x in cpg_list]
    ax3.plot(range(-1 * len(conservation), 0), obs2exp, color = 'red')
    gpc = [x[2] for x in cpg_list]
    ax3.bar(range(-1 * len(conservation), 0), gpc, color = 'red')
       
    # Find the highest y value to set y-axis height
    y_range.sort()
    highest_y = y_range[-1]
    tens_y = abs(highest_y/10)
    tens_y += 1

    # Set format of the plot(s)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
        
    num_cols = (len(greatest_hits)/2) + len(greatest_hits)%2
    ax1.legend(loc='best', ncol=num_cols, prop={'size':10})
    ax1.axhline(0, color = 'black')
    ax1.set_ylabel("Number of supporting motifs", labelpad = 0)
    ax1.set_yticks(range(-30, 31, 10))
    title_str = transcript_id + " Predicted TFBSs"
    ax1.set_title(title_str)

    ax2.set_yticks([0, 1])
    plt.setp(ax2.get_yticklabels(), fontsize=10)
    ax2.set_title("Conservation")
    
    ax3.axhline(0.6, color = 'black', alpha = 0.4)
##    ax3.set_xticks([-1000, -500, 0])
    ax3.set_yticks([0, 0.6, 1])
    plt.setp(ax3.get_yticklabels(), fontsize=10)
    ax3.set_title("CpG Observed/Expected Ratio")  

    plt.xlabel("Nucleotide position before TSS", labelpad=10)
    fig.savefig(top_folder_name_input + identifier + '.Promoterhisto' + name  + str(num_TFs_to_include) + '.' + str(promoter_len) + '.png', dpi=600, facecolor='white')

def x_greatest_hits_table_writer(top_x_greatest_hits, overlap, top_folder_name_input, name):
    """create a table of critical info for all of the top predicted TFs"""
    if overlap == True:
        output_table_name = top_folder_name_input + "Plotted_hits-Table" + name + "overlap.csv"
    if overlap == False:
        output_table_name = top_folder_name_input + "Plotted_hits-Table" + name + "no_overlap.csv"
            
    output_table = open(output_table_name, 'wb')
    writerUS=csv.writer(output_table) 
    writerUS.writerow(['binding_prot', 'start', 'end', 'score', 'support', 'strand', 'combined score', 'score ratio', 'motif'])

    dist_from_TSS_list = []
    for TF, values in top_x_greatest_hits.iteritems():
        for value in values:
            #value = [str(item) for item in value]
            dist_from_TSS_list.append([TF] + value)
            
    sorted_dist_from_TSS_list = sorted(dist_from_TSS_list, key=itemgetter(1), reverse=True)           
    for entry in sorted_dist_from_TSS_list:
        entry = [str(item) for item in entry]
        writerUS.writerow(entry)

################################################################################
# Output Correlations ##########################################################
################################################################################
        
def medisapiens_2_dict(csv_infile_name):
    """Take the medisapiens correlation table and create a dictionary where the gene names are the keys."""

    medisapiens_dict = {}                         
    csv_infile = open(csv_infile_name, 'rb')
    csv_reader = csv.reader(csv_infile, delimiter = ",")
    row_count = 0
    for row in csv_reader:
        if row_count != 0:
            corr_TF = row[1].lower()
            medisapiens_dict[corr_TF] = row[2:]
        row_count += 1

    csv_infile.close()
    return medisapiens_dict

def TF_name_cleaner(jaspar_TF):
    """Clean the list of Jaspar TFs so that they match those in the Medisapiens database."""

    cleaned_list = []    
    # Dict of JASPAR entry names to common gene names
    translate_dict = {'tcfe2a':'tcf3', 'rora_1':'rora', 'rora_2':'rora', 'zfp423':'znf423', 'tcfcp2l1':'tfcp2l1', 'mzf1_1-4':'mzf1', 'mzf1_5-13':'mzf1'}
    #special_cases = ['mzf1_5-13', 'mzf1_1-4', 'ewsr1-fli1', 'rar_dr5', 'rora_1', 'rora_2', 'gabpa', 'pou5f1', 'dux4', 'zfp423', 'hoxc9', 'e2f6', 'hnf1b', 'e2f6', 'pparg', 'hoxa5', 'sox2']

    if jaspar_TF in translate_dict:
        cleaned_list.append(translate_dict[jaspar_TF])
    else:
        cuts = jaspar_TF.split(':')
        for cut in cuts:
            if cut != "":
                cleaned_list.append(cut.split(" ")[0])
        
    return cleaned_list

def output_greatest_hits(greatest_hits_sorted, correlation_folder, top_folder_name_input, name, transcript_id):
    """Create an output table containing TFBSs sorted by combined affinity score.
    Any tissue specific MediSapiens correlation files present in the 'correlations' folder will be parsed.
    Correlations between the target gene and predicted TFBSs are included in the table.
    """

    output_table_name = top_folder_name_input + transcript_id + ".greatest_hits_sorted" + name +"csv"
    output_table = open(output_table_name, 'wb')
    writerUS=csv.writer(output_table) 
    header_row = ['binding_prot', 'start', 'end', 'score', 'support', 'strand', 'combined score', 'score ratio', 'motif']

    # create a new list which has a unique entry for each TF in combo TFBSs
    # these adjacent but separate entries will need to be regrouped in the output table
    output_greats = []
    for great in greatest_hits_sorted:
        TF = great[0].lower()
        cleaned_list = TF_name_cleaner(TF)

        for TF in cleaned_list:
            output_great = great[:]
            output_great[0] = TF
            output_greats.append(output_great)

    # retrieve any correlation files in the correlations folder
    if os.path.isdir(correlation_folder):
        correlation_files_high = [correlation_folder_name+f for f in os.listdir(correlation_folder_name) if os.path.isfile(correlation_folder_name+f) and 'high' in f and "~" not in f]
        correlation_files_low = [correlation_folder_name+f for f in os.listdir(correlation_folder_name) if os.path.isfile(correlation_folder_name+f) and 'low' in f and "~" not in f]
        correlation_files = correlation_files_high + correlation_files_low

        # Extract correlation values for each (target_protein-TF) pair for each tissue correlation file in './correlations'    
        for correlation_file in correlation_files:
            medisapiens_dict = medisapiens_2_dict(correlation_file)
            n_series = []

            for output_great in output_greats:
                TF = output_great[0].lower()
                if TF in medisapiens_dict:
                    medisapiens_corr = medisapiens_dict[TF][1]
                    medisapiens_n = medisapiens_dict[TF][5]
                    n_series.append(medisapiens_n)
                    output_great.append(medisapiens_corr)
                else:
                    output_great = great[:]
                    output_great.append('X')

            # if the first 20 correlations all have the same sample number, then display the header tissue with that sample number
            if len(set(n_series[:20])) == 1:
                n = n_series[0]
            else:
                n = 'lookup'
            file_name = correlation_file.replace(correlation_folder_name,"")
            corr_header = file_name.split(".")[1]
            corr_header = corr_header + '(' + str(n) + ')'
            header_row.append(corr_header)

    writerUS.writerow(header_row)
    for output_great in output_greats:
        writerUS.writerow(output_great)

    output_table.close()    

################################################################################
# Conservation #################################################################
################################################################################

def parse_bigfoot_pred(bigfoot_pred_infile):
    """Parse the bigfoot program's conservation output file"""

    with open(bigfoot_pred_infile) as f:
        conservation_values = f.read().splitlines()
    conservation_values = [float(x) for x in conservation_values]
    return conservation_values

def alignment_conservation(aligned_file_name):
    """Identify basic conservation of DNA sequence in the alignment."""

    alignment = AlignIO.read(aligned_file_name, "fasta")
    human_row = alignment[0]
    species_num = float(len(alignment))
    conservation = []
    for i in range(0, len(human_row)):
        
        human_nuc = human_row[i]
        if human_nuc != "-":
            alignment_col = alignment[:,i]
            common_char = collections.Counter(alignment_col).most_common()[0][0]
            char_count = collections.Counter(alignment_col).most_common()[0][1]
            if common_char != '-':
                col_conservation = char_count/species_num
            else:
                col_conservation = alignment_col.count(human_nuc)/species_num
            conservation.append(col_conservation)

    return conservation     

def CpG(aligned_file_name):
    """Score the CpG content of the human sequence over a 200 nt window."""

    alignment = AlignIO.read(aligned_file_name, "fasta")
    human_row = alignment[0]
    cpg_list = []

    # [1 C, 1 if G, 1 if CPG, CorG, num_cpg, obs2exp]
    for i in range(0, len(human_row)):     
        current_pos = human_row[i]
        if current_pos != '-':
            if i < len(human_row) - 1:
                next_pos = human_row[i+1]
            else:
                next_pos = False
             
            if current_pos == 'C' and next_pos == 'G':
                cpg_list.append([1, 0, 1])
                
            elif current_pos == 'C' and next_pos != 'G':
                cpg_list.append([1, 0, 0])

            elif current_pos == 'G':
                cpg_list.append([0, 1, 0])

            else:
                cpg_list.append([0, 0, 0])

    for i in range(0, len(cpg_list)):
        if i < 100:
            rolling_island = cpg_list[:i] + cpg_list[i:i+100]
        elif i > len(cpg_list) - 100:
            rolling_island = cpg_list[i - 100:i] + cpg_list[i:]
        else:
            rolling_island = cpg_list[i-100:i+100]

        Cs = sum([x[0] for x in rolling_island]) * 1.0
        Gs = sum([x[1] for x in rolling_island]) * 1.0
        CorG_ratio = (Cs+Gs)/len(rolling_island)
        num_cpg = sum([x[2] for x in rolling_island]) * 1.0
        obs = num_cpg/len(rolling_island)
        exp = (CorG_ratio/2)**2
        obs2exp = obs/exp
        cpg_list[i] = cpg_list[i] + [CorG_ratio, num_cpg, obs2exp]
        
    return cpg_list

################################################################################
# Variables/Thresholds/Input ###################################################
################################################################################
# load mono-nuc PFMs
TFBS_matrix_file_name = './TFBS_matrix.json'
TFBS_matrix_dict = load_json(TFBS_matrix_file_name)

# PFMs for each dinuc TF
dinuc_TFBS_matrix_dict = load_json('./dinucleotide_TFBS_PFM.json')

# Unique sites for each dinuc TF
dinuc_TFBS_cleaned_sites_dict = load_json('./cleaned_Jaspar_sites.json')

# names of TFs for which there is no site data, only PFM
non_dinuc_mammal_json = load_json('./non_dinuc_mammal_TFBSs.json')
non_dinuc_mammal_TFBSs = non_dinuc_mammal_json['non_dinuc_mammal_TFBSs']

# target Ensembl transcripts
# e.g. transcripts_dict = {'CA12-001':'ENST00000178638', 'CA12-002':'ENST00000344366', 'CA12-003':'ENST00000422263'}
transcripts_dict = {} 


# length of promoter to analyze
promoter_len = 1000

# Ensembl species group (http://rest.ensembl.org/documentation/info/compara_species_sets)
# (e.g. "fish", "sauropsids", "mammals", "primates")
species_group = "mammals"

# threshold for accepting predicted TFBSs present only in mono-nucleotide format
# (currently only TATA-Binding Protein (TBP))
human_percent_threshold = 0.50

# threshold for number of supporting species required to identify a predicted
# human TFBS as conserved
current_strand_threshold = 2

# number of TFs to include in greatest hits list and subsequent graph
num_TFs_to_include = 10

# (OPTIONAL) if bigfoot program has been used and produced phylogenetic footprinting values
bigfoot_pred_infile = ""

# (OPTIONAL) if expression correlation files have been retrieved from ist.medisapiens.com
correlation_folder = "./correlations"

################################################################################
# Process ######################################################################
################################################################################
for transcript_id, transcript in transcripts_dict.iteritems():
    identifier = transcript_id + '_' + transcript + "_" + species_group
    top_folder_name_input = "./" + identifier + '/'
    directoryCreator(top_folder_name_input)

    # retrieve transcript data (positions)
    transcript_dict, chromosome, chr_start, chr_end, strand, promoter_start, promoter_end = transfabulator(transcript)

    # retrieve alignment
    aligned_file_name = top_folder_name_input + "cleaned_alignment.fasta"
    ensembl_aligned_file_name = top_folder_name_input + "uncleaned_alignment.fasta"
    alignment = retrieve_genome_aligned(transcript_dict, promoter_len, species_group)
    fasta_writer(alignment, ensembl_aligned_file_name)

    # clean alignment
    alignment = find_complexity_beginning(alignment)
    alignment = selective_alignment(alignment)
    alignment = remove_non_ACGT(alignment)
    alignment = remove_gap_only(alignment)
    fasta_writer(alignment, aligned_file_name)
    
    # retrieve regulatory
    regulatory_decoded = retrieve_regulatory(chromosome, strand, promoter_start, promoter_end)
    converted_reg_dict = reg_position_translate(regulatory_decoded, promoter_start, promoter_end, strand, promoter_len)

    # identify TFBSs
    TFBSs_found_dict = TFBS_finder(alignment, TFBS_matrix_dict, dinuc_TFBS_matrix_dict, dinuc_TFBS_cleaned_sites_dict, non_dinuc_mammal_TFBSs, top_folder_name_input)
    TFBSs_found_dict_sorted, cluster_dict, local_hitrun_threshold = TFBSs_in_alignment(TFBSs_found_dict, aligned_file_name, current_strand_threshold)
    greatest_hits_all = found_promoter_table_writer(TFBSs_found_dict_sorted, top_folder_name_input, ".all.", 0, human_percent_threshold, local_hitrun_threshold)   

    # conservation
    # raw conservation scores
    conservation = alignment_conservation(aligned_file_name)
    # bigfoot conservation scores
    if os.path.isfile(bigfoot_pred_infile):
        conservation = parse_bigfoot_pred(bigfoot_pred_infile)
    # cpg
    cpg_list = CpG(aligned_file_name)

    # cluster identified hits conserved across multiple species
    greatest_hits_cluster = found_promoter_table_writer(cluster_dict, top_folder_name_input, ".clusters.", 2, human_percent_threshold, local_hitrun_threshold)
    top_x_greatest_hits, greatest_hits_sorted, top_x_greatest_hits_no_overlap = sort_greatest_hits(greatest_hits_cluster, num_TFs_to_include)

    # write output figures and tables
    plot_promoter(identifier, transcript_id, top_x_greatest_hits, greatest_hits_sorted, current_strand_threshold, promoter_len, top_folder_name_input, ".cluster.overlap.", num_TFs_to_include, converted_reg_dict, conservation, cpg_list)
    plot_promoter(identifier, transcript_id, top_x_greatest_hits_no_overlap, greatest_hits_sorted, current_strand_threshold, promoter_len, top_folder_name_input, ".cluster.NO_overlap.", num_TFs_to_include, converted_reg_dict, conservation, cpg_list)
    x_greatest_hits_table_writer(top_x_greatest_hits, True, top_folder_name_input, ".cluster.")
    x_greatest_hits_table_writer(top_x_greatest_hits_no_overlap, False, top_folder_name_input, ".cluster")
    output_greatest_hits(greatest_hits_sorted, correlation_folder, top_folder_name_input, ".cluster.", transcript_id)
