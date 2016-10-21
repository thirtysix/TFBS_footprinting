#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 2.7.0 ###########################################################
__version__ = "1.0.0b19"


# Libraries ####################################################################
import sys
import argparse
import textwrap
import os
import json
import time
import csv
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio import AlignIO

import httplib2
import math
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import*
from numpy import random as numpyrandom
from numpy import log10 as numpylog

from operator import itemgetter
import collections
from bisect import bisect_left
import itertools


################################################################################
# Functions ####################################################################
################################################################################
script_dir = os.path.dirname(__file__)
curdir = os.getcwd()

################################################################################
# Arguments ####################################################################
################################################################################

def get_args():
    """
    Retrieve arguments provided by the user.
    """

    # Instantiate the parser
    parser = argparse.ArgumentParser(
        prog="TFBS_analyzer2.py",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            TFBS Footprinting - Identification of conserved vertebrate transcription factor binding sites (TFBSs)')

            ------------------------------------------------------------------------------------------------------
            Example Usage:
                simplest:
                TFBS_analyzer2.py PATH_TO/sample_ids.txt

                all arguments:
                TFBS_analyzer2.py PATH_TO/sample_ids.txt -s homo_sapiens -g mammals -e low -pb 900 -pa 100 -l 5 -c 2 -tx 10 -o PATH_TO/Results/
            ------------------------------------------------------------------------------------------------------
            """))

    # Arguments
    parser.add_argument('t_ids_file', metavar='', type=str,
                        help='Required: Location of a file containing Ensembl target_species transcript ids (see sample file sample_ids.txt at https://github.com/thirtysix/TFBS_footprinting)")')
    parser.add_argument('--target_species', '-s', metavar='', choices = ['ailuropoda_melanoleuca', 'anas_platyrhynchos', 'anolis_carolinensis', 'astyanax_mexicanus',
                                                                         'bos_taurus', 'callithrix_jacchus', 'canis_familiaris', 'cavia_porcellus', 'chlorocebus_sabaeus',
                                                                         'choloepus_hoffmanni', 'danio_rerio', 'dasypus_novemcinctus', 'dipodomys_ordii', 'echinops_telfairi',
                                                                         'equus_caballus', 'erinaceus_europaeus', 'felis_catus', 'ficedula_albicollis', 'gadus_morhua',
                                                                         'gallus_gallus', 'gasterosteus_aculeatus', 'gorilla_gorilla', 'homo_sapiens', 'ictidomys_tridecemlineatus',
                                                                         'lepisosteus_oculatus', 'loxodonta_africana', 'macaca_mulatta', 'meleagris_gallopavo', 'microcebus_murinus',
                                                                         'mus_musculus', 'mus_spretus_spreteij', 'mustela_putorius_furo', 'myotis_lucifugus', 'nomascus_leucogenys',
                                                                         'ochotona_princeps', 'oreochromis_niloticus', 'oryctolagus_cuniculus', 'oryzias_latipes', 'otolemur_garnettii',
                                                                         'ovis_aries', 'pan_troglodytes', 'papio_anubis', 'pelodiscus_sinensis', 'poecilia_formosa', 'pongo_abelii',
                                                                         'procavia_capensis', 'pteropus_vampyrus', 'rattus_norvegicus', 'sorex_araneus', 'sus_scrofa',
                                                                         'taeniopygia_guttata', 'takifugu_rubripes', 'tarsius_syrichta', 'tetraodon_nigroviridis', 'tupaia_belangeri',
                                                                         'tursiops_truncatus', 'vicugna_pacos', 'xiphophorus_maculatus'],
                        type=str, default="homo_sapiens",
                        help='[default: "homo_sapiens"] - Target species (string), options are located at \
                        (https://github.com/thirtysix/TFBS_footprinting/blob/master/README.md#species).\
                        Conservation of TFs across other species will be based on identifying them in this species first.')
    parser.add_argument('--species_group', '-g', metavar='', choices = ["mammals", "primates", "fish", "sauropsids"], type=str, default="mammals",
                        help='("mammals", "primates", "sauropsids",  or "fish") [default: "mammals"] - Group of species (string) to identify conservation of TFs within.\
                        Your target species should be a member of this species group (e.g. "homo_sapiens" and "mammals" or "primates".\
                        The "primates" group does not have a low-coverage version.\
                        Groups and members are listed at (https://github.com/thirtysix/TFBS_footprinting/blob/master/README.md#species)')
    parser.add_argument('--coverage', '-e', metavar='',choices=("low", "high"), type=str, default="low",
                        help='("low" or "high") [default: "low"] - Which Ensembl EPO alignment of species to use.  The low coverage contains significantly more species and is recommended.\
                        The primate group does not have a low-coverage version.')
    parser.add_argument('--promoter_before_tss', '-pb', metavar='', choices = range(0, 100001), type=int, default=900,
                        help='(0-100,000) [default: 900] - Number (integer) of nucleotides upstream of TSS to include in analysis (0-100,000).')
    parser.add_argument('--promoter_after_tss', '-pa', metavar='', choices = range(0, 100001), type=int, default=100,
                        help='(0-100,000) [default: 100] - Number (integer) of nucleotides downstream of TSS to include in analysis.')
    parser.add_argument('--locality_threshold', '-l', metavar='', choices = range(0, 101), type=int, default=5,
                        help='(0-100) [default: 5] - Nucleotide distance (integer) upstream/downstream within which TF predictions in other species will be included to support a hit in the target species.')
    parser.add_argument('--conservation_min', '-c', metavar='', type=int, default=2,
                        help='(1-20)[default: 2] - Minimum number (integer) of species a predicted TF is found in, in alignment, to be considered conserved .')
    parser.add_argument('--top_x_tfs', '-tx', metavar='', choices = range(1, 21), type=int, default=10,
                        help='(1-20) [default: 10] - Number (integer) of unique TFs to include in output .svg figure.')
    parser.add_argument('--output_dir', '-o', metavar='', type=str, default=os.path.join(curdir, "tfbs_results"),
                        help=" ".join(['[default:', os.path.join(curdir, "tfbs_results"), '] - Full path of directory where result directories will be output.']))
    # Functionality to add later
    ##parser.add_argument('--pval', '-p', type=float, default=0.001, help='P-value (float) for determine score cutoff (range: 0.001 to 0.0000001) [default: 0.001]')
    ##parser.add_argument('--tf_ids', '-tfs', type=str, help='Optional: Location of a file containing a limited list of TFs to use in scoring alignment [default: all Jaspar TFs]')
    ##parser.add_argument('--noclean', '-nc', action = 'store_true', help='Optional: Don't clean retrieved alignment. Off by default.')
    args = parser.parse_args()
    transcript_ids_filename = args.t_ids_file
    species_group = args.species_group
    target_species = args.target_species
    coverage = args.coverage
    locality_threshold = args.locality_threshold
    strand_length_threshold = args.conservation_min
    promoter_before_tss = args.promoter_before_tss
    promoter_after_tss = args.promoter_after_tss
    top_x_tfs_count = args.top_x_tfs
    output_dir = args.output_dir
    
    return args, transcript_ids_filename, target_species, species_group, coverage, promoter_before_tss, promoter_after_tss, locality_threshold, strand_length_threshold, top_x_tfs_count, output_dir


################################################################################
# House-keeping functions ######################################################
################################################################################
def load_json(filename):
    with open(filename) as open_infile:
        json_data = json.load(open_infile)
    return json_data


def dump_json(filename, json_data):
    with open(filename, 'w') as open_outfile:
        json_data = json.dump(json_data, open_outfile)


def directory_creator(directory_name):
    """Create directory if it does not already exist"""
    if not os.path.isdir(directory_name):
        os.mkdir(directory_name)


def ensemblrest(query_type, options, output_type, ensembl_id=None):
    """Retrieve REST data from Ensembl using provided ID, query type, and options"""

    http = httplib2.Http()
    server = "http://rest.ensembl.org"
    
    full_query = server + query_type + ensembl_id + options
    logging.info(" ".join([time.strftime("%Y-%m-%d %H:%M:%S"), full_query]))

    if output_type == 'json':
        resp, json_data = http.request(full_query, method="GET")
        decoded_json = json.loads(json_data)
        ensemblrest_rate(resp)

        return decoded_json

    if output_type == 'fasta':
        resp, fasta_content = http.request(server+ext, method="GET", headers={"Content-Type":"text/x-fasta"})
        ensemblrest_rate(resp)

        return fasta_content


def ensemblrest_rate(resp):
    """read ensembl REST headers and determine if rate-limit has been exceeded, sleep appropriately if necessary"""
    if int(resp['x-ratelimit-remaining']) == 0:       
        if 'Retry-After' in resp:
            sleep_time = int(resp['Retry-After'])
            logging.warning(" ".join([time.strftime("%Y-%m-%d %H:%M:%S"), "Ensembl REST (Retry-After) requests sleeping for", str(sleep_time)]))
            sleep(sleep_time)
            
        else:
            sleep_time = 40
            logging.warning(" ".join([time.strftime("%Y-%m-%d %H:%M:%S"), "Ensembl REST requests sleeping for", str(sleep_time)]))
            sleep(sleep_time)


def parse_transcript_ids(transcript_ids_filename):
    """
    If user has provided a file with Ensembl transcript ids, parse these to a list.
    """

    with open(transcript_ids_filename, 'r') as transcript_ids_file:
        transcript_ids_list = transcript_ids_file.read().splitlines()
        transcript_ids_list = [x for x in transcript_ids_list if len(x)>0]

    return transcript_ids_list


################################################################################
# PWM analysis #################################################################
################################################################################
def pwm_maker(strand, motif_length, tf_motif, bg_nuc_freq_dict, neg_bg_nuc_freq_dict):
    """
    Make a PWM.
    """
    pwm = [[],[],[],[]]
    nuc_list = 'ACGT'
    # PWM according to http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647310/
    for i in range(0, motif_length):                   
        col = [tf_motif[0][i], tf_motif[1][i], tf_motif[2][i], tf_motif[3][i]]

        # number of sequences
        N = sum(col)

        # for each position (col) in the PFM.
        for j in range(0, len(tf_motif)):
            
            # nuc count at this position.
            nuc_count = tf_motif[j][i]

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
##            ppm[j].append(nuc_probability)
        
    pwm = pwm[:]

    return pwm


def PWM_scorer(seq,pwm,pwm_dict,pwm_type):
    """generate score for current seq given a pwm"""

    # set relevant variables based on whether the pwm is mono or dinucleotide
    if pwm_type == "mono":
        motif_dist = len(seq)
        span = 1

    if pwm_type == "dinuc":
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


def tfbs_finder(transcript_name, alignment, TFBS_matrix_dict, target_dir, pwm_score_threshold_dict, unaligned2aligned_index_dict, promoter_after_tss):
    """
    1.Convert PFM to PWM for each TF in the Jaspar dictionary.
    2.Score all positions in the cleaned sequence
    3.If score is greater than or equal to precomputed threshold, then keep, otherwise set to zero
    4.Return dictionary of pwm scores for [species][tf_name][strand]
    """
    start_time = time.time()
##    tfbss_found_dict_outfilename = os.path.join(target_dir, ".".join([transcript_id, "results", "all","json"]))
    tfbss_found_dict_outfilename = os.path.join(target_dir, "TFBSs_found.all.json")

    # Determine if the analysis has been done already, load results if so
    if os.path.isfile(tfbss_found_dict_outfilename):
        logging.info(" ".join([time.strftime("%Y-%m-%d %H:%M:%S"), "tfbss_found_dict already exists: loading"]))
        tfbss_found_dict = load_json(tfbss_found_dict_outfilename)
                                                    
    
    # If results don't already exist, time to party
    else:
        tfbss_found_dict = {}        
        align_chars = '-N .'
        
        for entry in alignment:
            species = entry['species']
            
            # Remove alignment/ambiguous characters from the sequences
            cleaned_seq = entry['seq']
            for char in align_chars:
                cleaned_seq = cleaned_seq.replace(char,"")
            entry_seqrecord = SeqRecord(Seq(cleaned_seq, alphabet=IUPAC.unambiguous_dna), id=species)
            forward_seq = str(entry_seqrecord.seq)
            reverse_seq = str(entry_seqrecord.seq.reverse_complement())
            seq_dict = {"+1": forward_seq, "-1":reverse_seq}

            # generate background frequencies of each mono-nucleotide for forward and reverse strands
            bg_nuc_freq_dict = {}
            neg_bg_nuc_freq_dict = {}
            nuc_list = 'ACGT'
            mononuc_pwm_dict = {"A":0,"C":1,"G":2,"T":3}
                                                    
            # calculate nucleotide frequencies for forward and reverse strands
            for nuc in nuc_list:
                bg_nuc_freq_dict[nuc] = (forward_seq.count(nuc) + 0.8) /(float(len(forward_seq)) + 0.8)
            for nuc in nuc_list:
                neg_bg_nuc_freq_dict[nuc] = (reverse_seq.count(nuc) + 0.8)/(float(len(reverse_seq)) + 0.8)


            # iterate through each tf_name and its motif
            for tf_name, tf_motif in TFBS_matrix_dict.iteritems():
                motif_length = len(tf_motif[0])
                
                if motif_length > 0:
                    # retrieve precomputed threshold and other information required for calculating the pvalue of the score
                    # set score threshold to zero if threshold less than zero
                    pwm_score_threshold_list = pwm_score_threshold_dict[tf_name][1]
                    pwm_score_threshold = pwm_score_threshold_list[0]
                    if pwm_score_threshold < 0:
                        pwm_score_threshold = 0
                    random_scores_len = pwm_score_threshold_dict[tf_name][0][0]
                    cutoff_index = pwm_score_threshold_dict[tf_name][0][1]
                    
                    # iterate through the forward and reverse strand sequences
                    for strand, seq in seq_dict.iteritems():
                        pwm = pwm_maker(strand, motif_length, tf_motif, bg_nuc_freq_dict, neg_bg_nuc_freq_dict)
                        
                        seq_length = len(seq)
                        # iterate through the nt sequence, extract a current frame based on the motif size and score
                        for i in range(0, len(seq) - motif_length):
                            current_frame = seq[i:i+motif_length]
                            current_frame_score = PWM_scorer(current_frame, pwm, mononuc_pwm_dict, 'mono')

                            # keep results that are above the precomputed threshold
                            if current_frame_score >= pwm_score_threshold:
                                current_frame_score = round(current_frame_score, 2)
                                current_frame_score_pvalue = determine_score_pvalue(current_frame_score, pwm_score_threshold_list, random_scores_len, cutoff_index)
                                hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end = start_end_found_motif(i, strand, seq_length, promoter_after_tss, motif_length)

                                # identify position in alignment from start of found motif in unaligned sequence
                                aligned_position = unaligned2aligned_index_dict[species][hit_loc_start]

                                # add to results dictionary by tf_name
                                if tf_name in tfbss_found_dict:
                                    tfbss_found_dict[tf_name].append([tf_name, species, current_frame, strand, hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end, current_frame_score, current_frame_score_pvalue, aligned_position])

                                else:
                                    tfbss_found_dict[tf_name] = [[tf_name, species, current_frame, strand, hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end, current_frame_score, current_frame_score_pvalue, aligned_position]]


        # sort results for each tf_name according to the position in alignment
        for tf_name, hits in tfbss_found_dict.iteritems():
            # ref-point
            tfbss_found_dict[tf_name] = sorted(hits, key = itemgetter(-1))
                    
        # save the results to file
        dump_json(tfbss_found_dict_outfilename, tfbss_found_dict)

    end_time = time.time()
    logging.info(" ".join([time.strftime("%Y-%m-%d %H:%M:%S"), "total time for this transcript:", str(end_time - start_time), "seconds"]))

    return tfbss_found_dict


def determine_score_pvalue(current_frame_score, pwm_score_threshold_list, random_scores_len, cutoff_index):
    """
    Determine the p-value of the PWM score of the current frame.
    Called only if score is above precomputed threshold.
    """
    # identify index of current frame score in the list of all random scores
    score_index = bisect_left(pwm_score_threshold_list, current_frame_score) + (random_scores_len + cutoff_index)

    # if the score/score_index is less than that of all random sequences, calculate pvalue
##    if score_index < len(pwm_score_threshold_list):
    if score_index < random_scores_len:

        current_frame_score_pvalue = 1 - (float(score_index)/random_scores_len) 

    # if the score/score_index equal to or greater than that of all random sequences, set minimum pvalue
    else:
        current_frame_score_pvalue = float(1)/random_scores_len

    return current_frame_score_pvalue
    
                            
def start_end_found_motif(i, strand, seq_length, promoter_after_tss, motif_length):
    """
    Determine the start/end positions of the found motif.
    """
    if strand == "+1":
        hit_loc_start = i
        hit_loc_before_TSS_start = i - seq_length + promoter_after_tss
        hit_loc_end = i + motif_length
        hit_loc_before_TSS_end = i - seq_length + motif_length + promoter_after_tss
    if strand == "-1":
        hit_loc_start = seq_length - i - motif_length
        hit_loc_before_TSS_start = (seq_length - i - motif_length) - seq_length + promoter_after_tss
        hit_loc_end = seq_length - i
        hit_loc_before_TSS_end = (seq_length - i) - seq_length + promoter_after_tss

    return hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end


def unaligned2aligned_indexes(cleaned_aligned_filename):
    """Create a dictionary for mapping aligned positions to unaligned positions"""

    with open(cleaned_aligned_filename, 'rU') as cleaned_aligned_file:
        aligned_entries_dict = SeqIO.to_dict(SeqIO.parse(cleaned_aligned_file, 'fasta'))

    unaligned2aligned_index_dict = {}

    for species, seqRecord in aligned_entries_dict.iteritems():
        unaligned2aligned_index_dict[species] = {}
        seq_str = str(seqRecord.seq)
        for aligned_index in range(len(seq_str)):
            if seq_str[aligned_index] != "-":
                unaligned_index = aligned_index - seq_str[:aligned_index].count("-")
                unaligned2aligned_index_dict[species][unaligned_index] = aligned_index

    return unaligned2aligned_index_dict


def find_clusters(target_species, tfbss_found_dict, cleaned_aligned_filename, strand_length_threshold, locality_threshold):
    """
    for each target species hit, find all hits in other species within the locality_threshold.
    identify the highest score for each species within the locality threshold
    create combined affinity score from the target species hit and those best scores from each species
    if two target species hits are within the locality threshold from one another, choose the hit which has the highest combined affinity score

    2. possibility
    """
    
    cluster_dict = {}
    for tf_name, hits in tfbss_found_dict.iteritems():
        cluster_dict[tf_name] = [] 
        clusters = []
        target_hit_seeds = []
        
        # identify all target_species hits for this tf_name
        for hit in hits:
            # ref-point
            if hit[1].lower() == target_species:
                target_hit_seeds.append(hit)

        # sort all target_species hits by some criteria: e.g. pwm score ### not necessary in this implementation
        # ref-point
        target_hit_seeds = sorted(target_hit_seeds, key = itemgetter(8), reverse=True)

        # for each target_species hit, find first all other hits in all other species within the locality threshold
        # select the best scoring hit from each species from these
        for target_hit_seed in target_hit_seeds:
            cluster = [target_hit_seed]
            target_species_strand = target_hit_seed[3]
            support_candidates = []
            cluster_added_species = []
                
            for hit in hits:
                # ref-point
                hit_strand = hit[3]
                # add a hit to support candidates if it is within the locality threshold, is on the same strand as the target_species seed, and is not a target_species seed
                if abs(target_hit_seed[-1] - hit[-1]) <= locality_threshold and hit_strand == target_species_strand and hit not in target_hit_seeds:
                    support_candidates.append(hit)

            # sort the support candidates by their scores
            support_candidates_score_sorted = sorted(support_candidates, key=itemgetter(8), reverse=True)

            # add to the cluster the best scoring support candidate from each species
            for support_candidate in support_candidates_score_sorted:
                support_candidate_species = support_candidate[1]
                if support_candidate_species not in cluster_added_species:
                    cluster.append(support_candidate)
                    cluster_added_species.append(support_candidate_species)

            clusters.append(cluster)

        # determine if each cluster meets length (strand) threshold
        for cluster in clusters:
            if len(cluster) >= strand_length_threshold:
                cluster_dict[tf_name].append(cluster)    
    
    # Calculate combined affinity score for clusters, attach to target_species hit
    for tf_name, clusters in cluster_dict.iteritems():
        for cluster in clusters:
            # ref-point
            combined_affinity_score = sum([hit[8] for hit in cluster])
            cluster[0].append(combined_affinity_score)
            # add the tfbs support (len cluster) to the target_species hit
            cluster[0].append(len(cluster))

    return cluster_dict


def cluster_table_writer(cluster_dict, target_dir, name, locality_threshold):
    """4. find best hits and output that list.
    take promoter results and output a .csv file containing all relevant information.
    only print to file those promoter hits which have support at that position."""

    output_table_name = os.path.join(target_dir, "TFBSs_found" + name  + "csv")

    with open(output_table_name, 'wb') as output_table:
        writerUS=csv.writer(output_table) 
        writerUS.writerow(['binding_prot', 'species', 'motif', 'strand', 'start', 'end', 'pre-TSS start', 'pre-TSS end', 'frame score', 'p-value', 'pos in align.', 'combined affinity score', 'support'])
        spacer_row = ["#############", "", "", "", "", "", "", "", "", "", ""]

        # for all clusters which have pass thresholds, write full cluster to .csv
        for tf_name, clusters in cluster_dict.iteritems():
            for cluster in clusters:
                for hit in cluster:
                    writerUS.writerow([str(x) for x in hit])
                writerUS.writerow(spacer_row)
                              

def target_species_hits_table_writer(sorted_clusters_target_species_hits_list, target_dir, name, locality_threshold):
    """4. find best hits and output that list.
    take promoter results and output a .csv file containing all relevant information.
    only print to file those promoter hits which have support at that position."""

    output_table_name = os.path.join(target_dir, "TFBSs_found" + name  + "csv")

    with open(output_table_name, 'wb') as output_table:
        writerUS=csv.writer(output_table) 
        writerUS.writerow(['binding_prot', 'species', 'motif', 'strand', 'start', 'end', 'pre-TSS start', 'pre-TSS end', 'frame score', 'p-value', 'pos in align.', 'combined affinity score', 'support'])

        # for all clusters which have pass thresholds, write full cluster to .csv
        for hit in sorted_clusters_target_species_hits_list:
            writerUS.writerow([str(x) for x in hit])



def sort_target_species_hits(cluster_dict):
    """
    Sort target_species hits which are part of a cluster by combined affinity score.
    """
    sorted_clusters_target_species_hits_dict = {}
    sorted_clusters_target_species_hits_list = []
    
    # create a sorted list of the target_species hits from clusters
    for tf_name, clusters in cluster_dict.iteritems():
        clusters_target_species_hits = []
            
        for cluster in clusters:
            clusters_target_species_hit = cluster[0]
            clusters_target_species_hits.append(clusters_target_species_hit)
            sorted_clusters_target_species_hits_list.append(clusters_target_species_hit)

    # ref-point     
        sorted_clusters_target_species_hits_dict[tf_name] = clusters_target_species_hits
        
    sorted_clusters_target_species_hits_list = sorted(sorted_clusters_target_species_hits_list, key=itemgetter(11), reverse = True)

    return sorted_clusters_target_species_hits_dict, sorted_clusters_target_species_hits_list


def top_x_greatest_hits(sorted_clusters_target_species_hits_list, top_x_tfs_count):
    """
    Identify the best scoring hits up to some threshold of number of tfs.
    Allows plotting more than one instance of a top tf, without increasing the total tf used count.
    e.g. 3 instances of KLF4 will count as only one tf used towards the top_x_tfs_count threshold.
    """
    # to keep track of how many tfs have been added
    top_x_tfs = []
    # to store the x greatest tfs and their locations
    top_x_greatest_hits_dict = {}


    # add all hits to single pool so top hits can be identified
    

    for sorted_clusters_target_species_hit in sorted_clusters_target_species_hits_list:
        # ref-point
        tf_name = sorted_clusters_target_species_hit[0]
        if (len(top_x_tfs) < top_x_tfs_count):

            # keep track of what & how many tfs have been added
            if tf_name not in top_x_tfs:
                top_x_tfs.append(tf_name)

            # add the hit to the top hits if the count threshold has not been met
            if tf_name in top_x_greatest_hits_dict:
                top_x_greatest_hits_dict[tf_name].append(sorted_clusters_target_species_hit)    
                
            else:
                top_x_greatest_hits_dict[tf_name] = [sorted_clusters_target_species_hit]

    return top_x_greatest_hits_dict


################################################################################
# Alignment Manipulation #######################################################
################################################################################
def retrieve_genome_aligned(species_group, target_species, chromosome, strand, promoter_start, promoter_end, coverage):
    """Takes as input target_species CCDS start position and size of promoter to be extracted.  Retrieves genome aligned,
    corresponding regions in all orthologs."""

    # Retrieve alignment if alignment FASTA does not already exist  
    query_type = "/alignment/block/region/"
    pre_options = target_species + "/" + chromosome + ":" + str(promoter_start) + "-" + str(promoter_end) + ":" + str(strand)   

    if coverage == "low":
        coverage_str = "EPO_LOW_COVERAGE"
    else:
        coverage_str = "EPO"
        
    options = pre_options + "?method=" + coverage_str + ";compact=1;content-type=application/json;species_set_group=" + species_group
    alignment_decoded = ensemblrest(query_type, options, 'json', "")

    if 'error' not in alignment_decoded:            
        # remove those entries which are computed ancestral species
        alignment = [x for x in alignment_decoded[0]['alignments'] if 'Ancestor' not in x['seq_region']]
    else:
        logging.warning(" ".join([time.strftime("%Y-%m-%d %H:%M:%S"), alignment_decoded]))
        alignment = []

    return alignment


def fasta_writer(alignment, outfile):
    """write ensembl JSON alignment to fasta file"""
    if not os.path.isfile(outfile):
        with open(outfile, "w") as aligned_file:
            for entry in alignment:
                record = SeqRecord(Seq(entry['seq'], alphabet = IUPAC.ambiguous_dna), id = entry['species'], description = "")
                SeqIO.write(record, aligned_file, 'fasta')


def remove_non_ACGT(alignment):
    """remove non alignment characters and ambiguous nucleotides.  should consider changing to replacing any non ACGT char to '-' """
    non_alignment_chars = " .N"
    for entry in alignment:
        for non_alignment_char in non_alignment_chars:
            entry['seq'] = entry['seq'].replace(non_alignment_char, '-')

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
                entry['seq'][i] = 'Z'
    for entry in alignment:
        entry['seq'] = "".join(entry['seq']).replace(u'Z',"")

    return alignment


def remove_duplicate_species(alignment):
    """
    If there are multiple entries for a single species in an alignment retrieved from Ensembl,
    keep the one which has more ACGT characters.
    """
    
    entry_ids = [x['species'] for x in alignment]
    duplicate_ids = list(set([x for x in entry_ids if entry_ids.count(x) > 1]))
    non_duplicate_alignment = [x for x in alignment if x['species'] not in duplicate_ids]
    for duplicate_id in duplicate_ids:
        duplicate_seqs = [x for x in alignment if x['species'] == duplicate_id]
        duplicate_seqs_lens = [x['seq'].count('-') for x in duplicate_seqs]
        sorted_duplicate_seqs_lens = duplicate_seqs_lens[:]
        sorted_duplicate_seqs_lens.sort()
        longest_seq = sorted_duplicate_seqs_lens[0]
        longest_seq_index = duplicate_seqs_lens.index(longest_seq)
        kept_seq = duplicate_seqs[longest_seq_index]
        if duplicate_id == target_species:
            non_duplicate_alignment = [kept_seq] + non_duplicate_alignment
        else:
            non_duplicate_alignment.append(kept_seq)     

    return non_duplicate_alignment


def selective_alignment(alignment):
    """Remove sequences from the alignment if they have less then 75% of the nucleotides of the target_species sequence."""

    target_species_entry = alignment[0]
    target_species_seq_2nd_half = target_species_entry['seq'][len(target_species_entry['seq'])/2:]
    target_species_seq_2nd_half = target_species_seq_2nd_half.replace("-","").replace("N","").replace(" ","").replace(".","")
    target_species_seq_2nd_half_len = len(target_species_seq_2nd_half)

    cleaned_alignment = []
    for entry in alignment:
        entry_seq_2nd_half = entry['seq'][len(entry['seq'])/2:]
        entry_seq_2nd_half = entry_seq_2nd_half.replace("-","").replace("N","").replace(" ","").replace(".","")      
        entry_seq_2nd_half_len = len(entry_seq_2nd_half)
        if float(entry_seq_2nd_half_len)/target_species_seq_2nd_half_len >= 0.75:
            cleaned_alignment.append(entry)

    return cleaned_alignment


def load_genome_aligned(aligned_filename):    
    """load previously retrieved alignment fasta file into dictionary"""
    with open(aligned_filename, 'r') as alignment_handle:
        alignment_list = list(SeqIO.parse(alignment_handle, 'fasta'))
    alignment = [{'seq': str(entry.seq), 'species':entry.id} for entry in alignment_list if '[' not in entry.id]
        
    return alignment


def alignment_tools(ensembl_aligned_filename, cleaned_aligned_filename, species_group, target_species, chromosome, strand, promoter_start, promoter_end, coverage):
    """
    Return cleaned alignment for further analysis.
    """

    # if cleaned alignment file doesn't exist, or the size is zero.
    if not os.path.isfile(cleaned_aligned_filename) or (os.path.isfile(cleaned_aligned_filename) and os.path.getsize(cleaned_aligned_filename) == 0):

        # If uncleaned Ensembl alignment file doesn't exist, or the size is zero: retrieve from Ensembl, write to file.
        if not os.path.isfile(ensembl_aligned_filename) or (os.path.isfile(ensembl_aligned_filename) and os.path.getsize(ensembl_aligned_filename) == 0):        
            alignment = retrieve_genome_aligned(species_group, target_species, chromosome, strand, promoter_start, promoter_end, coverage)
            fasta_writer(alignment, ensembl_aligned_filename)   

        # If uncleaned Ensembl file exists and size is not zero: clean, write to cleaned filename.
        if os.path.isfile(ensembl_aligned_filename) and (os.path.isfile(ensembl_aligned_filename) and os.path.getsize(ensembl_aligned_filename) > 0):          
            alignment = load_genome_aligned(ensembl_aligned_filename)
            alignment = remove_non_ACGT(alignment)
            alignment = remove_duplicate_species(alignment)
            alignment = selective_alignment(alignment)
            alignment = remove_gap_only(alignment)
            fasta_writer(alignment, cleaned_aligned_filename)

        # Uncleaned alignment file still doesn't exist (or size is zero): note in logfile.
        else:
            logging.warning(" ".join([time.strftime("%Y-%m-%d %H:%M:%S"), transcript_id, "- No ensembl alignment, or size is zero"]))
            alignment = []
            
    # Cleaned alignment file exists and size is not zero: load cleaned alignment.
    else:
        alignment = load_genome_aligned(cleaned_aligned_filename)

    return alignment


def transfabulator(transcript, transcript_dict_filename, promoter_before_tss, promoter_after_tss):
    """
    Given a transcript ID, retrieve Ensembl descriptive data for that transcript.
    Write json data to file.
    Based on user-defined values for target region (referenced to TSS),
    calculate genomic coordinates of target region.
    """

    # load transcript position data from json file if it already exists
    if os.path.isfile(transcript_dict_filename):
        logging.info(" ".join([time.strftime("%Y-%m-%d %H:%M:%S"), "transcript_dict already exists: loading"]))
        transcript_dict = load_json(transcript_dict_filename)

    # retrieve transcript position data from json file if it does not exist
    else:
        # Set parameters for retrieving Ensembl data via REST
        transcript_dict = {}
        query_type = '/overlap/id/'
        options = '?feature=transcript;content-type=application/json'

        # populate 'transcript_dict' dictionary with sub-dictionaries.
        # key[transcript_id] = {chromosome, strand, start, end} for each ensembl transcript id
        decoded_json_description = ensemblrest(query_type, options, 'json', transcript)
        for described_transcript in decoded_json_description:
            described_transcript_id = described_transcript['id']
            if described_transcript_id == transcript:
                transcript_dict[transcript] = described_transcript

        dump_json(transcript_dict_filename, transcript_dict)

    described_transcript = transcript_dict[transcript]

    # Extract position data
    chromosome = described_transcript['seq_region_name']
    chr_start = described_transcript['start']
    chr_end = described_transcript['end']
    strand = described_transcript['strand']
    transcript_name = described_transcript['external_name']

    if strand == 1:
        tss = chr_start
        #[promoter_start][promoter_end][TSS=chr_start][>GENE--->][chr_end]
        promoter_start = tss - promoter_before_tss
        promoter_end = tss - 1 + promoter_after_tss

    if strand == -1:
        tss = chr_end
        #[chr_start][<---GENE<][TSS=chr_end][promoter_start][promoter_end]
        promoter_start = tss + 1 - promoter_after_tss
        promoter_end = tss + promoter_before_tss

    return transcript_dict, transcript_name, chromosome, tss, strand, promoter_start, promoter_end


################################################################################
# Regulatory & Conservation Features ###########################################
################################################################################
def retrieve_regulatory(chromosome, strand, promoter_start, promoter_end, regulatory_decoded_filename, target_species):
    """Retrieve ensembl JSON data for regulatory features within the coordinates provided."""

    # determine if the regulatory data has already been retrieved, if so load, if not retrieve.
    if os.path.isfile(regulatory_decoded_filename):
        logging.info(" ".join([time.strftime("%Y-%m-%d %H:%M:%S"), "regulatory_decoded already exists: loading"]))
        regulatory_decoded = load_json(regulatory_decoded_filename)
        
    else:
        query_type = "/overlap/region/"
        pre_options = target_species + "/" + chromosome + ":" + str(promoter_start) + "-" + str(promoter_end) + ":" + str(strand)
        options = pre_options + "?feature=regulatory;content-type=application/json"
        regulatory_decoded = ensemblrest(query_type, options, 'json', "")
        dump_json(regulatory_decoded_filename, regulatory_decoded)
    
    return regulatory_decoded


def reg_position_translate(tss,regulatory_decoded,promoter_start,promoter_end,strand,promoter_before_tss,promoter_after_tss):
    """
    Convert positions of regulatory elements (json) to coordinates usable by the plot_promoter function.
    Requires supplying location relative to TSS (i.e. negative).
    For creating a left to right plot of the promoter, regardless of strand:
    Converted_reg_start is the leftmost regulatory position.
    Converted_reg_end is the rightmost regulatory position.
    """

    converted_reg_dict = {}
    for reg in regulatory_decoded:
        reg_id = reg['ID']
        reg_start = reg['start']
        reg_end = reg['end']
        description = reg['description']

        if strand == 1:
            #[promoter_start][reg_start][reg_end][promoter_end][chr_start][TSS>GENE--->][chr_end]
            converted_reg_start = (tss - reg_start) * -1
            converted_reg_end = (tss - reg_end) * -1
            if reg_start <= promoter_start:
                converted_reg_start = (-1 * promoter_before_tss + promoter_after_tss) + 0.001
            if reg_end >= promoter_end:
                converted_reg_end = promoter_after_tss - 0.001

        if strand == -1:
            #[chr_start][<---GENE<TSS][chr_end][promoter_start][reg_start][reg_end][promoter_end]
            converted_reg_start = (tss - reg_start)
            converted_reg_end = (tss - reg_end)

            if reg_start <= promoter_start:
                converted_reg_start = promoter_after_tss - 0.001        
            if reg_end >= promoter_end:
                converted_reg_end = (-1 * promoter_before_tss + promoter_after_tss) + 0.001
            
        converted_reg_dict[reg_id] = {'converted_start':converted_reg_start, 'converted_end':converted_reg_end, 'description':description}

    return converted_reg_dict


def alignment_conservation(aligned_filename):
    """Identify basic conservation of DNA sequence in the alignment."""

    alignment = AlignIO.read(aligned_filename, "fasta")
    target_species_row = alignment[0]
    species_num = float(len(alignment))
    conservation = []
    for i in range(0, len(target_species_row)):
        
        target_species_nuc = target_species_row[i]
        if target_species_nuc != "-":
            alignment_col = alignment[:,i]
            common_char = collections.Counter(alignment_col).most_common()[0][0]
            char_count = collections.Counter(alignment_col).most_common()[0][1]
            if common_char != '-':
                col_conservation = char_count/species_num
            else:
                col_conservation = alignment_col.count(target_species_nuc)/species_num
            conservation.append(col_conservation)

    return conservation


def CpG(aligned_filename):
    """Score the CpG content of the target_species sequence over a 200 nt window."""

    alignment = AlignIO.read(aligned_filename, "fasta")
    target_species_row = alignment[0]
    cpg_list = []

    # [1 C, 1 if G, 1 if CPG, CorG, num_cpg, obs2exp]
    for i in range(0, len(target_species_row)):     
        current_pos = target_species_row[i]
        if current_pos != '-':
            if i < len(target_species_row) - 1:
                next_pos = target_species_row[i+1]
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


def plot_promoter(alignment, alignment_len, promoter_before_tss, promoter_after_tss, transcript_name, top_x_greatest_hits_dict, target_dir, converted_reg_dict, conservation,cpg_list):
    """
    Plot the predicted TFBSs, onto a 5000 nt promoter graph, which possess support above the current strand threshold.
    ['binding_prot', 'species', 'motif', 'strand', 'start', 'end', 'pre-TSS start', 'pre-TSS end', 'frame score', 'p-value', 'pos in align.', 'combined affinity score', 'support']
    """

    fig = plt.figure(figsize=(10, 6))
    ax1 = plt.subplot2grid((10,1),(0,0), rowspan = 5, colspan = 11)
    ax2 = plt.subplot2grid((10,1),(6,0), sharex=ax1, rowspan = 2, colspan = 11)
    ax3 = plt.subplot2grid((10,1),(9,0), sharex=ax1, rowspan = 2, colspan = 11)

    # Generate data for each of the greatest_hits and plot corresponding bar
    color_series=['#FFB300','#803E75','#FF6800','#A6BDD7','#C10020','#CEA262','#817066','#007D34','#F6768E','#00538A','#FF7A5C','#53377A','#FF8E00','#B32851','#F4C800','#7F180D','#93AA00','#593315','#F13A13','#232C16']
    color_dict = {'CTCF':'#FF0000', 'TBP':'#FF00FF'}
    y_range = []
    labels_used = []

    # Set y-axis height based on number of entries in alignment
    tens_y = len(alignment)/10 + 1

    # Create a list sorted descending by number of supporting sequences, so that lower support hits that overlap can be seen.
    sorted_by_support_list = []
    for TF, great_hits in top_x_greatest_hits_dict.iteritems():
        for great_hit in great_hits:
            sorted_by_support_list.append(great_hit)
    # ref-point
    sorted_by_support_list = sorted(sorted_by_support_list, key=itemgetter(-1), reverse=True)

    # Choose color and plot TFs
    for x in sorted_by_support_list:
        tf_name = x[0]
        great_hit = x[1:]

        # choose a unique color for each tf_name
        if tf_name not in color_dict: 
            pick = numpyrandom.randint(0, len(color_series) - 1)
            picked_color = color_series[pick]
            color_series.remove(picked_color)
            color_dict[tf_name] = picked_color
        else:
            picked_color = color_dict[tf_name]

        # if the label has been used, set label to "", otherwise labels will repeat in legend
        if tf_name in labels_used:
            lab = ""
        else:
            lab = tf_name
            labels_used.append(tf_name)
            
        edge_color = picked_color
            
        x_series = []
        y_series = []
        # ref-point
        binding_site_start = great_hit[5]
        binding_site_end = great_hit[6]
        binding_support = great_hit[-1]
        binding_strand = int(great_hit[2])
        TF_center_point = float(binding_site_start + binding_site_end)/2
        TF_width = abs(binding_site_start - binding_site_end)
        x_series.append(TF_center_point)
        y_series.append(binding_support * binding_strand)
        y_range.append(binding_support)                    
        ax1.bar(x_series, y_series, facecolor = picked_color, edgecolor = edge_color, linewidth=1, alpha=0.9, align = 'center', width = TF_width, label = lab)

    # Plot bars for Ensembl regulatory information
    # All will be horizontally plotted in some shade of red
    if len(converted_reg_dict) > 0:
        alpha_gradient = 1.0
        alpha_gradient_dict = {1:0, 2:0.5, 3:0.25, 4:0.225}
        reg_height = 1
        reg_height = (tens_y * 1.0)/4 
        for reg_id, data in converted_reg_dict.iteritems():
            converted_start = int(data['converted_start'])
            converted_end = int(data['converted_end'])
            # limit length to first two words so legend isn't overrun
            description = data['description']
            reg_x_series = []
            reg_y_series = []
            center_point = float(converted_start + converted_end)/2
            reg_x_series.append(center_point)
            reg_y_series.append(reg_height)
            reg_x_series.append(center_point)
            reg_y_series.append(reg_height * -1)
            reg_width = abs(converted_start - converted_end)
            ax1.bar(reg_x_series, reg_y_series, facecolor='red', edgecolor='red', alpha=alpha_gradient, align = 'center', width=reg_width, label=description)
            alpha_gradient -= alpha_gradient_dict[len(converted_reg_dict)]
            reg_height += 0.5  
                            
    # Conservation plot
    ax2.plot(range(-1 * alignment_len + promoter_after_tss, promoter_after_tss), conservation, color='0.55')
    
    # CpG plot
    # [1 C, 1 if G, 1 if CPG, CorG, num_cpg, obs2exp]
    obs2exp = [x[5] for x in cpg_list]
    ax3.plot(range(-1 * alignment_len + promoter_after_tss, promoter_after_tss), obs2exp, color = 'red')
    gpc = [x[2] for x in cpg_list]
    gpc=[]
    top_obs2exp = ax3.get_ylim()[-1]
    
    for x in cpg_list:
        if x[2] == 0:
            gpc.append(x[2])
        else:
            if top_obs2exp <= 1:
                gpc.append(1)
            else:
                gpc.append(top_obs2exp)

    ax3.set_yticks(range(0,2))
    ax3.bar(range(-1 * alignment_len + promoter_after_tss, promoter_after_tss), gpc, color = 'black')

##    # Find the highest y value to set y-axis height    
    # Set format of the plot(s)
    # Hide x-ticks for first two plots
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    
    # ax1 predicted TFBSs
    num_cols = 6
    ax1.legend(bbox_to_anchor=[0., 1.05, 1.0, .102], loc='center', ncol=num_cols, prop={'size':9}, mode="expand", borderaxespad=0.)
    ax1.axhline(0, color = 'black')
    ax1.set_ylabel("Number of supporting motifs", labelpad = 0)
    ax1.set_yticks(range(-1 * ((tens_y)*10), ((tens_y)*10)+1, 10))
    title_str = transcript_name + " Predicted TFBSs"
    fig.suptitle(title_str, x=.05, rotation='vertical', fontsize=18)
    xtick_jump = 10 ** (int(numpylog(promoter_before_tss + promoter_after_tss)) - 1)
    ax3.set_xticks(range(-1 * promoter_before_tss, promoter_after_tss + 1, xtick_jump))

    # ax2 conservation
    ax2.set_yticks([0, 1])
    plt.setp(ax2.get_yticklabels(), fontsize=10)
    ax2.set_title("Conservation")

    # ax3 CpG islands
    ax3.axhline(0.6, color = 'black', alpha = 0.4)
    ax3.set_yticks([0, 0.6, 1])
    plt.setp(ax3.get_yticklabels(), fontsize=10)
    ax3.set_title("CpG Observed/Expected Ratio")  
    plt.xlabel("Nucleotide position before TSS", labelpad=10)

    plt.xlim([-1 * promoter_before_tss, promoter_after_tss + 1])
    # produce .svg figure
    fig.savefig(os.path.join(target_dir, os.path.basename(target_dir) + '.Promoterhisto'  + '.svg'), facecolor='white')
    plt.clf()
    plt.close()


################################################################################
# Initiating Variables #########################################################
################################################################################

################################################################################
# Execution ####################################################################
################################################################################

def main():
    """
    All the things.
    """

    print("Executing tfbs_footprinter version %s." % __version__)
    
    args, transcript_ids_filename, target_species, species_group, coverage, promoter_before_tss, promoter_after_tss, locality_threshold, strand_length_threshold, top_x_tfs_count, output_dir = get_args()

    # Create directory for results
    directory_creator(output_dir)

    # analysis variables
    # dictionary of thresholds for each TF
    pwm_score_threshold_dict = load_json(os.path.join(script_dir, 'data/above_cutoff_scores_dict.10000000.0.001.json'))

    # load mono-nuc PFMs
    TFBS_matrix_filename = os.path.join(script_dir, 'data/pwms.json')
    TFBS_matrix_dict = load_json(TFBS_matrix_filename)
    TFBS_matrix_dict = {k.upper():v for k,v in TFBS_matrix_dict.iteritems()}


    total_time_start = time.time()
    logging.basicConfig(filename=os.path.join(output_dir, 'TFBS_analyzer2.log'),level=logging.DEBUG)
    logging.info(" ".join([time.strftime("%Y-%m-%d %H:%M:%S"), str(vars(args))]))

    transcript_ids_list = parse_transcript_ids(transcript_ids_filename)
    for transcript_id in transcript_ids_list:
        logging.info(" ".join([time.strftime("%Y-%m-%d %H:%M:%S"), transcript_id]))
        decoded_json_description = ensemblrest('/archive/id/', '?content-type=application/json', 'json', transcript_id)

        if 'error' in decoded_json_description:
            logging.warning(" ".join([time.strftime("%Y-%m-%d %H:%M:%S"), decoded_json_description['error']]))
            continue
        
        if 'error' not in decoded_json_description:
            target_dir_name = "_".join([transcript_id, species_group])
            target_dir = os.path.join(output_dir, target_dir_name)
            directory_creator(target_dir)
            transcript_dict_filename = os.path.join(target_dir, "transcript_dict.json")
            transcript_dict, transcript_name, chromosome, tss, strand, promoter_start, promoter_end = transfabulator(transcript_id, transcript_dict_filename, promoter_before_tss, promoter_after_tss)

            # filenames for alignment and ensembl regulatory data
            ensembl_aligned_filename = os.path.join(target_dir, "alignment_uncleaned.fasta")
            cleaned_aligned_filename = os.path.join(target_dir, "alignment_cleaned.fasta")
            alignment = alignment_tools(ensembl_aligned_filename, cleaned_aligned_filename, species_group, target_species, chromosome, strand, promoter_start, promoter_end, coverage)

            if len(alignment) > 0:
                target_species_row = alignment[0]
                alignment_len = len(target_species_row['seq'].replace('-',''))

                # retrieve regulatory
                regulatory_decoded_filename = os.path.join(target_dir, "regulatory_decoded.json")
                regulatory_decoded = retrieve_regulatory(chromosome, strand, promoter_start, promoter_end, regulatory_decoded_filename, target_species)
                regulatory_decoded = load_json(regulatory_decoded_filename)
                converted_reg_dict = reg_position_translate(tss,regulatory_decoded,promoter_start,promoter_end,strand,promoter_before_tss,promoter_after_tss)

                # conservation
                conservation = alignment_conservation(cleaned_aligned_filename)
                cpg_list = CpG(cleaned_aligned_filename)

                # create index of aligned to unaligned positions
                unaligned2aligned_index_dict = unaligned2aligned_indexes(cleaned_aligned_filename)

                # score alignment for tfbss
                tfbss_found_dict = tfbs_finder(transcript_name, alignment, TFBS_matrix_dict, target_dir, pwm_score_threshold_dict, unaligned2aligned_index_dict, promoter_after_tss)

                # sort through scores, identify hits in target_species supported in other species
                cluster_dict = find_clusters(target_species, tfbss_found_dict, cleaned_aligned_filename, strand_length_threshold, locality_threshold)

                # write cluster entries to .csv
                cluster_table_writer(cluster_dict, target_dir, ".clusters.", locality_threshold)

                # sort the target_species hits supported by other species
                sorted_clusters_target_species_hits_dict, sorted_clusters_target_species_hits_list = sort_target_species_hits(cluster_dict)
                target_species_hits_table_writer(sorted_clusters_target_species_hits_list, target_dir, ".sortedclusters.", locality_threshold)
                
                # extract the top x target_species hits supported by other species
                top_x_greatest_hits_dict = top_x_greatest_hits(sorted_clusters_target_species_hits_list, top_x_tfs_count)

                # plot the top x target_species hits
                if len(top_x_greatest_hits_dict) > 0:
                    plot_promoter(alignment, alignment_len, promoter_before_tss, promoter_after_tss, transcript_name, top_x_greatest_hits_dict, target_dir, converted_reg_dict, conservation, cpg_list)


    total_time_end = time.time()
    logging.info(" ".join([time.strftime("%Y-%m-%d %H:%M:%S"), "total time for", str(len(transcript_ids_list)), "transcripts:", str(total_time_end - total_time_start), "seconds"]) + "\n\n")
