#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 2.7.0 ###########################################################
__version__ = "1.0.0b46"


# Libraries ####################################################################
import sys
import signal
import wget
import tarfile
import argparse
import textwrap
import os
import json
import msgpack
import time
import csv
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import AlignIO

import socket
import httplib2
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import mpl
from numpy import random as numpyrandom
from decimal import Decimal

from operator import itemgetter
from bisect import bisect_left
from bisect import bisect_right

################################################################################
# Description ##################################################################
################################################################################
"""
Thresholds are allowed to be negative, as based on whole genome scoring.

Improvements to be made:
1) Change FANTOM CAGE correlation analysis to only include top CAGE peaks,
instead of all CAGES.  Should reduce size and focus on most relevant CAGEs.
2) Account for GTRD metacluster peak count.  Switch to GTRD peaks by TF if the
number of TFs available becomes larger than those offered by JASPAR.  Similarly,
if all JASPAR TFs become present in GTRD database, then switch from metaclusters
to peaks by TF.
"""


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
        prog="tfbs_footprinter",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            TFBS Footprinting - Identification of conserved vertebrate transcription factor binding sites (TFBSs).
            See https://github.com/thirtysix/TFBS_footprinting for additional usage instructions.

            ------------------------------------------------------------------------------------------------------
            Example Usage:
            
                simplest:
                tfbs_footprinter PATH_TO/sample_ensembl_ids.txt

                all arguments:
                tfbs_footprinter -t PATH_TO/sample_ensembl_ids.txt -tfs PATH_TO/sample_jaspar_tf_ids.txt -pb 900 -pa 100 -tx 10 -p 0.01 -update

                run the sample analysis:
                Option #1: tfbs_footprinter -t PATH_TO/sample_analysis/sample_analysis_list.csv
                Option #2: tfbs_footprinter -t PATH_TO/sample_analysis/sample_ensembl_ids.txt

                update the experimental data files (not needed often):
                tfbs_footprinter -update

                Results will be output to the current directory in a created directory named "tfbs_results"
            ------------------------------------------------------------------------------------------------------
            """))

##    # retrieve most recent list of available Ensembl species
##    available_species_query_type_low_coverage = "/info/compara/species_sets/EPO_low_coverage"
##    available_species_query_type = "/info/compara/species_sets/EPO"
##    available_species_options = "?content-type=application/json"
##
##    available_species_dict_list_low_coverage = [d["species_set"] for d in ensemblrest(available_species_query_type_low_coverage, available_species_options, "json", "", log=False)]
##    available_species_dict_list = [d["species_set"] for d in ensemblrest(available_species_query_type, available_species_options, "json", "", log=False)]
##    combined_available_species = [sp for sp_sl in available_species_dict_list_low_coverage + available_species_dict_list for sp in sp_sl]
##    available_species = list(set(combined_available_species))
##    available_species.sort()
    
    # species groups
    species_groups = ["mammals", "primates", "sauropsids", "fish"]

    # Arguments
##    parser.add_argument('t_ids_file', metavar='', type=str,
##                        help='Required: Location of a file containing Ensembl target_species transcript ids.  Input options are either a text file of Ensembl transcript ids or a .csv file with individual values set for each parameter.')
    parser.add_argument('--t_ids_file', '-t', metavar='', type=str,
                        help='Required for running an analysis.  Location of a file containing Ensembl target_species transcript ids.  Input options are either a text file of Ensembl transcript ids or a .csv file with individual values set for each parameter.')    
    parser.add_argument('--tf_ids_file', '-tfs', metavar='', type=str, default = None, help='Optional: Location of a file containing a limited list of Jaspar TFs to use in scoring alignment \
                                                                                                (see sample file tf_ids.txt at https://github.com/thirtysix/TFBS_footprinting) [default: all Jaspar TFs]')
##    parser.add_argument('--target_species', '-s', metavar='', choices = available_species,
##                        type=str, default="homo_sapiens",
##                        help='[default: "homo_sapiens"] - Target species (string), options are located at \
##                        (https://github.com/thirtysix/TFBS_footprinting/blob/master/README.md#species).\
##                        Conservation of TFs across other species will be based on identifying them in this species first.')
    
##    parser.add_argument('--species_group', '-g', metavar='', choices = species_groups, type=str, default="mammals",
##                        help='("mammals", "primates", "sauropsids",  or "fish") [default: "mammals"] - Group of species (string) to identify conservation of TFs within.\
##                        Your target species should be a member of this species group (e.g. "homo_sapiens" and "mammals" or "primates").\
##                        The "primates" group does not have a low-coverage version.\
##                        Groups and members are listed at (https://github.com/thirtysix/TFBS_footprinting/blob/master/README.md#species)')
##
##    parser.add_argument('--coverage', '-e', metavar='',choices=("low", "high"), type=str, default="low",
##                        help='("low" or "high") [default: "low"] - Which Ensembl EPO alignment of species to use.  The low coverage contains significantly more species and is recommended.\
##                        The primate group does not have a low-coverage version.')

    parser.add_argument('--promoter_before_tss', '-pb', metavar='', choices = range(-10000, 100001), type=int, default=900,
                        help='(0-100,000) [default: 900] - Number (integer) of nucleotides upstream of TSS to include in analysis.  If this number is negative the start point will be downstream of the TSS, the end point will then need to be further downstream.')

    parser.add_argument('--promoter_after_tss', '-pa', metavar='', choices = range(-10000, 100001), type=int, default=100,
                        help='(0-100,000) [default: 100] - Number (integer) of nucleotides downstream of TSS to include in analysis.  If this number is negative the end point will be upstream of the TSS.  The start point will then need to be further upstream.')

    parser.add_argument('--top_x_tfs', '-tx', metavar='', choices = range(1, 21), type=int, default=10,
                        help='(1-20) [default: 10] - Number (integer) of unique TFs to include in output .svg figure.')

##    parser.add_argument('--output_dir', '-o', metavar='', type=str, default=os.path.join(curdir, "tfbs_results"),
##                        help=" ".join(['[default:', os.path.join(curdir, "tfbs_results"), '] - Full path of directory where result directories will be output.  Make sure that the root directory already exists.']))

    # for now pvalue refers to the PWM score, in the future it will need to relate to the combined affinity score
    parser.add_argument('--pval', '-p', type=float, default=0.01, help='P-value (float) for determine score cutoff (range: 0.1 to 0.0000001) [default: 0.01]')

    parser.add_argument('--exp_data_update', '-update', action="store_true", help='Download the latest experimental data files for use in analysis.  Will run automatically if the "data" directory does not already exist (e.g. first usage).')

    # Functionality to add later
    ##parser.add_argument('--noclean', '-nc', action = 'store_true', help='Optional: Don't clean retrieved alignment. Off by default.')

    # pre-processing the arguments
    args = parser.parse_args()
    args_lists = []
    transcript_ids_filename = args.t_ids_file
    exp_data_update = args.exp_data_update

    if transcript_ids_filename:
        filename, file_extension = os.path.splitext(transcript_ids_filename)

        if file_extension == ".csv" or file_extension == ".tsv":
            # If the user has provided a .csv file with the required parameters defined for each Ensembl transcript id
            # this can be parsed to run unique analyses for each.
            
            parsed_arg_lines = file_to_datalist(transcript_ids_filename)[1:]
            
            for i, parsed_arg_line in enumerate(parsed_arg_lines):
                if len(parsed_arg_line) < 6:
                    print "Incomplete arguments in input file on line", i
                    
                else:
                    transcript_id, target_tfs_filename, promoter_before_tss, promoter_after_tss, top_x_tfs_count, pval = parsed_arg_line

                    # promoter_before_tss/promoter_after_tss
                    try:
                        promoter_before_tss = int(promoter_before_tss)
                    except:
                        print "Entered promoter before TSS", promoter_before_tss, "in line", i, "is not an integer.  Defaulting to 900."
                        promoter_before_tss = 900

                    try:
                        promoter_after_tss = int(promoter_after_tss)
                    except:
                        print "Entered promoter after TSS", promoter_after_tss, "in line", i, "is not an integer.  Defaulting to 100."
                        promoter_after_tss = 100

                    # top_x_tfs_count
                    try:
                        top_x_tfs_count = int(top_x_tfs_count)
                    except:
                        print "Entered top x tfs count", top_x_tfs_count, "in line", i, "is not an integer.  Defaulting to 10."
                        top_x_tfs_count = 10

                    # p-value
                    try:
                        pval = float(pval)
                    except:
                        print "Entered p-value threshold", pval, "in line", i, "is not float.  Defaulting to 0.01."
                        pval = 0.01

                    # update exp data
                    exp_data_update = False
                    
                    parsed_cleaned_arg_line = [transcript_id, target_tfs_filename, promoter_before_tss, promoter_after_tss, top_x_tfs_count, pval]
                    args_lists.append([args, transcript_ids_filename] + parsed_cleaned_arg_line)

        else:
            # If the analysis does not require setting the parameters individually for each Ensembl transcript id then build
            # build a list which has all of the parameters set as the same, in this way there can be a similar input format
            # as a .tsv, and standardized handling in the rest of the analysis.
            target_tfs_filename = args.tf_ids_file
            promoter_before_tss = args.promoter_before_tss
            promoter_after_tss = args.promoter_after_tss
            top_x_tfs_count = args.top_x_tfs
            pval = args.pval
            exp_data_update = args.exp_data_update
            
            transcript_ids_list = parse_transcript_ids(transcript_ids_filename)
            for transcript_id in transcript_ids_list:
                args_list = [args, transcript_ids_filename, transcript_id, target_tfs_filename, promoter_before_tss, promoter_after_tss, top_x_tfs_count, pval]
                args_lists.append(args_list)

    return args_lists, exp_data_update


################################################################################
# House-keeping functions ######################################################
################################################################################
def signal_handler(signal, frame):
        print('You have manually stopped tfbs_footprinter with Ctrl+C')
        sys.exit(0)


def load_json(filename):
    if os.path.exists(filename):
        with open(filename) as open_infile:
            return json.load(open_infile)
    else:
        return None
        
    return json_data


def dump_json(filename, json_data):
    with open(filename, 'w') as open_outfile:
        json_data = json.dump(json_data, open_outfile)


def load_msgpack(object_filename):
    """unpack a msgpack file to object."""

    if os.path.exists(object_filename):
        with open(object_filename, 'r') as object_file:
            return msgpack.unpack(object_file, max_array_len=200000, use_list=False)
    else:
        return None


def save_msgpack(msgpack_obj, msgpack_filename):
    """Save msgpack object to file."""
    
    with open(msgpack_filename, 'w') as msgpack_file:
        msgpack.pack(msgpack_obj, msgpack_file, use_bin_type=True) 


def directory_creator(directory_name):
    """
    Create directory if it does not already exist.
    """

    if not os.path.isdir(directory_name):
        os.mkdir(directory_name)


def is_online():
    """
    Test if the system is online.
    This breaks when TFBS_footprinter outlasts Google.
    """

    REMOTE_SERVER = "www.google.com"
    try:
        host = socket.gethostbyname(REMOTE_SERVER)
        s = socket.create_connection((host, 80), 2)
        return True

    except:
        logging.info(" ".join(["System does not appear to be connected to the internet."]))
        return False
        

def ensemblrest(query_type, options, output_type, ensembl_id=None, log=False):
    """
    Retrieve REST data from Ensembl using provided ID, query type, and options.
    """

    http = httplib2.Http()
    server = "http://rest.ensembl.org"
    
    full_query = server + query_type + ensembl_id + options

    if log:
        logging.info(" ".join(["Ensembl REST query made:", full_query]))

    success = False
    try_count = 0
    max_tries = 10
    fail_sleep_time = 120
    decoded_json = {}
    fasta_content = []

    if output_type == 'json':
        if success == False and try_count < max_tries:
            try:
                resp, json_data = http.request(full_query, method="GET")
                decoded_json = json.loads(json_data)
                ensemblrest_rate(resp)
                try_count += 1
                success = True
                return decoded_json

            except:
                logging.info(" ".join(["Ensembl REST query unsuccessful, attempt:", "/".join([str(try_count), str(max_tries)]), "Sleeping:", str(fail_sleep_time), "seconds.", full_query]))
                try_count += 1
                print " ".join(["Ensembl REST query unsuccessful, attempt:", "/".join([str(try_count), str(max_tries)]), "Sleeping:", str(fail_sleep_time), "seconds.", "See logfile for query."])
                time.sleep(fail_sleep_time)

        # return empty decoded_json if max tries has elapsed
        return decoded_json
    
    if output_type == 'fasta':
        if success == False and try_count < max_tries:
            try:
                resp, fasta_content = http.request(server, method="GET", headers={"Content-Type":"text/x-fasta"})
                ensemblrest_rate(resp)
                try_count += 1
                success = True
                return fasta_content

            except:    
                logging.info(" ".join(["Ensembl REST query unsuccessful, attempt:", "/".join([str(try_count), str(max_tries)]), "Sleeping:", str(fail_sleep_time), "seconds.", full_query]))
                try_count += 1
                print " ".join(["Ensembl REST query unsuccessful, attempt:", "/".join([str(try_count), str(max_tries)]), "Sleeping:", str(fail_sleep_time), "seconds.", "See logfile for query."])
                time.sleep(fail_sleep_time)

        # return empty fasta_content if max tries has elapsed
        return fasta_content


def ensemblrest_rate(resp):
    """
    Read ensembl REST headers and determine if rate-limit has been exceeded, sleep appropriately if necessary.
    """

    if int(resp['x-ratelimit-remaining']) == 0:       
        if 'Retry-After' in resp:
            sleep_time = int(resp['Retry-After'])
            logging.warning(" ".join(["Ensembl REST (Retry-After) requests sleeping for:", str(sleep_time)]))
            time.sleep(sleep_time)
            
        else:
            sleep_time = 60
            logging.warning(" ".join(["Ensembl REST requests sleeping for:", str(sleep_time)]))
            time.sleep(sleep_time)


def parse_transcript_ids(transcript_ids_filename):
    """
    If user has provided a file with Ensembl transcript ids, parse these to a list.
    """

    with open(transcript_ids_filename, 'r') as transcript_ids_file:
        transcript_ids_list = transcript_ids_file.read().splitlines()
        transcript_ids_list = [x for x in transcript_ids_list if len(x)>0]

    return transcript_ids_list


def parse_tf_ids(target_tfs_filename):
    """
    If user has provided a file with Ensembl transcript ids, parse these to a list.
    """

    with open(target_tfs_filename, 'r') as target_tfs_file:
        target_tfs_list = target_tfs_file.read().splitlines()
        target_tfs_list = [x.upper() for x in target_tfs_list if len(x)>0]

    return target_tfs_list


def file_to_datalist(data_filename):
    """
    Starting with a filename, import and convert data to a list.
    Attempted to use csv.Sniffer() but it was inconsistent in delimiter detection.
    Instead, limit users to comma or tab-separated input files.
    """

    # accept only tab
    filename, file_extension = os.path.splitext(data_filename)
    file_extension = file_extension.lower()
    
    ext_delimiter_dict = {".csv":",", ".tsv":"\t"}
    
    with open(data_filename, 'r') as data_file:
        csv_reader = csv.reader(data_file, delimiter = ext_delimiter_dict[file_extension])
        all_data = list(csv_reader)

    return all_data
        

def compare_tfs_list_jaspar(target_tfs_list, TFBS_matrix_dict):
    """
    If user has provided a file containing Jaspar TF ids,
    compare candidate entries to those in the loaded dictionary of Jaspar PWMs.
    """

    jaspar_dict_keys = TFBS_matrix_dict.keys()
    erroneous = list(set(target_tfs_list) - set(jaspar_dict_keys))

    target_tfs_list = list(set(jaspar_dict_keys).intersection(target_tfs_list))
    if len(erroneous) > 0:
        logging.warning(" ".join(["the following tf ids are not in the Jaspar database:", ", ".join(erroneous)]))

    return target_tfs_list


def experimentalDataUpdater(exp_data_update):
    """
    Update the experimental data by downloading it from the Amazon repository.
    Only activates if the user specifically calls for an update, or the data directory does not exist.
    """

    experimental_data_dir = os.path.join(script_dir, 'data')
    experimental_data_present = False

    # test if data dir exists
    if not os.path.exists(experimental_data_dir):
        directory_creator(experimental_data_dir)
        exp_data_update = True
        print "Data dir doesn't exist"

    # if the data dir exists, check to see that all of the required file patterns are present
    else:
        required_data_file_patterns = ["pwms.json", "all_tfs_thresholds", "all_pwms_loglikelihood_dict"]
        experimental_data_filenames = [x for x in os.listdir(experimental_data_dir) if os.path.isfile(os.path.join(experimental_data_dir, x))]
        all_patterns_matched = all([any([required_data_file_pattern in experimental_data_filename
                                          for experimental_data_filename in experimental_data_filenames])
                                         for required_data_file_pattern in required_data_file_patterns])
                                    
        if all_patterns_matched:
            experimental_data_present = True
        else:
            exp_data_update = True

    # perform an update of the base data dir
    if exp_data_update:
        experimental_data_url = "https://s3.us-east-2.amazonaws.com/tfbssexperimentaldata/data.tar.gz"
        experimental_data_down_loc = os.path.join(script_dir,'data.tar.gz')
##        current_versions_file = os.path.join(experimental_data_dir, "experimental_data.current_versions.json")
        print "Downloading the most current experimental data"
        logging.info(" ".join(["Downloading most current experimental data."]))

        try:
            wget.download(experimental_data_url, out=experimental_data_down_loc)
            tar = tarfile.open(experimental_data_down_loc)
            tar.extractall(experimental_data_dir)
            experimental_data_present = True

        except:
##            logging.warning(" ".join(["Error in downloading experimental data.  Check your internet connection, and make sure the transcript id is of the format 'ENST00000378357' so the correct species data can be downloaded"]))
            logging.warning(" ".join(["Error in downloading experimental data.  Check your internet connection."]))
            experimental_data_present = False

    return experimental_data_present

        
def experimentaldata(target_species):
    """
    Retrieve the experimental data for non-human species by downloading it from the Amazon repository.
    Activates if the current installation Data directory does not include data for the target species.
    Example location: https://s3.us-east-2.amazonaws.com/tfbssexperimentaldata/acanthochromis_polyacanthus.tar.gz
    """

    experimental_data_dir = os.path.join(script_dir, 'data')
    experimental_data_species_dir = os.path.join(experimental_data_dir, target_species)
    experimental_data_species_dir_tar = os.path.join(experimental_data_dir, ".".join([target_species, "tar.gz"]))

    if not os.path.exists(experimental_data_species_dir): 	
        aws_server = "https://s3.us-east-2.amazonaws.com"
        experimental_data_species_url = "/".join([aws_server, "tfbssexperimentaldata", ".".join([target_species, "tar.gz"])])
        print experimental_data_species_url
        logging.info("Downloading most current experimental data for %s." %target_species)

        try:
            wget.download(experimental_data_species_url, out=experimental_data_species_dir_tar)
            tar = tarfile.open(experimental_data_species_dir_tar)
            tar.extractall(experimental_data_dir)

        except:
            logging.warning(" ".join(["Error in downloading experimental data.  Check your internet connection."]))


def experimentalDataUpdater_beta():
    """
    Update the experimental data by downloading it from the Amazon repository.
    Using a file which contains an dictionary of the most up to date exp. data filenames,
    activates if data directory does not exist, if the data directory does not contain the most recent files,
    or if it has been >= 60 days since last update.
    This version of the updater is perhaps too error prone for this stage of development.
    """
    
    ## download experimental data if not already present or if it is outdated
    current_version_url = "https://s3.us-east-2.amazonaws.com/tfbssexperimentaldata/experimental_data.current_versions.json"
    experimental_data_url = "https://s3.us-east-2.amazonaws.com/tfbssexperimentaldata/data.tar.gz"
    experimental_data_down_loc = os.path.join(script_dir,'data.tar.gz')
    experimental_data_dir = os.path.join(script_dir, 'data')
    current_versions_file = os.path.join(experimental_data_dir, "experimental_data.current_versions.json")
    update_required = False
    
    if not os.path.exists(experimental_data_dir):
        directory_creator(experimental_data_dir)
        update_required = True
    else:
        if os.path.exists(current_versions_file):
            # check if all current versions are in the data dir
            current_versions = load_json(current_versions_file)
            current_versions_filenames = current_versions.values()
            owned_versions_filenames = os.listdir(experimental_data_dir)
            missing_files = [x for x in current_versions_filenames if x not in owned_versions_filenames]
            if len(missing_files) > 0:
                update_required = True

            # check if 60 days has passed since last check of versions
            if 'last_checked' in current_versions:
                current_versions_last_checked = current_versions['last_checked']
                if (time.time() - current_versions_last_checked)/(3600*24) >= 60:
                    update_required = True
        else:
            update_required = True

    # download the most current experimental data
    if update_required:
        print "Downloading the most current experimental data"
        logging.info(" ".join(["Downloading most current experimental data."]))

        try:
            wget.download(current_version_url, out=experimental_data_dir)
            wget.download(experimental_data_url, out=experimental_data_down_loc)
            tar = tarfile.open(experimental_data_down_loc)
            tar.extractall(experimental_data_dir)

        except:
            logging.warning(" ".join(["Error in downloading experimental data.  Check your internet connection."]))

##        os.remove(experimental_data_down_loc)
        
    # update the current versions file last checked time with current time 
    if os.path.exists(current_versions_file):
        current_versions = load_json(current_versions_file)
        current_versions['last_checked'] = time.time()
        dump_json(current_versions_file, current_versions)


def species_specific_data(target_species, chromosome, species_specific_data_dir):
    """
    Many datasets are species-specific.  If the current target species has species-specific datasets, load them.
    Altered to allow for multiple versions of data.  The matching files are sorted and those occuring later in the list, are presumably later editions.
    """

    logging.info(" ".join(["Species-specific data: loading"]))
    # load GERP locations
    gerp_conservation_locations_dict = {}
    species_group = ""
    gerp_data_dir = os.path.join(species_specific_data_dir, "gerp_data")
    gerp_conservation_locations_dict_filenames = [os.path.join(gerp_data_dir, x) for x in os.listdir(gerp_data_dir) if "gerp_conservation.locations_dict" in x and target_species in x]
    if len(gerp_conservation_locations_dict_filenames) > 0:
        gerp_conservation_locations_dict_filenames.sort()
        gerp_conservation_locations_dict_filename = gerp_conservation_locations_dict_filenames[-1]
        gerp_conservation_locations_dict = load_msgpack(gerp_conservation_locations_dict_filename)
        species_group = gerp_conservation_locations_dict_filename.split(".")[3]
    
    # load GERP conservation weights
    gerp_conservation_weight_dict = {}
    gerp_conservation_weight_dict_filenames = [os.path.join(gerp_data_dir, x) for x in os.listdir(gerp_data_dir) if "gerp_conservation.weight_dict" in x and target_species in x]
    if len(gerp_conservation_weight_dict_filenames) > 0:
        gerp_conservation_weight_dict_filenames.sort()
        gerp_conservation_weight_dict_filename = gerp_conservation_weight_dict_filenames[-1]
        gerp_conservation_weight_dict = load_msgpack(gerp_conservation_weight_dict_filename)

    # load human CAGE locs occuring near promoters
    cage_dict = {}
    cage_data_dir = os.path.join(species_specific_data_dir, "cage_data")
    if os.path.exists(cage_data_dir):
        cage_dict_filename = os.path.join(cage_data_dir, ".".join(["homo_sapiens", "CAGE", "Chr"+chromosome.upper(), "2000", "promoters", "genomic_coords", "msg"]))
        if os.path.exists(cage_dict_filename):
            cage_dict = load_msgpack(cage_dict_filename)

    # load CAGE locs occuring near promoters of TFs
    TF_cage_dict = {}
    cage_data_dir = os.path.join(species_specific_data_dir, "cage_data")
    if os.path.exists(cage_data_dir):
        TF_cage_dict_filename = os.path.join(cage_data_dir, ".".join(["homo_sapiens", "jasparTFs", "CAGE", "complete_protein_dict", "json"]))
        if os.path.exists(TF_cage_dict_filename):
            TF_cage_dict = load_json(TF_cage_dict_filename)

    # load CAGE dist weights
    cage_dist_weights_dict = {}
    if os.path.exists(cage_data_dir):
        cage_dist_weights_dict_filenames = [os.path.join(cage_data_dir, x) for x in os.listdir(cage_data_dir) if "cage_dist_weights" in x and target_species in x]
        if len(cage_dist_weights_dict_filenames) > 0:
            cage_dist_weights_dict_filenames.sort()
            cage_dist_weights_dict_filename = cage_dist_weights_dict_filenames[-1]
            cage_dist_weights_dict = load_json(cage_dist_weights_dict_filename)

    # load CAGE correlations
    cage_correlations_dict = {}
    if os.path.exists(cage_data_dir):
        cage_correlations_dict_filenames = [os.path.join(cage_data_dir, x) for x in os.listdir(cage_data_dir) if "rekeyed_combined_cage_corr_dict" in x and target_species in x]
        if len(cage_correlations_dict_filenames) > 0:
            cage_correlations_dict_filenames.sort()
            cage_correlations_dict_filename = cage_correlations_dict_filenames[-1]
            cage_correlations_dict = load_msgpack(cage_correlations_dict_filename)

    # load CAGE correlation weights
    cage_corr_weights_dict = {}
    if os.path.exists(cage_data_dir):
        cage_corr_weights_dict_filenames = [os.path.join(cage_data_dir, x) for x in os.listdir(cage_data_dir) if "cage_corr_weights" in x and target_species in x]
        if len(cage_corr_weights_dict_filenames) > 0:
            cage_corr_weights_dict_filenames.sort()
            cage_corr_weights_dict_filename = cage_corr_weights_dict_filenames[-1]
            cage_corr_weights_dict = load_json(cage_corr_weights_dict_filename)
            cage_corr_weights_dict = {float(k):v for k,v in cage_corr_weights_dict.iteritems()}

    # load CAGE keys
    cage_keys_dict = {}
    if os.path.exists(cage_data_dir):
        cage_keys_dict_filenames = [os.path.join(cage_data_dir, x) for x in os.listdir(cage_data_dir) if "cage_ids_key_dict" in x and target_species in x]
        if len(cage_keys_dict_filenames) > 0:
            cage_keys_dict_filenames.sort()
            cage_keys_dict_filename = cage_keys_dict_filenames[-1]
            cage_keys_dict = load_json(cage_keys_dict_filename)
        

    # load JASPAR tfs to Ensembl transcript ids
    jasparTFs_transcripts_dict_filenames = [os.path.join(species_specific_data_dir, x) for x in os.listdir(species_specific_data_dir) if "jasparTFs.transcripts.single_protein" in x and target_species in x]
    if len(jasparTFs_transcripts_dict_filenames) > 0:
        jasparTFs_transcripts_dict_filenames.sort()
        jasparTFs_transcripts_dict_filename = jasparTFs_transcripts_dict_filenames[-1]
        jasparTFs_transcripts_dict = load_json(jasparTFs_transcripts_dict_filename)
    else:
        jasparTFs_transcripts_dict = {}

    # load CpG score weights
    cpg_obsexp_weights_dict_filenames = [os.path.join(species_specific_data_dir, x) for x in os.listdir(species_specific_data_dir) if "cpg_obsexp_weights" in x and target_species in x]
    if len(cpg_obsexp_weights_dict_filenames) > 0:
        cpg_obsexp_weights_dict_filenames.sort()
        cpg_obsexp_weights_dict_filename = cpg_obsexp_weights_dict_filenames[-1]
        cpg_obsexp_weights_dict = load_json(cpg_obsexp_weights_dict_filename)
        cpg_obsexp_weights_dict = {float(k):float(v) for k,v in cpg_obsexp_weights_dict.iteritems()}
        cpg_obsexp_weights_dict_keys = cpg_obsexp_weights_dict.keys()
        cpg_obsexp_weights_dict_keys.sort()
    else:
        cpg_obsexp_weights_dict = {}
        cpg_obsexp_weights_dict_keys = []

    # load GTEx variants
    gtex_variants = {}
    gtex_data_dir = os.path.join(species_specific_data_dir, "gtex_data")
    if os.path.exists(gtex_data_dir):
        gtex_chrom_dict_filename = os.path.join(gtex_data_dir, ".".join([target_species, "gtex_v7", "Chr"+chromosome.upper(), "min_unique", "eqtls", "grch38","msg"]))
        if os.path.exists(gtex_chrom_dict_filename):
            gtex_variants = load_msgpack(gtex_chrom_dict_filename)

    # load GTEx weights
    gtex_weights_dict = {}
    if os.path.exists(gtex_data_dir):
        gtex_weights_dict_filenames = [os.path.join(gtex_data_dir, x) for x in os.listdir(gtex_data_dir) if "gtex_weights" in x and target_species in x]
        if len(gtex_weights_dict_filenames) > 0:
            gtex_weights_dict_filenames.sort()
            gtex_weights_dict_filename = gtex_weights_dict_filenames[-1]
            gtex_weights_dict = load_json(gtex_weights_dict_filename)
            gtex_weights_dict = {float(k):float(v) for k,v in gtex_weights_dict.iteritems()}

    # load human meta clusters from GTRD project
    gtrd_metaclusters_dict = {}
    gtrd_data_dir = os.path.join(species_specific_data_dir, "gtrd_data")
    if os.path.exists(gtrd_data_dir):
        gtrd_metaclusters_chrom_dict_filename = os.path.join(gtrd_data_dir, ".".join([target_species, "metaclusters", "interval", "Chr"+chromosome.upper(), "clipped", "ordered", "tupled", "msg"]))
        if os.path.exists(gtrd_metaclusters_chrom_dict_filename):
            gtrd_metaclusters_dict = load_msgpack(gtrd_metaclusters_chrom_dict_filename)
            
    # load metacluster overlap weights
    metacluster_overlap_weights_dict = {}
    if os.path.exists(gtrd_data_dir):
        metacluster_overlap_weights_dict_filenames = [os.path.join(gtrd_data_dir, x) for x in os.listdir(gtrd_data_dir) if "metaclusters_overlap_weights_dict" in x and target_species in x]
        if len(metacluster_overlap_weights_dict_filenames) > 0:        
            metacluster_overlap_weights_dict_filenames.sort()
            metacluster_overlap_weights_dict_filename = metacluster_overlap_weights_dict_filenames[-1]
            metacluster_overlap_weights_dict = load_json(metacluster_overlap_weights_dict_filename)
            metacluster_overlap_weights_dict = {float(k):float(v) for k,v in metacluster_overlap_weights_dict.iteritems()}

    # load human ATAC-Seq from Encode project
    atac_seq_dict = {}
    atac_seq_data_dir = os.path.join(species_specific_data_dir, "atac_data")
    if os.path.exists(atac_seq_data_dir):
        atac_seq_chrom_dict_filename = os.path.join(atac_seq_data_dir, ".".join([target_species, "atac-seq", "Chr"+chromosome.upper(), "msg"]))
        if os.path.exists(atac_seq_chrom_dict_filename):
            atac_seq_dict = load_msgpack(atac_seq_chrom_dict_filename)

    # load ATAC-Seq dist weights
    atac_dist_weights_dict = {}
    if os.path.exists(atac_seq_data_dir):
        atac_dist_weights_dict_filenames = [os.path.join(atac_seq_data_dir, x) for x in os.listdir(atac_seq_data_dir) if "atac_dist_weights" in x and target_species in x]
        if len(atac_dist_weights_dict_filenames) > 0:
            atac_dist_weights_dict_filenames.sort()
            atac_dist_weights_dict_filename = atac_dist_weights_dict_filenames[-1]
            atac_dist_weights_dict = load_json(atac_dist_weights_dict_filename)

    return gerp_conservation_locations_dict, gerp_conservation_weight_dict, species_group, cage_dict, TF_cage_dict, cage_dist_weights_dict, cage_correlations_dict, cage_corr_weights_dict, cage_keys_dict, jasparTFs_transcripts_dict, atac_dist_weights_dict, metacluster_overlap_weights_dict, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, gtex_variants, gtex_weights_dict, gtrd_metaclusters_dict, atac_seq_dict


def overlap_range(x,y):
    """
    Identify an overlap between two lists of two numbers.
    """

    x.sort()
    y.sort()
    
    return range(max(x[0], y[0]), min(x[-1], y[-1])+1)

################################################################################
# PWM analysis #################################################################
################################################################################
def pwm_maker(strand, motif_length, tf_motif, bg_nuc_freq_dict, neg_bg_nuc_freq_dict):
    """
    Make a PWM from a nucleotide frequency table.
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


def PWM_scorer(seq, pwm, pwm_dict, pwm_type):
    """
    Generate score for current seq given a pwm.
    """

    # set relevant variables based on whether the pwm is mono or dinucleotide
##    if pwm_type == "mono":
##    motif_dist = len(seq)
##    span = 1

##    if pwm_type == "dinuc":
##        motif_dist = len(seq) - 1
##        span = 2
    
    seq_score = 0.0
    for i in range(0, len(seq)):
##        seq_score += pwm[pwm_dict[seq[i:i+1]][i]]
        seq_score += pwm[pwm_dict[seq[i:i+1]]][i]
##        nuc = seq[i:i+1]
##        row = pwm_dict[nuc]
##        score = pwm[row][i]
##        seq_score += score
        
##    # account for sequences which are non-standard code    
##    non_standard_dict = {'R':['A','G'],
##                         'Y':['C','T'],
##                         'S':['G','C'],
##                         'W':['A','T'],
##                         'K':['G','T'],
##                         'M':['A','C'],
##                         'B':['C','G','T'],
##                         'D':['A','G','T'],
##                         'H':['A','C','T'],
##                         'V':['A','C','G'],
##                         'B':['C','G','T']}
##
##    possible_seqs = []
##    non_standards_in_seq = [x for x in non_standard_dict.keys() if x in seq]
##    if len(non_standards_in_seq) > 0:
##        for non_standard_in_seq in non_standards_in_seq:
##            variant_chars = non_standard_dict[non_standard_in_seq]
##            for variant_char in variant_chars:
##                variant_seq = seq.replace(non_standard_in_seq, variant_char)
##                possible_seqs.append(variant_seq)
##    else:
##        possible_seqs.append(seq)
##    
##    # iterate through candidate sequence, and score each mono or dinucleotide
##    possible_seq_scores = []
##    for possible_seq in possible_seqs:
##        seq_score = 0.0
##        for i in range(0, motif_dist):
##            nuc = possible_seq[i:i+span]
##            row = pwm_dict[nuc]
##            score = pwm[row][i]
##            seq_score += score
##        possible_seq_scores.append(seq_score)
##
##    # use the best score achieved
##    possible_seq_scores.sort()
##    seq_score = possible_seq_scores[-1]
    
    return seq_score


def tfbs_finder(transcript_name, alignment, target_tfs_list, TFBS_matrix_dict, target_dir, pwm_score_threshold_dict, all_pwms_loglikelihood_dict, unaligned2aligned_index_dict, promoter_after_tss, pval):
    """
    1.Convert PFM to PWM for each TF in the Jaspar dictionary.
    2.Score all positions in the cleaned sequence
    3.If score is greater than or equal to precomputed threshold, then keep, otherwise set to zero
    4.Return dictionary of pwm scores for [species][tf_name][strand]
    """
    start_time = time.time()
    
    tfbss_found_dict_outfilename = os.path.join(target_dir, "TFBSs_found.all.json")

    # Determine if the analysis has been done already, load results if so
    if os.path.isfile(tfbss_found_dict_outfilename):
        logging.info(" ".join(["tfbss_found_dict already exists: loading"]))
        tfbss_found_dict = load_json(tfbss_found_dict_outfilename)
                                                    
    # If results don't already exist, time to party
    else:
        tfbss_found_dict = {}        
        align_chars = '-N .'
        mononuc_pwm_dict = {"A":0,"C":1,"G":2,"T":3}
        
        entry = alignment[0]
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

        # Use empirical data from whole human genome, but limited to within -2000/+200 bp of TSSs
        # http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0033204
        # http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0033204.s023
    ##            bg_nuc_freq_dict = {'A':0.247, 'C':0.251, 'G':0.254, 'T':0.248}
    ##            neg_bg_nuc_freq_dict = {'A':0.247, 'C':0.251, 'G':0.254, 'T':0.248}

        # https://arxiv.org/pdf/q-bio/0611041.pdf
        # empirical data from complete genome
        bg_nuc_freq_dict = {'A':0.292, 'C':0.207, 'G':0.207, 'T':0.292}
        neg_bg_nuc_freq_dict = {'A':0.292, 'C':0.207, 'G':0.207, 'T':0.292}
        
        # iterate through each tf_name and its motif
        for tf_name in target_tfs_list:
            if tf_name in TFBS_matrix_dict:
                tf_motif = TFBS_matrix_dict[tf_name]
                motif_length = len(tf_motif[0])
                
                if motif_length > 0:
                    # retrieve precomputed threshold and other information required for calculating the pvalue of the score
                    tf_pwm_score_threshold_dict = pwm_score_threshold_dict[tf_name]                    
                    pvals_scores_list = [[k,v] for k,v in tf_pwm_score_threshold_dict.iteritems()]
                    pvals_scores_list_sorted = sorted(pvals_scores_list, key=itemgetter(1))
                    scores_list_sorted = [x[1] for x in pvals_scores_list_sorted]
                    pvals_list = [x[0] for x in pvals_scores_list_sorted]
                    pvals_list.sort()

                    if pval in tf_pwm_score_threshold_dict:
                        tf_pwm_score_threshold = tf_pwm_score_threshold_dict[pval]
                    else:
                        if pval == 1:
                            tf_pwm_score_threshold = -100000
                        else:
                            pval = pvals_list[bisect_left(pvals_list, pval)]
                            tf_pwm_score_threshold = tf_pwm_score_threshold_dict[pval]

                    tfbss_found_dict[tf_name] = []
                    
                    # iterate through the forward and reverse strand sequences
                    for strand, seq in seq_dict.iteritems():
                        pwm = pwm_maker(strand, motif_length, tf_motif, bg_nuc_freq_dict, neg_bg_nuc_freq_dict)
                        
                        seq_length = len(seq)
                        # iterate through the nt sequence, extract a current frame based on the motif size and score
                        for i in range(0, seq_length - motif_length):
                            current_frame = seq[i:i+motif_length]
                            current_frame_score = PWM_scorer(current_frame, pwm, mononuc_pwm_dict, 'mono')

                            # keep results that are above the precomputed threshold
    ##                            if current_frame_score >= tf_pwm_score_threshold:
                            if current_frame_score >= tf_pwm_score_threshold:    
                                current_frame_score = round(current_frame_score, 2)
                                
                                pval_index = bisect_left(scores_list_sorted, current_frame_score)
                                if pval_index >= len(pvals_scores_list_sorted):
                                    pval_index = -1
    ##                                elif pval_index == 0:
    ##                                    pval_index == 0
                                else:
                                    pval_index -= 1

                                # account for pvalues which are larger than the current largest, so that they can be distinguished appropriately in the result table.
                                if pval_index < 0:
                                    current_frame_score_pvalue = ">" + str(pvals_scores_list_sorted[0][0])
                                else:
                                    current_frame_score_pvalue = str(pvals_scores_list_sorted[pval_index][0])
                                    
                                hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end = start_end_found_motif(i, strand, seq_length, promoter_after_tss, motif_length)

                                # identify position in alignment from start of found motif in unaligned sequence
                                aligned_position = unaligned2aligned_index_dict[species][hit_loc_start]

                                tfbss_found_dict[tf_name].append([tf_name, current_frame, strand, hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end, current_frame_score, current_frame_score_pvalue])
                                # add to results dictionary by tf_name
##                                if tf_name in tfbss_found_dict:
##                                    tfbss_found_dict[tf_name].append([tf_name, species, current_frame, strand, hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end, current_frame_score, current_frame_score_pvalue, aligned_position])
##                                    tfbss_found_dict[tf_name].append([tf_name, current_frame, strand, hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end, current_frame_score, current_frame_score_pvalue])
##
##                                else:
##                                    tfbss_found_dict[tf_name] = [[tf_name, species, current_frame, strand, hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end, current_frame_score, current_frame_score_pvalue, aligned_position]]
##                                    tfbss_found_dict[tf_name] = [[tf_name, current_frame, strand, hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end, current_frame_score, current_frame_score_pvalue]]
        # sort results for each tf_name according to the position in alignment
##        for tf_name, hits in tfbss_found_dict.iteritems():
            # ref-point
##            tfbss_found_dict[tf_name] = sorted(hits, key = itemgetter(10))
##            tfbss_found_dict[tf_name] = sorted(hits, key = itemgetter(3))
            
##        # save the results to file
##        dump_json(tfbss_found_dict_outfilename, tfbss_found_dict)

    end_time = time.time()
    logging.info(" ".join(["total time for tfbs_finder() for this transcript:", str(end_time - start_time), "seconds"]))

    return tfbss_found_dict
    
                            
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
    """
    Create a dictionary for mapping aligned positions to unaligned positions.
    """

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


def find_clusters(ens_gene_id, chr_start, chr_end, alignment, target_species, chromosome, tfbss_found_dict, cleaned_aligned_filename, converted_gerps_in_promoter, gerp_conservation_weight_dict, converted_cages, converted_metaclusters_in_promoter, converted_atac_seqs_in_promoter, converted_eqtls, gtex_weights_dict, transcript_id, cage_dict, TF_cage_dict, cage_dist_weights_dict, atac_dist_weights_dict, metacluster_overlap_weights_dict, cpg_list, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, jasparTFs_transcripts_dict, cage_keys_dict, cage_correlations_dict, cage_corr_weights_dict, gtex_variants, gene_len):
    """
    For each target species hit:
    Identify the highest score for each species within the locality threshold.
    Create combined affinity score from the target species hit and those best scores from each species.
    If two target species hits are within the locality threshold from one another, choose the hit which has the highest combined affinity score.
    """
    start_time = time.time()
        
    cluster_dict = {}
    
    for tf_name, hits in tfbss_found_dict.iteritems():
        
        if len(hits) > 0:
            cluster_dict[tf_name] = []
            tf_len = len(hits[0][1])
            target_cages, tf_cages = cage_correlations_summing_preparation(transcript_id, cage_dict, TF_cage_dict, tf_name, jasparTFs_transcripts_dict)        
            eqtl_occurrence_log_likelihood = eqtl_overlap_likelihood(converted_eqtls, chr_start, chr_end, tf_len, gene_len, gtex_variants, ens_gene_id)

            for hit in hits:
                # ref-point
                combined_affinity_score = 0
                target_species_hit = hit
                target_species_pwm_score = target_species_hit[7]
                species_weights_sum = 0
                cage_weights_sum = 0
                eqtls_weights_sum = 0
                atac_weights_sum = 0
                metacluster_weights_sum = 0
                corr_weight_sum = 0
                tf_len = len(hit[1])

                # datasets only available for homo sapiens
                if target_species == "homo_sapiens":
                    cage_weights_sum = cage_weights_summing(transcript_id, target_species_hit, cage_dist_weights_dict, converted_cages)
                    eqtls_weights_sum = eqtls_weights_summing(eqtl_occurrence_log_likelihood, ens_gene_id, target_species_hit, converted_eqtls, gtex_weights_dict, chr_start, chr_end, gtex_variants, tf_len, gene_len)
                    atac_weights_sum = atac_weights_summing(transcript_id, target_species_hit, atac_dist_weights_dict, converted_atac_seqs_in_promoter)
                    metacluster_weights_sum = metacluster_weights_summing(transcript_id, target_species_hit, metacluster_overlap_weights_dict, converted_metaclusters_in_promoter)
    ##                corr_weight_sum, cage_correlations_hit_tf_dict = cage_correlations_summing(target_species_hit, transcript_id, target_cages, tf_cages, jasparTFs_transcripts_dict, cage_keys_dict, cage_correlations_dict, cage_corr_weights_dict, cage_correlations_hit_tf_dict)
                    corr_weight_sum = cage_correlations_summing(target_species_hit, transcript_id, target_cages, tf_cages, jasparTFs_transcripts_dict, cage_keys_dict, cage_correlations_dict, cage_corr_weights_dict)

                species_weights_sum = gerp_weights_summing(target_species, transcript_id, chromosome, target_species_hit, gerp_conservation_weight_dict, converted_gerps_in_promoter)    
                cpg_weight = cpg_weights_summing(transcript_id, target_species_hit, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, cpg_list)

                # calculate the complete score (combined affinity)
                experimental_weights = [species_weights_sum, cage_weights_sum, eqtls_weights_sum, atac_weights_sum, metacluster_weights_sum, cpg_weight, corr_weight_sum]
                combined_affinity_score += sum(experimental_weights) + target_species_pwm_score
                combined_affinity_score = round(combined_affinity_score, 2)

                hit.append(combined_affinity_score)
            
                experimental_weights_rounded = [round(x, 2) for x in experimental_weights]
                hit += experimental_weights_rounded
                cluster_dict[tf_name].append(hit)
                
##                if tf_name in cluster_dict:
##                    cluster_dict[tf_name].append([hit])
##                else:
##                    cluster_dict[tf_name] = [hit]

##    for tf_name, hits in cluster_dict.iteritems():
##        cluster_dict[tf_name] = thits)

    total_time = time.time() - start_time
    logging.info(" ".join(["total time for find_clusters() for this transcript:", str(total_time), "seconds"]))
    
##    return cluster_dict, cage_correlations_hit_tf_dict
    return cluster_dict


##def alignment_info_content(aligned_filename):
##    """
##    DEPRECATED: Now using GERP scores for conservation analysis.
##    The BioPython alignInfo.information_content module is slow,
##    and uses generic nucleotide frequencies (C,G = 0.4; A,T = 0.1).
##    Identify information content (IC) of all locations in the alignment.
##    This dictionary of ICs can then be used for putative TFBSs.
##    *FIX: Use an existing cleaned alignment object instead of reading from file.
##    """
##
##    bg_nuc_freq_dict = {'A':0.292, 'C':0.207, 'G':0.207, 'T':0.292}
##    alignment = AlignIO.read(aligned_filename, "fasta")
##    target_species_row = alignment[0]
##    info_content_dict = {}
##    pseudo_count = 1
##
##    for i in range(0, len(target_species_row)):   
##        alignment_col = alignment[:,i]
##        char_counts = collections.Counter(alignment_col)
##        chars_count = float(sum(char_counts.values()))
##        info_content = sum([((char_count+pseudo_count)/chars_count) * math.log(((char_count+pseudo_count)/chars_count)/bg_nuc_freq_dict[char]) for char, char_count in char_counts.iteritems() if char in "ACGT"])
##        info_content_dict[i] = info_content
##    
##    return info_content_dict


##def alignment_summary(alignment):
##    """
##    DEPRECATED: Now using GERP scores for conservation analysis.
##    DEPRECATED: Now using alignment_info_content().
##    Uses Biopython's information_content module, which is very slow.
##    The current implementation replaces the previous one which called this
##    function for every putative TFBS, which can greatly outnumber the length of the
##    alignment if the p-value threshold is set low.  Current implementation scores
##    the whole alignment at the beginning and stores IC content to a dict.
##    Build an alignment object.  Generate an alignment summary.
##    Score for information content at each position in the alignment.
##    """
##
##    # build alignment for analysis of conservation via information content
##    if len(alignment) > 1:
##        msl_list = []
##        for entry in alignment:
##            record = SeqRecord(Seq(entry['seq'], alphabet = generic_dna), id = entry['species'], description = "")
##            msl_list.append(record)
##
##        # Generate an alignment summary.
##        msl = MultipleSeqAlignment(msl_list)
##        msl_summary = AlignInfo.SummaryInfo(msl)
##
##    else:
##        msl_summary = None
##    
##    # Score for information content at each position in the alignment. 
##    info_content_dict = {}
##    if msl_summary != None:
##        for i in range(0, len(alignment[0]['seq'])-1):
####            if i%100 ==0:
####                print i
##            info_content = msl_summary.information_content(i, i+1, chars_to_ignore = ['N'])
##            info_content_dict[i] = info_content        
##
##    return info_content_dict


##def conservation_information_content(target_species_hit, info_content_dict):
##    """
##    DEPRECATED: Now using GERP scores for conservation analysis.
##    Sum pre-calculated values, for the information content at each location
##    in the alignment, across for this location in the alignment.
##    """
##
##    if len(info_content_dict) != 0:
##        start = target_species_hit[4]
##        end = target_species_hit[5]
##
##        info_content = sum([info_content_dict[x] for x in range(start, end)])
##
##    else:
##        info_content = 0
##
##    return info_content
    

##def alignment_summary(alignment):
##    """
##    Build an alignment object.  Generate an alignment summary.
##    """
##
##    if len(alignment) > 1:
##        # build alignment for analysis of conservation via information content
##        msl_list = []
##        for entry in alignment:
##            record = SeqRecord(Seq(entry['seq'], alphabet = generic_dna), id = entry['species'], description = "")
##            msl_list.append(record)
##
##        msl = MultipleSeqAlignment(msl_list)
##        msl_summary = AlignInfo.SummaryInfo(msl)
##
##    else:
##        msl_summary = None
##
##    return msl_summary
##        
##
##def conservation_information_content(target_species_hit, msl_summary):
##    """
##    For the target hit, extract the slice of the alignment where it occurs and
##    calculate the information content.
##    """
##
####    if len(alignment) > 1:
####        # build alignment for analysis of conservation via information content
####        msl_list = []
####        for entry in alignment:
####            record = SeqRecord(Seq(entry['seq'], alphabet = generic_dna), id = entry['species'], description = "")
####            msl_list.append(record)
####
####        msl = MultipleSeqAlignment(msl_list)
####        msl_summary = AlignInfo.SummaryInfo(msl)
##        # extract location of this predicted binding site
##    ##    target_species_hit = cluster[0]
##    ##    target_species_hit = hit
##
##    if msl_summary != None:
##        start = target_species_hit[4]
##        end = target_species_hit[5]
##
##        print msl_summary
##        info_content = msl_summary.information_content(start, end, chars_to_ignore = ['N'])
##
##    else:
##        info_content = 0
##
##    return info_content
##    


def eqtl_overlap_likelihood(converted_eqtls, chr_start, chr_end, tf_len, gene_len, gtex_variants, ens_gene_id):
    """
    Likelihood of eQTL occurrence.
    """

    eqtl_occurrence_log_likelihood = 0
    
    if ens_gene_id in gtex_variants:
        eqtl_occurrence_log_likelihood_dict = {}
    ##    if len(converted_eqtls) > 0:        
        transcript_len  = float(chr_end - chr_start)

        # determine size of search space, and probability of observing an eQTL in this gene.
        # GTEx searches for variants which occur over the span of the gene + 1,000,000 nt upstream+downstream.
        eqtl_search_space = 2000000 + gene_len
        associated_gtx_eqtls = gtex_variants[ens_gene_id]
        variant_count = len(associated_gtx_eqtls) * 1.
        eqtl_occurrence_log_likelihood = -1 * math.log(((tf_len * variant_count)/(eqtl_search_space-tf_len)) * (transcript_len/gene_len), 2)

    return eqtl_occurrence_log_likelihood
    

def eqtls_weights_summing(eqtl_occurrence_log_likelihood, ens_gene_id, target_species_hit, converted_eqtls, gtex_weights_dict, chr_start, chr_end, gtex_variants, tf_len, gene_len):
    """
    Identify if any of the eQTLs associated with this gene overlap this predicted TFBS.
    Retrieve the log-likelihood scores for all of them.
    Fix.
    """
  
    eqtl_weights = []

    if len(converted_eqtls) > 0:        
##        transcript_len  = float(chr_end - chr_start)

        # Only needs to be calculated once, therefore break into new function (eqtl_overlap_likelihood) so that it is not performed for each of the TFBS predictions #
        # Requires the tf_len, which will need to be stored in a dict instead of drawing it from the current hit in order to avoid repetition #
        # determine size of search space, and probability of observing an eQTL in this gene.
        # GTEx searches for variants which occur over the span of the gene + 1,000,000 nt upstream+downstream.
##        eqtl_search_space = 2000000 + gene_len
##        associated_gtx_eqtls = gtex_variants[ens_gene_id]
##        variant_count = len(associated_gtx_eqtls) * 1.
##        eqtl_occurrence_log_likelihood = -1 * math.log(((tf_len * variant_count)/(eqtl_search_space-tf_len)) * (transcript_len/gene_len), 2)

        # determine the weight score for likelihood of this magnitude eQTL.
        motif_start = target_species_hit[5]
        motif_end = target_species_hit[6]

        for converted_eqtl in converted_eqtls:
            converted_eqtl_start = converted_eqtl[0]
            converted_eqtl_end = converted_eqtl[1]
            converted_eqtl_score_mag = abs(converted_eqtl[2])

            overlap = overlap_range([motif_start, motif_end], [converted_eqtl_start, converted_eqtl_end])
##            if motif_start<=converted_eqtl_start<=motif_end or motif_start<=converted_eqtl_end<=motif_end:

            if len(overlap) > 0:
                eqtl_weight = gtex_weights_dict[converted_eqtl_score_mag]
                eqtl_weights.append(eqtl_weight + eqtl_occurrence_log_likelihood)

    eqtl_weights_sum = sum(eqtl_weights)

    return eqtl_weights_sum


def cage_correlations_summing_preparation(transcript_id, cage_dict, TF_cage_dict, tf_name, jasparTFs_transcripts_dict):
    """
    Extract transcript relevant cages (target_cages) once for each TF under analysis.  This could be reduced,
    to just once for the whole analysis of the current transcript (better that the last version where it was once per hit).
    Extract TF relevant data just once per TF, instead of for every hit.
    Identify the CAGE ides which are associated with the promoter of the TF currently under analysis.
    The correlation between these and the target gene will be extracted and correlating log-weight will be summed.
    """

    target_cages = []
    tf_cages = []

    # current transcript (target) cages
    if transcript_id in cage_dict:
        target_cages = cage_dict[transcript_id]

    if tf_name in TF_cage_dict:
        tf_cages = TF_cage_dict[tf_name]
##    # JASPAR tfs are often hetero multimers
##    # therefore we should parse the individual proteins and identify transcripts for each
##    split_tf_names = clean_jaspar_names([tf_name])
##    tf_transcripts = []
##    for split_tf_name in split_tf_names:
##        if split_tf_name in jasparTFs_transcripts_dict:
##            tf_transcripts += jasparTFs_transcripts_dict[split_tf_name]
##
##    # for each JASPAR transcript, compile associated FANTOM CAGEs
##    for tf_transcript in tf_transcripts:
##        tf_cages += cage_dict[tf_transcript]
##
##    # CAGEs may be shared by multiple transcripts with TSSs in close proximity
##    tf_cages = list(set([x[0] for x in tf_cages]))

    return target_cages, tf_cages


##def cage_correlations_summing(target_species_hit, transcript_id, target_cages, tf_cages, jasparTFs_transcripts_dict, cage_keys_dict, cage_correlations_dict, cage_corr_weights_dict, cage_correlations_hit_tf_dict):
def cage_correlations_summing(target_species_hit, transcript_id, target_cages, tf_cages, jasparTFs_transcripts_dict, cage_keys_dict, cage_correlations_dict, cage_corr_weights_dict):
    """
    Extract correlation values between CAGEs associated with a predicted TFBS protein,
    and CAGEs associated with the current gene.
    """

    corr_weights_ls = []
    corr_weight_sum = 0  

    # cages for all transcripts of the predicted TFBS's proteins
    tf_name = target_species_hit[0]

    for tf_cage in tf_cages:
        if tf_cage in cage_keys_dict:
            tf_cage_key = cage_keys_dict[tf_cage]
            for target_cage_list in target_cages:
                target_cage = target_cage_list[0]
                target_cage_key = cage_keys_dict[target_cage]

                if tf_cage_key in cage_correlations_dict:
                    if target_cage_key in cage_correlations_dict[tf_cage_key]:
                        cage_correlation = cage_correlations_dict[tf_cage_key][target_cage_key]
                        cage_corr_weight = cage_corr_weights_dict[abs(cage_correlation)]
                        corr_weight_sum += cage_corr_weight
                        corr_weights_ls.append(cage_corr_weight)

        else:
            print tf_cage

    if len(corr_weights_ls) > 0:
        corr_weights_ls.sort()
        corr_weight_sum = corr_weights_ls[-1]

##    cage_correlations_hit_tf_dict[tf_name] = corr_weight_sum
    
##    return corr_weight_sum, cage_correlations_hit_tf_dict
    return corr_weight_sum
           

def cage_weights_summing(transcript_id, target_species_hit, cage_dist_weights_dict, converted_cages):
    """
    Retrieve a log-likelihood score for this from the pre-existing dictionary.
    """

    cage_weights = []

    # ref-point
##    if transcript_id in cage_dict:
##    motif_start = target_species_hit[5]
##    motif_end = target_species_hit[6]
##    motif_midpoint = (motif_end + motif_start)/2
    motif_midpoint = (target_species_hit[5] + target_species_hit[6])/2
    
    for converted_cage in converted_cages:
##        transcript_cage_start = converted_cage[0]
##        transcript_cage_end = converted_cage[1]
##        cage_midpoint = (transcript_cage_start + transcript_cage_end)/2
        cage_midpoint = (converted_cage[0] + converted_cage[1])/2
        motif_cage_dist = str(abs(motif_midpoint - cage_midpoint))
        
        if motif_cage_dist in cage_dist_weights_dict:
            cage_weight = cage_dist_weights_dict[motif_cage_dist]
            cage_weights.append(cage_weight)              

    cage_weights_sum = sum(cage_weights)

    return cage_weights_sum


def atac_weights_summing(transcript_id, target_species_hit, atac_dist_weights_dict, converted_atac_seqs_in_promoter):
    """
    Identify ATAC-Seq peaks which are near a putative TFBS.
    Retrieve a log-likelihood score for this from the pre-existing dictionary.
    """

    atac_weights = []
    
    motif_start = target_species_hit[5]
    motif_end = target_species_hit[6]
    motif_midpoint = (motif_end + motif_start)/2
    
    for converted_atac in converted_atac_seqs_in_promoter:            
        transcript_atac_start = converted_atac[0]
        transcript_atac_end = converted_atac[1]
        atac_midpoint = (transcript_atac_start + transcript_atac_end)/2

        motif_atac_dist = str(abs(motif_midpoint - atac_midpoint))
        
        if motif_atac_dist in atac_dist_weights_dict:
            atac_weight = atac_dist_weights_dict[motif_atac_dist]
            atac_weights.append(atac_weight)              

    atac_weights_sum = sum(atac_weights)

    return atac_weights_sum


##def metacluster_weights_summing(transcript_id, target_species_hit, metacluster_dist_weights_dict, converted_metaclusters_in_promoter):
##    """
##    Generate a log-likelihood score for a putative TFBS based on the distances
##    to the nearest metaclusters.
##    """
##
##
##    metacluster_weights = []
##
####    if transcript_id in metacluster_dict:
####        transcript_metaclusters = metacluster_dict[transcript_id]        
##    motif_start = target_species_hit[6]
##    motif_end = target_species_hit[7]
##    motif_midpoint = (motif_end + motif_start)/2
##    
##    for converted_metacluster in converted_metaclusters_in_promoter:            
##        transcript_metacluster_start = converted_metacluster[0]
##        transcript_metacluster_end = converted_metacluster[1]
##        metacluster_midpoint = (transcript_metacluster_start + transcript_metacluster_end)/2
##
##        motif_metacluster_dist = str(abs(motif_midpoint - metacluster_midpoint))
##        
##        if motif_metacluster_dist in metacluster_dist_weights_dict:
##            metacluster_weight = metacluster_dist_weights_dict[motif_metacluster_dist]
##            metacluster_weights.append(metacluster_weight)              
##
##    metacluster_weights_sum = sum(metacluster_weights)
##
##    return metacluster_weights_sum


def metacluster_weights_summing(transcript_id, target_species_hit, metacluster_overlap_weights_dict, converted_metaclusters_in_promoter):
    """
    Identify the number of metaclusters which overlap this putative TFBS.
    Retrieve a log-likelihood score for this from the pre-existing dictionary.
    """  

    num_ovelapping_metaclusters = 0
    metacluster_weights_sum = 0

    # ref-point
    motif_start = target_species_hit[5]
    motif_end = target_species_hit[6]

    for converted_metacluster in converted_metaclusters_in_promoter:
        transcript_metacluster_start = converted_metacluster[0]
        transcript_metacluster_end = converted_metacluster[1]

        # doesn't work clusters are bigger than motifs
##        if transcript_metacluster_start<=motif_start<=transcript_metacluster_end or transcript_metacluster_start<=motif_end<=transcript_metacluster_end:
        overlap = overlap_range([motif_start, motif_end], [transcript_metacluster_start, transcript_metacluster_end])
        if len(overlap)>0:
            num_ovelapping_metaclusters += 1

    if num_ovelapping_metaclusters in metacluster_overlap_weights_dict:
        metacluster_weights_sum = metacluster_overlap_weights_dict[num_ovelapping_metaclusters]
    else:
        print "metacluster overlap sum not in weight dict"
        logging.warning(" ".join(["metacluster overlap sum not in weight dict"]))

    return metacluster_weights_sum


def gerp_weights_summing(target_species, transcript_id, chromosome, target_species_hit, gerp_conservation_weight_dict, converted_gerps_in_promoter):
    """
    Identify the gerps which are near this predicted TFBS.
    Retrieve a log-likelihood score for this distance from the pre-existing dictionary.
    """  

    # ref-point
    motif_start = target_species_hit[5]
    motif_end = target_species_hit[6]
    tf_len = len(target_species_hit[1])
    gerp_weights_sum = 0
    dists = []
    for converted_gerp_in_promoter in converted_gerps_in_promoter:
        converted_gerp_in_promoter_start = converted_gerp_in_promoter[0]
        converted_gerp_in_promoter_end = converted_gerp_in_promoter[1]

        dist = distance_solve([motif_start, motif_end], [converted_gerp_in_promoter_start, converted_gerp_in_promoter_end])
        dists.append(dist)
        dists.sort()
        best_dist = dists[0]

        # should only overlap be used?
##        if best_dist == 0:
        if best_dist <=1000:
            if target_species == "homo_sapiens":
                if best_dist in gerp_conservation_weight_dict[chromosome][str(tf_len)]:
                    gerp_weights_sum = gerp_conservation_weight_dict[chromosome][str(tf_len)][best_dist]

            else:
                if best_dist in gerp_conservation_weight_dict:
                    gerp_weights_sum = gerp_conservation_weight_dict[best_dist]

    return gerp_weights_sum


def cpg_weights_summing(transcript_id, target_species_hit, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, cpg_list):
    """
    Retrieve a CpG weight score based on the CpG obs/exp of the midpoint of the
    current predicted TFBS.
    """

    if len(cpg_obsexp_weights_dict_keys)>0:
        # retrieve locations and CpG obs/exp score for the midpoint of this predicted TFBS
        # ref-point
        motif_start = target_species_hit[3]
        motif_end = target_species_hit[4]
        motif_midpoint = (motif_end + motif_start)/2
        cpg_obsexp = cpg_list[motif_midpoint][-1]

        # extract the weight for the obsexp which is just less than the current obsexp
        next_lesser_obsexp_index = bisect_left(cpg_obsexp_weights_dict_keys, cpg_obsexp)
        if next_lesser_obsexp_index != len(cpg_obsexp_weights_dict_keys):  
            next_lesser_obsexp = cpg_obsexp_weights_dict_keys[next_lesser_obsexp_index]
        else:
            next_lesser_obsexp = cpg_obsexp_weights_dict_keys[-1]
        
        cpg_weight = cpg_obsexp_weights_dict[next_lesser_obsexp]
        
        return cpg_weight

    else:
        return 0

def clean_jaspar_names(uncleaned_jaspar_ids):
    """
    Clean names of jaspar transcription factor names.
    MSX3 <- lost in humans.
    RHOX11 <- only present in 3 species.
    DUX <- mouse only gene.
    EWSR1 <- didn't end up in the Ensembl BioMart export.
    MIX-A <- jaspar says present in xenopus laevis, but not even showing
    in Ensembl as a gene for any species.
    """

    special_dict = {"EWSR1-FLI1" : ["EWSR1","FLI1"]}
    names_list = []

    # split the combined names
    for uncleaned_jaspar_id in uncleaned_jaspar_ids:
        uncleaned_jaspar_id = uncleaned_jaspar_id.upper()
        split_names = uncleaned_jaspar_id.split("::")
        for name in split_names:
            names_list.append(name)

    # replace variants
    for i, name in enumerate(names_list):
        names_list[i] = name.replace("(VAR.2)","").replace("(VAR.3)","")

    tmp_list = []
    for i, name in enumerate(names_list):
        if name in special_dict:
            tmp_list += special_dict[name]
        else:
            tmp_list.append(name)

    names_list = list(set(tmp_list))
    names_list.sort()

    return names_list
                              

def target_species_hits_table_writer(sorted_clusters_target_species_hits_list, target_dir, name):
    """
    Write results to table for only target species.
    """

    output_table_name = os.path.join(target_dir, "TFBSs_found" + name  + "csv")

    with open(output_table_name, 'wb') as output_table:
        writerUS=csv.writer(output_table) 
        writerUS.writerow(['binding prot.', 'motif', 'strand', 'start', 'end', 'TSS-relative start', 'TSS-relative end', 'frame score', 'p-value', 'combined\naffinity\nscore', 'species\nweights\nsum', 'cage\nweights\nsum', 'eqtls\nweights\nsum', 'atac\nweights\nsum', 'metacluster\nweights\nsum', 'cpg\nweight', 'corr.\nweight\nsum'])

        # for all clusters which have pass thresholds, write full cluster to .csv
        # ref-point
        for hit in sorted_clusters_target_species_hits_list:
            pval_str = hit[8]
            if ">" not in pval_str:
                hit[8] = "{0:.3e}".format(Decimal(hit[8]))
##                hit[8] = ">" + "{0:.3e}".format(Decimal(hit[8].strip(">")))
##            else:
##                hit[8] = "{0:.3e}".format(Decimal(hit[8]))
            writerUS.writerow([str(x) for x in hit])


def sort_target_species_hits(cluster_dict):
    """
    Sort target_species hits which are part of a cluster by combined affinity score.
    """
    sorted_clusters_target_species_hits_list = []

    for tf_name, hits in cluster_dict.iteritems():
        for hit in hits:
            sorted_clusters_target_species_hits_list.append(hit)
##            sorted_clusters_target_species_hits_list += hit
            
    # ref-point
    sorted_clusters_target_species_hits_list = sorted(sorted_clusters_target_species_hits_list, key=itemgetter(9), reverse = True)
    return sorted_clusters_target_species_hits_list


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
def retrieve_genome_aligned(target_species, chromosome, strand, promoter_start, promoter_end):
    """
    Takes as input target_species CCDS start position and size of promoter to be extracted.  Retrieves genome aligned,
    corresponding regions in all orthologs.
    ***Alignment no longer necessary as the use of pre-computed GERP scores means that no conservation calculation needs to be made.
    Additionally, the newest implementation of Ensembl REST alignment will not retrieve target species sequence,
    if there is no alignment with other species at that location.
    For example, a request for homo_sapiens alignment from chr1:1-10,000 will only return the locations where an alignment
    exists with the target species group (e.g. mammals_low).  This may only exist at chr1:2000-9500.
    """

##    # Retrieve alignment if alignment FASTA does not already exist  
##    query_type = "/alignment/block/region/"
##    pre_options = target_species + "/" + chromosome + ":" + str(promoter_start) + "-" + str(promoter_end) + ":" + str(strand)   
##
##    if coverage == "low":
##        coverage_str = "EPO_LOW_COVERAGE"
##    else:
##        coverage_str = "EPO"
##        
##    options = pre_options + "?method=" + coverage_str + ";compact=1;content-type=application/json;species_set_group=" + species_group
##    alignment_decoded = ensemblrest(query_type, options, 'json', "", log=True)
##
##    if 'error' not in alignment_decoded:            
##        # remove those entries which are computed ancestral species
##        alignment = [x for x in alignment_decoded[0]['alignments'] if 'Ancestor' not in x['seq_region']]
##    else:
##        logging.warning(" ".join([alignment_decoded['error']]))
##        if alignment_decoded['error'].lower() == "no alignment available for this region":
    query_type = "/sequence/region/"
    pre_options = target_species + "/" + chromosome + ":" + str(promoter_start) + "-" + str(promoter_end) + ":" + str(strand)
    options = pre_options + "?content-type=application/json"
    target_only_decoded = ensemblrest(query_type, options, 'json', "", log=True)
    if 'seq' in target_only_decoded:
        target_only_decoded['species'] = target_species
        alignment = [target_only_decoded]

    else:
        alignment = []

    return alignment


def fasta_writer(alignment, outfile):
    """
    Write ensembl JSON alignment to fasta file.
    """
    
    if not os.path.isfile(outfile) or (os.path.isfile(outfile) and os.path.getsize(outfile) == 0):
        with open(outfile, "w") as aligned_file:
            for entry in alignment:
                record = SeqRecord(Seq(entry['seq'], alphabet = IUPAC.ambiguous_dna), id = entry['species'], description = "")
                SeqIO.write(record, aligned_file, 'fasta')


def remove_non_ACGT(alignment):
    """
    Remove non alignment characters and ambiguous nucleotides.  should consider changing to replacing any non ACGT char to '-'.
    """

    # account for sequences which are non-standard code    
    non_standard_dict = {'R':['A','G'],
                         'Y':['C','T'],
                         'S':['G','C'],
                         'W':['A','T'],
                         'K':['G','T'],
                         'M':['A','C'],
                         'B':['C','G','T'],
                         'D':['A','G','T'],
                         'H':['A','C','T'],
                         'V':['A','C','G'],
                         'B':['C','G','T']}
    
    non_alignment_chars = " .N"
    for entry in alignment:
        for non_alignment_char in non_alignment_chars:
            entry['seq'] = entry['seq'].replace(non_alignment_char, '-')

        for multi_char, replacement_list in non_standard_dict.iteritems():
                entry['seq'] = entry['seq'].replace(non_alignment_char, replacement_list[0])

    return alignment


def remove_gap_only(alignment):
    """
    Find columns in the alignment where the entire column is '-',
        replace the '-' with 'P', then remove the '*'.
    """

    if len(alignment) > 0:
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


def remove_duplicate_species(alignment, target_species):
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


##def selective_alignment(alignment):
##    """
##    Remove sequences from the alignment if they have less then 75% of the nucleotides of the target_species sequence.
##    Work needed: identify scenarios where length of target sequence affects alignment and subsequent scoring.
##    """
##
##    target_species_entry = alignment[0]
##    target_species_seq_2nd_half = target_species_entry['seq'][len(target_species_entry['seq'])/2:]
##    target_species_seq_2nd_half = target_species_seq_2nd_half.replace("-","").replace("N","").replace(" ","").replace(".","")
##    target_species_seq_2nd_half_len = len(target_species_seq_2nd_half)
##
##    cleaned_alignment = []
##    if target_species_seq_2nd_half_len > 0:
##        for entry in alignment:
##            entry_seq_2nd_half = entry['seq'][len(entry['seq'])/2:]
##            entry_seq_2nd_half = entry_seq_2nd_half.replace("-","").replace("N","").replace(" ","").replace(".","")      
##            entry_seq_2nd_half_len = len(entry_seq_2nd_half)
##            if float(entry_seq_2nd_half_len)/target_species_seq_2nd_half_len >= 0.75:
##                cleaned_alignment.append(entry)
##
##    return cleaned_alignment


def load_genome_aligned(aligned_filename):    
    """
    Load previously retrieved alignment fasta file into dictionary.
    """
    
    with open(aligned_filename, 'r') as alignment_handle:
        alignment_list = list(SeqIO.parse(alignment_handle, 'fasta'))
    alignment = [{'seq': str(entry.seq), 'species':entry.id} for entry in alignment_list if '[' not in entry.id]
        
    return alignment


def alignment_tools(ensembl_aligned_filename, cleaned_aligned_filename, target_species, chromosome, strand, promoter_start, promoter_end):
    """
    Return cleaned alignment for further analysis.
    """

    # if cleaned alignment file doesn't exist, or the size is zero.
    if not os.path.isfile(cleaned_aligned_filename) or (os.path.isfile(cleaned_aligned_filename) and os.path.getsize(cleaned_aligned_filename) == 0):

        # If uncleaned Ensembl alignment file doesn't exist, or the size is zero: retrieve from Ensembl, write to file.
        if not os.path.isfile(ensembl_aligned_filename) or (os.path.isfile(ensembl_aligned_filename) and os.path.getsize(ensembl_aligned_filename) == 0):
            alignment = retrieve_genome_aligned(target_species, chromosome, strand, promoter_start, promoter_end)
            fasta_writer(alignment, ensembl_aligned_filename)   

        # If uncleaned Ensembl file exists and size is not zero: clean, write to cleaned filename.
        if os.path.isfile(ensembl_aligned_filename) and (os.path.isfile(ensembl_aligned_filename) and os.path.getsize(ensembl_aligned_filename) > 0):          
            alignment = load_genome_aligned(ensembl_aligned_filename)
            alignment = remove_non_ACGT(alignment)
            alignment = remove_duplicate_species(alignment, target_species)
##            # analysis is now based on conservation around individual hits, so removing sequences based on completeness is wasteful
##            alignment = selective_alignment(alignment)
            alignment = remove_gap_only(alignment)
            fasta_writer(alignment, cleaned_aligned_filename)

        # Uncleaned alignment file still doesn't exist (or size is zero): note in logfile.
        else:
            logging.warning(" ".join(["No ensembl alignment, or size is zero"]))
            alignment = []
            
    # Cleaned alignment file exists and size is not zero: load cleaned alignment.
    else:
        alignment = load_genome_aligned(cleaned_aligned_filename)

    return alignment


def test_transcript_id(decoded_json_description, transcript_id):
    """     
    Test if the dictionary for the target transcript id:
    Indicates it is a transcript.
    Doesn't contain errors/Is complete.
    """

    transcript_id_pass = False

    if 'error' not in decoded_json_description:
        if 'object_type' in decoded_json_description:
            if decoded_json_description['object_type'].lower() == 'transcript':
                transcript_id_pass = True
            else:
                logging.warning(" ".join([transcript_id, "This input does not appear to be a valid Ensembl transcript ID.  Ensembl REST defines it as:", decoded_json_description['object_type']]))
        else:
            logging.warning(" ".join([transcript_id, "This input does not appear to be a valid Ensembl transcript ID.  Please check it for errors."]))
    else:
        logging.warning(" ".join(["Ensembl REST responds with error:", decoded_json_description['error']]))

    return transcript_id_pass


def transfabulator(transcript, transcript_dict_filename):
    """
    Given a transcript ID, retrieve Ensembl descriptive data for that transcript.
    """

    retrieve_transcript_data = True
    # load transcript position data from json file if it already exists
    if os.path.isfile(transcript_dict_filename):
        if os.path.getsize(transcript_dict_filename) != 0:
            logging.info(" ".join(["transcript_dict already exists: loading"]))
            decoded_json_description = load_json(transcript_dict_filename)
            retrieve_transcript_data = False

    # retrieve transcript position data from json file if it does not exist
    # or contains no data.
    if retrieve_transcript_data:
        # Set parameters for retrieving Ensembl data via REST
        query_type = '/lookup/id/'
        options = '?feature=transcript;content-type=application/json'

        # populate 'transcript_dict' dictionary with sub-dictionaries.
        # key[transcript_id] = {chromosome, strand, start, end} for each ensembl transcript id
        decoded_json_description = ensemblrest(query_type, options, 'json', transcript, log=True)
        decoded_json_description = {k.lower():v for k,v in decoded_json_description.iteritems()}
        
    return decoded_json_description


def transcript_data_retrieve(decoded_json_description, transcript_dict_filename, promoter_before_tss, promoter_after_tss):
    """
    The retrieved transcript data is assumed complete, extract important data.
    Write json data to file.
    Based on user-defined values for target region (referenced to TSS),
    calculate genomic coordinates of target region.
    """

    # Extract position data
    chromosome = decoded_json_description['seq_region_name']
    chr_start = decoded_json_description['start']
    chr_end = decoded_json_description['end']
    strand = decoded_json_description['strand']
    ens_gene_id = decoded_json_description['parent']
    target_species = decoded_json_description['species']
    if 'display_name' in decoded_json_description:
        transcript_name = decoded_json_description['display_name']
    else:
        transcript_name = ""
    dump_json(transcript_dict_filename, decoded_json_description)

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

    return target_species, transcript_name, ens_gene_id, chromosome, tss, strand, promoter_start, promoter_end, chr_start, chr_end


def gene_data_retrieve(ens_gene_id):
    """
    Retrieve gene data for the parent gene of the target transcript.
    """

    # determine likelihood of overlapping an eQTL at all.
    # Set parameters for retrieving Ensembl data via REST
    
    query_type = '/lookup/id/'
    options = '?feature=transcript;content-type=application/json'

    # retrieve_gene_len
    decoded_json_description = ensemblrest(query_type, options, 'json', ens_gene_id, log=True)
    decoded_json_description = {k.lower():v for k,v in decoded_json_description.iteritems()}
    gene_start = decoded_json_description['start']
    gene_end = decoded_json_description['end']
    gene_len = gene_end - gene_start

    return gene_len

################################################################################
# Regulatory & Conservation Features ###########################################
################################################################################
def retrieve_regulatory(chromosome, strand, promoter_start, promoter_end, regulatory_decoded_filename, target_species):
    """
    Retrieve ensembl JSON data for regulatory features within the coordinates provided.
    """

    # determine if the regulatory data has already been retrieved, if so load, if not retrieve.
    if os.path.isfile(regulatory_decoded_filename):
        logging.info(" ".join(["regulatory_decoded already exists: loading"]))
        regulatory_decoded = load_json(regulatory_decoded_filename)
        
    else:
        query_type = "/overlap/region/"
        pre_options = target_species + "/" + chromosome + ":" + str(promoter_start) + "-" + str(promoter_end) + ":" + str(strand)
        options = pre_options + "?feature=regulatory;content-type=application/json"
        regulatory_decoded = ensemblrest(query_type, options, 'json', "", log=True)

        # rename the Ensembl regulatory elements so that they don't overtake the space available for names.
        for reg in regulatory_decoded:
            if "description" in reg:
                if reg["description"] == "Transcription factor binding site":
                    reg["description"] = "Pred. TFBS"
                if reg["description"] == "Open chromatin region":
                    reg["description"] = "Open chromatin"
                if reg["description"] == "Predicted promoter":
                    reg["description"] = "Pred. promoter"

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
        if "error" not in reg:
            reg_id = reg['id']
            reg_start = reg['start']
            reg_end = reg['end']
            description = reg['description']

            if strand == 1:
                #[promoter_start][reg_start][reg_end][promoter_end][chr_start][TSS>GENE--->][chr_end]
                converted_reg_start = (tss - reg_start) * -1
                converted_reg_end = (tss - reg_end) * -1
                if reg_start <= promoter_start:
                    converted_reg_start = (-1 * promoter_before_tss)
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


##def alignment_conservation(aligned_filename):
##    """
##    Identify basic conservation of DNA sequence in the alignment.
##    """
##
##    alignment = AlignIO.read(aligned_filename, "fasta")
##    target_species_row = alignment[0]
##    species_num = float(len(alignment))
##    conservation = []
##    for i in range(0, len(target_species_row)):
##        
##        target_species_nuc = target_species_row[i]
##        if target_species_nuc != "-":
##            alignment_col = alignment[:,i]
##            common_char = collections.Counter(alignment_col).most_common()[0][0]
##            char_count = collections.Counter(alignment_col).most_common()[0][1]
##            if common_char != '-':
##                col_conservation = char_count/species_num
##            else:
##                col_conservation = alignment_col.count(target_species_nuc)/species_num
##            conservation.append(col_conservation)
##
##    return conservation


def CpG(aligned_filename):
    """
    Score the CpG content of the target_species sequence over a 200 nt window.
    """

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
            rolling_island = cpg_list[i-100:i] + cpg_list[i:]
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


def cage_position_translate(transcript_id,tss,cage_dict,promoter_start,promoter_end,strand,promoter_before_tss,promoter_after_tss):
    """
    Convert the CAGE data genome positions to those which can be mapped into the final figure.
    """

    converted_cages = []
    
    if transcript_id in cage_dict:
        cages = cage_dict[transcript_id]
        
        for cage in cages:
            cage_strand = cage[2]
            cage_start = cage[3]
            cage_end = cage[4]
            cage_desc = cage[6]

            if cage_strand == "+":
                #[promoter_start][cage_start][cage_end][promoter_end][chr_start][TSS>GENE--->][chr_end]
                converted_cage_start = (tss - cage_start) * -1
                converted_cage_end = (tss - cage_end) * -1

            if cage_strand == "-":
                #[chr_start][<---GENE<TSS][chr_end][promoter_start][cage_start][cage_end][promoter_end]
                converted_cage_start = (tss - cage_start)
                converted_cage_end = (tss - cage_end)

            converted_cage = [converted_cage_start, converted_cage_end, cage_desc]
            converted_cages.append(converted_cage)

        converted_cages = sorted(converted_cages, key=itemgetter(2))

    return converted_cages


def gtex_position_translate(ens_gene_id,gtex_variants,tss,promoter_start,promoter_end,strand,promoter_before_tss,promoter_after_tss):
    """
    Convert the GTEx data genome positions to those which can be mapped into the final figure.
    Reduce to those that are within range of the promoter before/after tss.
    """

    converted_eqtls = []
    if ens_gene_id in gtex_variants:
        eqtls = gtex_variants[ens_gene_id]

        for eqtl in eqtls:
            if len(eqtl) == 2:
                loc = eqtl[0]
                eqtl_length = 1
                eqtl_effect = eqtl[1]
            else:
                loc = eqtl[0]
                eqtl_length = eqtl[1]
                eqtl_effect = eqtl[2]

##            overlap = overlap_range([promoter_start, promoter_end], [loc, loc+eqtl_length])
            
##            if len(overlap) >0:
            if promoter_start<=loc<=promoter_end or promoter_start<=loc+eqtl_length<=promoter_end:
                #[promoter_start][eqtl_start][eqtl_end][promoter_end][chr_start][TSS>GENE--->][chr_end]
                if strand == 1:
                    converted_eqtl_start = (tss - loc) * -1
                    converted_eqtl_end = (tss - loc + eqtl_length) * -1

                #[chr_start][<---GENE<TSS][chr_end][promoter_start][eqtl_start][eqtl_end][promoter_end]
                if strand == -1:
                    converted_eqtl_start = (tss - loc)
                    converted_eqtl_end = (tss - loc + eqtl_length)

                # save to final list
                converted_eqtl = [converted_eqtl_start, converted_eqtl_end, eqtl_effect]
                converted_eqtls.append(converted_eqtl)

    return converted_eqtls


def distance_solve(r1, r2):
     # sort the two ranges such that the range with smaller first element
     # is assigned to x and the bigger one is assigned to y
     x, y = sorted((r1, r2))

     #now if x[1] lies between x[0] and y[0](x[1] != y[0] but can be equal to x[0])
     #then the ranges are not overlapping and return the differnce of y[0] and x[1]
     #otherwise return 0 
     if x[1] < y[0]:
        return y[0] - x[1]
     return 0
    

def gerp_positions_translate(target_dir, gerp_conservation_locations_dict, chromosome, strand, promoter_start, promoter_end, tss):
    """
    Identify GERP constrained conservation locations which occur within the defined promoter region.
    Convert positions of GERP elements (json) to coordinates usable by the plot_promoter function.
    Requires supplying location relative to TSS (i.e. negative).
    For creating a left to right plot of the promoter, regardless of strand:
    Converted_reg_start is the leftmost regulatory position.
    Converted_reg_end is the rightmost regulatory position.
    """

    potential_gerps_in_promoter = []
    gerps_in_promoter = []

    # because a prediction can occur at the start/end of a defined promoter
    extended_range = 1000
    
    if chromosome in gerp_conservation_locations_dict:
##        left_most_index = bisect_left([x[0] for x in gerp_conservation_locations_dict[chromosome]], promoter_start - extended_range)
##        right_most_index = bisect_right([x[0]+x[1] for x in gerp_conservation_locations_dict[chromosome]], promoter_end + extended_range)
                
##        potential_gerps_in_promoter = gerp_conservation_locations_dict[chromosome][left_most_index-1:right_most_index+1]
        for potential_gerp_in_promoter in gerp_conservation_locations_dict[chromosome]:
##        for potential_gerp_in_promoter in potential_gerps_in_promoter:
##            overlap = overlap_range([promoter_start - extended_range, promoter_end + extended_range], [potential_gerp_in_promoter[0], potential_gerp_in_promoter[0]+potential_gerp_in_promoter[1]])
            if promoter_start - extended_range<=potential_gerp_in_promoter[0]<=promoter_end + extended_range or promoter_start - extended_range<=potential_gerp_in_promoter[0]+potential_gerp_in_promoter[1]<=promoter_end + extended_range:
##            if len(overlap) > 0:
                gerps_in_promoter.append(potential_gerp_in_promoter)
             

    # convert the positions of the in-promoter metaclusters to tss-relative
    converted_gerps_in_promoter = []
    for gerp_in_promoter in gerps_in_promoter:
        gerp_start = gerp_in_promoter[0]
        gerp_end = gerp_start + gerp_in_promoter[1]

        if strand == 1:
            converted_gerp_start = (tss - gerp_start) * -1
            converted_gerp_end = (tss - gerp_end) * -1
        if strand == -1:
            converted_gerp_start = (tss - gerp_start)
            converted_gerp_end = (tss - gerp_end)

        converted_gerp = [converted_gerp_start, converted_gerp_end]
        converted_gerps_in_promoter.append(converted_gerp)

    return converted_gerps_in_promoter


def gtrd_positions_translate(target_dir, gtrd_metaclusters_dict, chromosome, strand, promoter_start, promoter_end, tss):
    """
    Identify GTRD metaclusters which occur within the defined promoter region.
    Convert positions of metaclusters (json) to coordinates usable by the plot_promoter function.
    Requires supplying location relative to TSS (i.e. negative).
    For creating a left to right plot of the promoter, regardless of strand:
    Converted_reg_start is the leftmost regulatory position.
    Converted_reg_end is the rightmost regulatory position.
    """

    potential_metaclusters_in_promoter = []        
    promoter_start_millions = promoter_start/1000000
    promoter_end_millions = promoter_end/1000000

    # retrieve the metacluster peaks on which the chrom that the transcript is found
    # if the millions place is the same for each then the metaclusters come from a single
    # subdict entry
    if promoter_start_millions == promoter_end_millions:
        if promoter_start_millions in gtrd_metaclusters_dict:        
            potential_metaclusters_in_promoter += gtrd_metaclusters_dict[promoter_start_millions]

    # have to account for the possibility that this location spans a millions place
    # e.g. from 999,000 - 1,001,000
    else:
        if promoter_start_millions in gtrd_metaclusters_dict:
            potential_metaclusters_in_promoter += gtrd_metaclusters_dict[promoter_start_millions]

        if promoter_end_millions in gtrd_metaclusters_dict:
            potential_metaclusters_in_promoter += gtrd_metaclusters_dict[promoter_end_millions]

    # identify if the metacluster occurs within user-defined promoter region
    metaclusters_in_promoter = []    
    for potential_metacluster in potential_metaclusters_in_promoter:
        metacluster_start = potential_metacluster[0]
        metacluster_end = metacluster_start + potential_metacluster[1]
        metacluster_peak_count = potential_metacluster[2]
        
        overlap = overlap_range([promoter_start, promoter_end], [metacluster_start, metacluster_end])
        
        if len(overlap) > 0:
            metaclusters_in_promoter.append(potential_metacluster)                

    # convert the positions of the in-promoter metaclusters to tss-relative
    converted_metaclusters_in_promoter = []
##    gtrd_outfilename = os.path.join(target_dir, os.path.basename(target_dir) + '.gtrd.txt')
##    with open(gtrd_outfilename, 'w') as gtrd_outfile:
    for metacluster_in_promoter in metaclusters_in_promoter:
        metacluster_start = metacluster_in_promoter[0]
        metacluster_end = metacluster_start + metacluster_in_promoter[1]
        metacluster_peak_count = metacluster_in_promoter[2]

        if strand == 1:
            converted_metacluster_start = (tss - metacluster_start) * -1
            converted_metacluster_end = (tss - metacluster_end) * -1
        if strand == -1:
            converted_metacluster_start = (tss - metacluster_start)
            converted_metacluster_end = (tss - metacluster_end)

        converted_metacluster = [converted_metacluster_start, converted_metacluster_end, metacluster_peak_count]
        converted_metaclusters_in_promoter.append(converted_metacluster)

    return converted_metaclusters_in_promoter


def atac_pos_translate(atac_seq_dict, chromosome, strand, promoter_start, promoter_end, tss):
    """
    Identify merged ATAC-Seq peaks which occur within the defined promoter region.
    Convert positions of ATAC-Seq peaks (json) to coordinates usable by the plot_promoter function.
    Requires supplying location relative to TSS (i.e. negative).
    For creating a left to right plot of the promoter, regardless of strand:
    Converted_reg_start is the leftmost regulatory position.
    Converted_reg_end is the rightmost regulatory position.
    """

    potential_atac_seqs_in_promoter = []
    
##    chromosome = "chr" + chromosome.lower()
##    if chromosome in atac_seq_dict:
    promoter_start_millions = promoter_start/1000000
    promoter_end_millions = promoter_end/1000000

    # retrieve the ATAC-Seq peaks on which the chrom that the transcript is found
    # if the millions place is the same for each then the atac-seqs come from a single subdict entry
    if promoter_start_millions == promoter_end_millions:
        if promoter_start_millions in atac_seq_dict:
            potential_atac_seqs_in_promoter += atac_seq_dict[promoter_start_millions]

    # have to account for the possibility that this location spans a millions place
    # e.g. from 999,000 - 1,001,000
    else:
        if promoter_start_millions in atac_seq_dict:
            potential_atac_seqs_in_promoter += atac_seq_dict[promoter_start_millions]

        if promoter_end_millions in atac_seq_dict:
            potential_atac_seqs_in_promoter += atac_seq_dict[promoter_end_millions]
            
    # identify if the ATAC-Seq peak occurs within user-defined promoter region
    atac_seqs_in_promoter = []    
    for potential_atac_seq in potential_atac_seqs_in_promoter:
        
        atac_seq_start = potential_atac_seq[0]
        atac_seq_end = potential_atac_seq[0] + potential_atac_seq[1]
        atac_seq_score = potential_atac_seq[2]

        overlap = overlap_range([promoter_start, promoter_end], [atac_seq_start, atac_seq_end])
        
        if len(overlap) > 0:
            atac_seqs_in_promoter.append([atac_seq_start, atac_seq_end, atac_seq_score])            

    # convert the positions of the in-promoter atac_seqs to tss-relative
    converted_atac_seqs_in_promoter = []
    for atac_seq_in_promoter in atac_seqs_in_promoter:
        atac_seq_start = atac_seq_in_promoter[0]
        atac_seq_end = atac_seq_in_promoter[1]
        atac_seq_score = atac_seq_in_promoter[2]

        if strand == 1:
            converted_atac_seq_start = (tss - atac_seq_start) * -1
            converted_atac_seq_end = (tss - atac_seq_end) * -1
        if strand == -1:
            converted_atac_seq_start = (tss - atac_seq_start)
            converted_atac_seq_end = (tss - atac_seq_end)

        converted_atac_seq = [converted_atac_seq_start, converted_atac_seq_end, atac_seq_score]
        converted_atac_seqs_in_promoter.append(converted_atac_seq)

    return converted_atac_seqs_in_promoter
    

##def plot_promoter(target_species, transcript_id, species_group, alignment, alignment_len, promoter_before_tss, promoter_after_tss, transcript_name, top_x_greatest_hits_dict, target_dir, converted_reg_dict, converted_gerps_in_promoter, cpg_list, converted_cages, converted_metaclusters_in_promoter, converted_atac_seqs_in_promoter, converted_eqtls, cage_correlations_hit_tf_dict):
def plot_promoter(target_species, transcript_id, species_group, alignment, alignment_len, promoter_before_tss, promoter_after_tss, transcript_name, top_x_greatest_hits_dict, target_dir, converted_reg_dict, converted_gerps_in_promoter, cpg_list, converted_cages, converted_metaclusters_in_promoter, converted_atac_seqs_in_promoter, converted_eqtls):
    """
    Plot the predicted TFBSs, onto a 5000 nt promoter graph, which possess support above the current strand threshold.
    ['binding_prot', 'species', 'motif', 'strand', 'start', 'end', 'TSS-relative start', 'TSS-relative end', 'frame score', 'p-value', 'pos in align.', 'combined affinity score', 'support']
    """

    if target_species == "homo_sapiens":
        fig = plt.figure(figsize=(10, 6))
        ax1 = plt.subplot2grid((20,1),(0,0), rowspan = 6, colspan = 11)
        ax8 = plt.subplot2grid((20,1),(6,0), rowspan = 2, colspan = 11)
        ax2 = plt.subplot2grid((20,1),(8,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax3 = plt.subplot2grid((20,1),(10,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax4 = plt.subplot2grid((20,1),(12,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax5 = plt.subplot2grid((20,1),(14,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax6 = plt.subplot2grid((20,1),(16,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax7 = plt.subplot2grid((20,1),(18,0), sharex=ax1, rowspan = 2, colspan = 11)

        # Set format of the plot(s)
        # Hide x-ticks for all plots except the lowest
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax3.get_xticklabels(), visible=False)
        plt.setp(ax4.get_xticklabels(), visible=False)
        plt.setp(ax5.get_xticklabels(), visible=False)
        plt.setp(ax6.get_xticklabels(), visible=False)
        plt.setp(ax8.get_xticklabels(), visible=False)

        # plt + ax labels
        ax1.text(1.02,.5,'Predicted TFBSs', verticalalignment='center', transform=ax1.transAxes, rotation='vertical', fontsize=8)
        ax1.set_ylabel("Combined Affinity Score", fontsize = 8, labelpad = 0)
        ax1.text(1.005,0.99,'+ strand', verticalalignment='top', transform=ax1.transAxes, rotation='vertical', fontsize=6)
        ax1.text(1.005,.01,'- strand', verticalalignment='bottom', transform=ax1.transAxes, rotation='vertical', fontsize=6)
        ax2.text(1.01,.5,'GERP\nConserv.\n'+species_group, verticalalignment='center', transform=ax2.transAxes, rotation='vertical', fontsize=5)
        ax3.text(1.01,.5,'CpG\nObs/Exp', verticalalignment='center', transform=ax3.transAxes, rotation='vertical', fontsize=6)
        ax4.text(1.01,.5,'eQTLs', verticalalignment='center', transform=ax4.transAxes, rotation='vertical', fontsize=6)
        ax5.text(1.01,.5,'TFBS\nMeta\nClusters', verticalalignment='center', transform=ax5.transAxes, rotation='vertical', fontsize=6)
        ax6.text(1.01,.5,'ATAC-Seq', verticalalignment='center', transform=ax6.transAxes, rotation='vertical', fontsize=6)
        ax7.text(1.01,.5,'CAGE\nPeaks\n(TSSs)', verticalalignment='center', transform=ax7.transAxes, rotation='vertical', fontsize=6)
        ax8.text(1.01,.5,'TF\nExpress.\nCorr.', verticalalignment='center', transform=ax8.transAxes, rotation='vertical', fontsize=6)

    ### as of now the data for non-human species is limited to predicted TFBSs, conservation, and CpG
    else:
        fig = plt.figure(figsize=(10, 6))
        ax1 = plt.subplot2grid((10,1),(0,0), rowspan = 6, colspan = 11)
        ax2 = plt.subplot2grid((10,1),(6,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax3 = plt.subplot2grid((10,1),(8,0), sharex=ax1, rowspan = 2, colspan = 11)

        # Set format of the plot(s)
        # Hide x-ticks for all plots except the lowest
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)

        # plt + ax labels
        ax1.text(1.02,.5,'Predicted TFBSs', verticalalignment='center', transform=ax1.transAxes, rotation='vertical', fontsize=8)
        ax1.set_ylabel("Combined Affinity Score", fontsize = 8, labelpad = 0)
        ax1.text(1.005,0.99,'+ strand', verticalalignment='top', transform=ax1.transAxes, rotation='vertical', fontsize=6)
        ax1.text(1.005,.01,'- strand', verticalalignment='bottom', transform=ax1.transAxes, rotation='vertical', fontsize=6)
        ax2.text(1.01,.5,'GERP\nConserv.\n'+species_group, verticalalignment='center', transform=ax2.transAxes, rotation='vertical', fontsize=5)
        ax3.text(1.01,.5,'CpG\nObs/Exp', verticalalignment='center', transform=ax3.transAxes, rotation='vertical', fontsize=6)


    # plot title
    title_str = target_species+"\n"+" ".join([transcript_name, transcript_id])
    fig.text(0.065, 0.5, title_str, horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, rotation='vertical', fontsize=14)
    mpl.rcParams['axes.linewidth'] = 1.1
    plt.xlabel("Nucleotide position before TSS", labelpad=5)
    
    # Generate data for each of the greatest_hits and plot corresponding bar
    color_series=['#FFB300','#803E75','#FF6800','#A6BDD7','#C10020','#CEA262','#817066','#007D34','#F6768E','#00538A','#FF7A5C','#53377A','#FF8E00','#B32851','#F4C800','#7F180D','#93AA00','#593315','#F13A13','#232C16']
    color_dict = {'CTCF':'#FF0000', 'TBP':'#FF00FF'}
    y_range = []
    labels_used = []

    # Create a list sorted descending by number of supporting sequences, so that lower support hits that overlap can be seen.
    sorted_by_ca_list = []
    for TF, great_hits in top_x_greatest_hits_dict.iteritems():
        for great_hit in great_hits:
            sorted_by_ca_list.append(great_hit)

    # ref-point
##    sorted_by_ca_list = sorted(sorted_by_ca_list, key=itemgetter(12), reverse=True)
    sorted_by_ca_list = sorted(sorted_by_ca_list, key=itemgetter(9), reverse=True)
    ### AX1: Predicted TFBSs
    for sorted_great_hit in sorted_by_ca_list:
        tf_name = sorted_great_hit[0]

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
##        binding_site_start = sorted_great_hit[6]
##        binding_site_end = sorted_great_hit[7]
##        combined_affinity = sorted_great_hit[12]
##        binding_strand = int(sorted_great_hit[3])
        #'binding prot.', 'motif', 'strand', 'start', 'end', 'TSS-relative start', 'TSS-relative end', 'frame score', 'p-value', 'combined\naffinity\nscore', 'species\nweights\nsum', 'cage\nweights\nsum', 'eqtls\nweights\nsum', 'atac\nweights\nsum', 'metacluster\nweights\nsum', 'cpg\nweight', 'corr.\nweight\nsum']
        binding_site_start = sorted_great_hit[5]
        binding_site_end = sorted_great_hit[6]
        combined_affinity = sorted_great_hit[9]
        binding_strand = int(sorted_great_hit[2])

        TF_center_point = float(binding_site_start + binding_site_end)/2
        TF_width = abs(binding_site_start - binding_site_end)
        x_series.append(TF_center_point)                  
        y_series.append(combined_affinity * binding_strand)
        y_range.append(combined_affinity)
        ax1.bar(x_series, y_series, facecolor = picked_color, edgecolor = edge_color, linewidth=1, alpha=0.9, align = 'center', width = TF_width, label = lab)

    # Set y-axis height based on number of entries in alignment
    y_range.sort()
    tens_y = int(y_range[-1])/10 + 1

    # Ensembl regulatory information
    # All will be horizontally plotted in some shade of red
    if len(converted_reg_dict) > 0:
        alpha_gradient = 1.0
##        alpha_gradient_dict = {1:0, 2:0.5, 3:0.25, 4:0.225}
        alpha_gradient_dict = {1:0}
        for i in range(2,100):
            alpha_gradient_dict[i] = 1./i
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

    ax1.axhline(0, color = 'black', linewidth=0.5)
    
    ### AX2: Add GERP conservation bars 
    for converted_gerp_in_promoter in converted_gerps_in_promoter:
        converted_gerp_start = converted_gerp_in_promoter[0]
        converted_gerp_end = converted_gerp_in_promoter[1]
        alpha_gradient = 1
        gerp_height = 1
        
        gerp_x_series = []
        gerp_y_series = []
        gerp_midpoint = float(converted_gerp_start + converted_gerp_end)/2
        gerp_x_series.append(gerp_midpoint)
        gerp_y_series.append(gerp_height)

        gerp_width = abs(converted_gerp_start - converted_gerp_end)
        ax2.bar(gerp_x_series, gerp_y_series, facecolor='black', edgecolor='black', alpha=alpha_gradient, align = 'center', width=gerp_width)

    ax2.set_yticks([0, 1])
    plt.setp(ax2.get_yticklabels(), fontsize=6)
    ax2.set_ylim(0, 1)

    
    ### AX3: CpG plot
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

    ax3.bar(range(-1 * alignment_len + promoter_after_tss, promoter_after_tss), gpc, color = 'black')
    if top_obs2exp <1:
        top_obs2exp = 1
    ax3.set_ylim(0, top_obs2exp)
    ax3.set_yticks([0, 0.6, 1])
    ax3.set_yticklabels([0, 0.6, 1], va='center')
    plt.setp(ax3.get_yticklabels(), fontsize=6)
    ax3.axhline(0.6, color = 'black', alpha = 0.4)

    ### human-specific experimental data
    if target_species == "homo_sapiens":

        ### AX7: CAGE plot
        cage_height = 1
        cage_labels = []
        for converted_cage in converted_cages:
            converted_cage_start = converted_cage[0]
            converted_cage_end = converted_cage[1]
            description = converted_cage[2]
            if ".." in description:
                description = ""
            cage_x_series = []
            cage_y_series = []
            cage_center_point = float(converted_cage_start + converted_cage_end)/2
            cage_x_series.append(cage_center_point)
            cage_y_series.append(cage_height)
            
            cage_width = abs(converted_cage_start - converted_cage_end)
            ax7.bar(cage_x_series, cage_y_series, facecolor='black', edgecolor='black', align = 'center', width=cage_width, label=description)

            # add label for the CAGE peak
            if -1 * promoter_before_tss <= converted_cage_start <= promoter_after_tss + 1 or -1 * promoter_before_tss <= converted_cage_end <= promoter_after_tss + 1:
                plt.text(cage_center_point, cage_height, description, color="red", rotation = 270, fontsize=5, horizontalalignment='center', verticalalignment='top')

        ax7.axes.get_yaxis().set_visible(False)

        ### AX5: GTRD plot
        gtrd_height = 1
        for converted_metacluster_in_promoter in converted_metaclusters_in_promoter:
            converted_metacluster_start = converted_metacluster_in_promoter[0]
            converted_metacluster_end = converted_metacluster_in_promoter[1]
            metacluster_peak_count = converted_metacluster_in_promoter[2]
            alpha_gradient = 0.5 + (metacluster_peak_count/1220.0)/2
            
            gtrd_x_series = []
            gtrd_y_series = []
            gtrd_center_point = float(converted_metacluster_start + converted_metacluster_end)/2
            gtrd_x_series.append(gtrd_center_point)
            gtrd_y_series.append(gtrd_height)

            gtrd_width = abs(converted_metacluster_start - converted_metacluster_end)
            ax5.bar(gtrd_x_series, gtrd_y_series, facecolor='black', edgecolor='black', alpha=alpha_gradient, align = 'center', width=gtrd_width)
        ax5.axes.get_yaxis().set_visible(False)

        ### AX6: ATAC-Seq plot
        for converted_atac_seq_in_promoter in converted_atac_seqs_in_promoter:
            converted_atac_seq_start = converted_atac_seq_in_promoter[0]
            converted_atac_seq_end = converted_atac_seq_in_promoter[1]
            atac_seq_peak_score = converted_atac_seq_in_promoter[2]
            alpha_gradient = 0.5 + atac_seq_peak_score/93.234864
            
            gtrd_x_series = []
            gtrd_y_series = []
            gtrd_midpoint = float(converted_atac_seq_start + converted_atac_seq_end)/2
            gtrd_x_series.append(gtrd_midpoint)
            gtrd_y_series.append(gtrd_height)

            gtrd_width = abs(converted_atac_seq_start - converted_atac_seq_end)
            ax6.bar(gtrd_x_series, gtrd_y_series, facecolor='black', edgecolor='black', alpha=alpha_gradient, align = 'center', width=gtrd_width)
        ax6.axes.get_yaxis().set_visible(False)

        ### AX4: eQTLs plot
        colors = ["green", "red"]
        magnitudes = []
        for converted_eqtl in converted_eqtls:
            converted_eqtl_start, converted_eqtl_end, converted_eqtl_mag = converted_eqtl
            if -1 * promoter_before_tss <= converted_eqtl_start <= promoter_after_tss + 1 or -1 * promoter_before_tss <= converted_eqtl_end <= promoter_after_tss + 1:
                eqtl_midpoint = float(converted_eqtl_start + converted_eqtl_end)/2
                eqtl_width = abs(converted_eqtl_start - converted_eqtl_end)
                eqtl_x_series = []
                eqtl_y_series = []
                eqtl_x_series.append(eqtl_midpoint)
                eqtl_y_series.append(converted_eqtl_mag)
                magnitudes.append(converted_eqtl_mag)
                if converted_eqtl_mag > 0:
                    c = colors[0]
                else:
                    c = colors[1]
                ax4.bar(eqtl_x_series, eqtl_y_series, facecolor=c, edgecolor=c, align = 'center', width=eqtl_width)
    ##            # arrow does not format properly, perhaps due to size.  y value starts not at 0, and arrow wraps over itself.
    ##            ax4.arrow(eqtl_midpoint, 0, 0, converted_eqtl_mag, color=c, length_includes_head = True, lw=10, width=0.01)

        ax4_yticks = [-1,0,1]
        if len(magnitudes) > 0:
            magnitudes.sort()
            ax4_yticks = [math.floor(magnitudes[0]), 0, math.ceil(magnitudes[-1])]
        ax4.set_yticks(ax4_yticks)
        ax4.set_yticklabels(ax4_yticks, va='center')
        plt.setp(ax4.get_yticklabels(), fontsize=6)
        ax4.axhline(0.0, color = 'black', alpha = 0.4)

        ### AX8: cage_correlations
        # rebuild dict to have just the top correlation
        # ref-point
##        plot_tfs_corrs_colors = [(tf_name, cage_correlations_hit_tf_dict[tf_name], color_dict[tf_name]) if tf_name in cage_correlations_hit_tf_dict else (tf_name, 0, color_dict[tf_name]) for tf_name in top_x_greatest_hits_dict]
        plot_tfs_corrs_colors = [(tf_name, hits_list[0][16],  color_dict[tf_name]) for tf_name, hits_list in top_x_greatest_hits_dict.iteritems()]
                                 
        plot_tfs_corrs_colors_sorted = sorted(plot_tfs_corrs_colors, key=itemgetter(1), reverse=True)
        ax8.bar(range(0, len(plot_tfs_corrs_colors_sorted)), [x[1] for x in plot_tfs_corrs_colors_sorted], color=[x[2] for x in plot_tfs_corrs_colors_sorted], edgecolor = "none")
        ax8.set_ylim(0, plot_tfs_corrs_colors_sorted[0][1]+1)
        ax8.set_xlim(-1, len(top_x_greatest_hits_dict))
        ax8.set_yticks([0, math.ceil(plot_tfs_corrs_colors_sorted[0][1])+1])
        plt.setp(ax8.get_yticklabels(), fontsize=6)   


    ## set ticks
    # based on 100's
    if y_range[-1] <= 100:
        for falling_y_thresh in range(100, -1, -10):
            if y_range[-1] < falling_y_thresh:
                y_thresh = falling_y_thresh
        ax1.set_yticks(range(-1* y_thresh, y_thresh+1, 10))
    else:
        ax1.set_yticks(range(-1 * (((tens_y*10)/100)+1)*100, (((tens_y*10)/100)+2)*100, 100))

##    if y_range[-1] <= 50:
##        ax1.set_yticks(range(-50, 50+1, 10))
##    else:
##        ax1.set_yticks(range(-1 * (((tens_y*10)/100)+1)*100, (((tens_y*10)/100)+2)*100, 100))

    ylabs=ax1.get_yticks().tolist()
    ylabs=[abs(x) for x in ylabs]
    ax1.set_yticklabels(ylabs)
    plt.setp(ax1.get_yticklabels(), fontsize=8)

    # Misc    
    plt.xlim([-1 * promoter_before_tss, promoter_after_tss + 1])

    # legend
    num_cols = 6
    ax1.legend(bbox_to_anchor=[0., 1.1, 1.0, .102], loc='center', ncol=num_cols, prop={'size':8}, mode="expand", borderaxespad=0.)
                      
    # produce .svg figure
    plt.subplots_adjust(hspace=0.40)
    fig.savefig(os.path.join(target_dir, os.path.basename(target_dir) + '.Promoterhisto'  + '.svg'), facecolor='white', bbox_inches='tight')
    plt.clf()
    plt.close()
    
##    # variable x-ticks
##    dist = promoter_before_tss + promoter_after_tss
##    rough_interval = dist/10
##    power = int(np.log10(rough_interval))
##    xtick_jump = (rough_interval/(10**power)) * 10**power
##    ax3.set_xticks(range(-1 * promoter_before_tss, promoter_after_tss + 1, xtick_jump))




################################################################################
# Initiating Variables #########################################################
################################################################################

################################################################################
# Execution ####################################################################
################################################################################
signal.signal(signal.SIGINT, signal_handler)

def main():
    """
    All the things.
    """
    total_time_start = time.time()
    print("Executing tfbs_footprinter version %s." % __version__)

    # Create directory for results
    output_dir = os.path.join(curdir, "tfbs_results")
    directory_creator(output_dir)

    # begin timing and logging
    logging.basicConfig(filename=os.path.join(os.path.dirname(output_dir), 'TFBS_footprinter.log'), level=logging.INFO, format='%(asctime)s:    [%(levelname)s]    %(message)s')
    logging.info(" ".join(["***NEW SET OF ANALYSES HAS BEGUN***"]))

    if is_online():
        args_lists, exp_data_update = get_args()

        # if experimental data dir does not exist or user has requested an exp data update, then update.
        experimental_data_present = experimentalDataUpdater(exp_data_update)

        if experimental_data_present:
            if len(args_lists) > 0:                
                # analysis variables
                # non-species-specific
                # dictionary of thresholds for each TF
                ##    # updated version, which requires the presence of a current versions file
                ##    pwm_score_threshold_dict_filename = os.path.join(experimental_data_dir, current_versions["jaspar_thresholds"])
                pwm_score_threshold_dict_filename = os.path.join(script_dir, 'data/all_tfs_thresholds.jaspar_2018.1.json')        
                pwm_score_threshold_dicta = load_json(pwm_score_threshold_dict_filename)
                pwm_score_threshold_dict = {}
                for k,v in pwm_score_threshold_dicta.iteritems():
                    pwm_score_threshold_dict[k] = {float(kk):float(vv) for kk,vv in v.iteritems()}

                # load mono-nuc PFMs
                TFBS_matrix_filename = os.path.join(script_dir, 'data/pwms.json')
                TFBS_matrix_dict = load_json(TFBS_matrix_filename)
                TFBS_matrix_dict = {k.upper():v for k,v in TFBS_matrix_dict.iteritems()}

                # load JASPAR PWM score weights
                all_pwms_loglikelihood_dict_filename = os.path.join(script_dir, 'data/all_pwms_loglikelihood_dict.reduced.msg')
                all_pwms_loglikelihood_dict = load_msgpack(all_pwms_loglikelihood_dict_filename)

                last_target_species = None
                
            for i, args_list in enumerate(args_lists):
                args, transcript_ids_filename, transcript_id, target_tfs_filename, promoter_before_tss, promoter_after_tss, top_x_tfs_count, pval = args_list
                print "Ensembl transcript id:", transcript_id
                logging.info(" ".join(["***ANALYSIS OF A NEW TRANSCRIPT HAS BEGUN:", transcript_id]))
                logging.info(" ".join(["Arguments used in this run:", str(args_list)]))

                # target dir naming
                start_end = "("+"_".join([str(promoter_before_tss), str(promoter_after_tss)])+")"
                target_dir_name = "_".join([transcript_id+start_end, str(pval)])
                target_dir = os.path.join(output_dir, target_dir_name)

                # check if results have been created for this query ***check multiple files also (e.g. output table and figure)***.
                cluster_dict_filename = os.path.join(target_dir, "cluster_dict.msg")
                if not os.path.exists(cluster_dict_filename):

                    # identify if the target transcript id exists in Ensembl
                    transcript_dict_filename = os.path.join(target_dir, "transcript_dict.json")

                    decoded_json_description = transfabulator(transcript_id, transcript_dict_filename)
                    transcript_id_pass = test_transcript_id(decoded_json_description, transcript_id)

                    # parse target transcript id data from successful retrieval, and continue
                    if transcript_id_pass:
                        # create target output dir 
                        directory_creator(target_dir)
                        logging.info(" ".join(["Results will be output to:", target_dir]))
                        target_species, transcript_name, ens_gene_id, chromosome, tss, strand, promoter_start, promoter_end, chr_start, chr_end = transcript_data_retrieve(decoded_json_description, transcript_dict_filename, promoter_before_tss, promoter_after_tss)
                        gene_len = gene_data_retrieve(ens_gene_id)
                        
                        # species-specific
                        species_specific_data_dir = os.path.join(script_dir, 'data', target_species)
                        experimentaldata(target_species)
                        if target_species != last_target_species or chromosome != last_chromosome:
                            if os.path.exists(species_specific_data_dir):
                                gerp_conservation_locations_dict, gerp_conservation_weight_dict, species_group, cage_dict, TF_cage_dict, cage_dist_weights_dict, cage_correlations_dict, cage_corr_weights_dict, cage_keys_dict, jasparTFs_transcripts_dict, atac_dist_weights_dict, metacluster_overlap_weights_dict, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, gtex_variants, gtex_weights_dict, gtrd_metaclusters_dict, atac_seq_dict = species_specific_data(target_species, chromosome, species_specific_data_dir)
                            last_target_species = target_species
                            last_chromosome = chromosome
                        
                        # load target tfs
                        if target_tfs_filename == "" or target_tfs_filename == None:
                            target_tfs_filename = None
                            target_tfs_list = TFBS_matrix_dict.keys()

                        if target_tfs_filename != None:
                            target_tfs_list = parse_tf_ids(target_tfs_filename)
                            target_tfs_list = compare_tfs_list_jaspar(target_tfs_list, TFBS_matrix_dict)

                        # filenames for alignment and ensembl regulatory data
                        ensembl_aligned_filename = os.path.join(target_dir, "alignment_uncleaned.fasta")
                        cleaned_aligned_filename = os.path.join(target_dir, "alignment_cleaned.fasta")
                        alignment = alignment_tools(ensembl_aligned_filename, cleaned_aligned_filename, target_species, chromosome, strand, promoter_start, promoter_end)

                        # continue if there is an alignment from Ensembl, and after cleaning
                        if len(alignment) > 0:

                            target_species_row = alignment[0]
                            alignment_len = len(target_species_row['seq'].replace('-',''))

                            # retrieve regulatory
                            regulatory_decoded_filename = os.path.join(target_dir, "regulatory_decoded.json")
                            regulatory_decoded = retrieve_regulatory(chromosome, strand, promoter_start, promoter_end, regulatory_decoded_filename, target_species)
                            converted_reg_dict = reg_position_translate(tss,regulatory_decoded,promoter_start,promoter_end,strand,promoter_before_tss,promoter_after_tss)

                            # conservation
                            converted_gerps_in_promoter = gerp_positions_translate(target_dir, gerp_conservation_locations_dict, chromosome, strand, promoter_start, promoter_end, tss)

                            # identify information content of each column of the alignment
                            cpg_list = CpG(cleaned_aligned_filename)

                            # identify CAGEs in proximity to Ensembl TSS, convert for plotting
                            converted_cages = cage_position_translate(transcript_id,tss,cage_dict,promoter_start,promoter_end,strand,promoter_before_tss,promoter_after_tss)

                            # identify eQTLs in proximity to Ensembl TSS, convert for plotting
                            converted_eqtls = gtex_position_translate(ens_gene_id,gtex_variants,tss,promoter_start,promoter_end,strand,promoter_before_tss,promoter_after_tss)

                            # GTRD metaclusters
                            converted_metaclusters_in_promoter = gtrd_positions_translate(target_dir, gtrd_metaclusters_dict, chromosome, strand, promoter_start, promoter_end, tss)

                            # ATAC-seq data
                            converted_atac_seqs_in_promoter = atac_pos_translate(atac_seq_dict, chromosome, strand, promoter_start, promoter_end, tss)

                            # create index of aligned to unaligned positions
                            unaligned2aligned_index_dict = unaligned2aligned_indexes(cleaned_aligned_filename)

                            cluster_dict_json_filename = os.path.join(target_dir, "cluster_dict.json")
                            if not os.path.exists(cluster_dict_filename):
                                # score alignment for tfbss
                                tfbss_found_dict = tfbs_finder(transcript_name, alignment, target_tfs_list, TFBS_matrix_dict, target_dir, pwm_score_threshold_dict, all_pwms_loglikelihood_dict, unaligned2aligned_index_dict, promoter_after_tss, pval)
                                
                                # sort through scores, identify hits in target_species supported in other species
                                cluster_dict = find_clusters(ens_gene_id, chr_start, chr_end, alignment, target_species, chromosome, tfbss_found_dict, cleaned_aligned_filename, converted_gerps_in_promoter, gerp_conservation_weight_dict,  converted_cages, converted_metaclusters_in_promoter, converted_atac_seqs_in_promoter, converted_eqtls, gtex_weights_dict, transcript_id, cage_dict, TF_cage_dict, cage_dist_weights_dict, atac_dist_weights_dict, metacluster_overlap_weights_dict, cpg_list, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, jasparTFs_transcripts_dict, cage_keys_dict, cage_correlations_dict, cage_corr_weights_dict, gtex_variants, gene_len)                            
                                tfbss_found_dict.clear()
                                dump_json(cluster_dict_json_filename, cluster_dict)

                            else:
                                cluster_dict = load_json(cluster_dict_filename)

                            # sort the target_species hits supported by other species
                            sorted_clusters_target_species_hits_list = sort_target_species_hits(cluster_dict)
                            target_species_hits_table_writer(sorted_clusters_target_species_hits_list, target_dir, ".sortedclusters.")
                            
                            # extract the top x target_species hits supported by other species
                            top_x_greatest_hits_dict = top_x_greatest_hits(sorted_clusters_target_species_hits_list, top_x_tfs_count)

                            # plot the top x target_species hits
                            if len(top_x_greatest_hits_dict) > 0:
                                plot_promoter(target_species, transcript_id, species_group, alignment, alignment_len, promoter_before_tss, promoter_after_tss, transcript_name, top_x_greatest_hits_dict, target_dir, converted_reg_dict, converted_gerps_in_promoter, cpg_list, converted_cages, converted_metaclusters_in_promoter, converted_atac_seqs_in_promoter, converted_eqtls)

            total_time_end = time.time()
            logging.info(" ".join(["Total time for", str(len(args_lists)), "transcripts:", str(total_time_end - total_time_start), "seconds"]) + "\n\n")

    else:
        print "System does not appear to be connected to the internet.  Exiting TFBS_footprinter."
