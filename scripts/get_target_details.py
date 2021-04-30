"""
get_target_details.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with chembl data
"""
import sys
import argparse
import logging
import pprint
import json
from chembl import utils
from chembl import core_utils

#TODO: add flag for verbose
#TODO: add flag to switch between getting results from molecules vs targets

parser = argparse.ArgumentParser(prog='get_target_details.py',
    description='download the target details from chembl and deliver as json dictionary with key as TARGET_CHEMBL_ID.')
parser.add_argument('-i', '--id', dest='chembl_id', help='a single target chembl id')
parser.add_argument('-f', '--file', dest='in_file', help='txt file containing list of chembl_ids')
parser.add_argument('-o', '--output', dest='out_file', help='file to write csv or json results -- otherwise stdout')
parser.add_argument('-v', '--verbose', default=False, action="store_true", help='verbose output')

args = parser.parse_args()

# create logger
logger = logging.getLogger()
logger.level = logging.INFO if args.verbose else logging.ERROR
stream_handler = logging.StreamHandler(sys.stderr)
logger.addHandler(stream_handler)

# process the out_file argument
out_file = args.out_file if args.out_file else ""

# read in the ids -- single or from list
id_list = list()
if (args.in_file):
    id_list = core_utils.read_list_from_file(args.in_file)
if (args.chembl_id):
    id_list.append(args.chembl_id)

if not len(id_list):
    logger.error('stopping. Do not have any chembl_ids defined. Either -i or -f must be defined (or both)')
    exit()

query_results = utils.get_target_details(id_list)
if not len(query_results):
    logger.error("Stopping. Did not return any target_details for {} id(s)".format(len(id_list)))
    exit()

targets_json = dict()
for tgt in query_results:
    targets_json[tgt['target_chembl_id']] = tgt

if (out_file):
    core_utils.write_json(out_file, targets_json)
else:
    print(json.dumps(targets_json, indent=2))
