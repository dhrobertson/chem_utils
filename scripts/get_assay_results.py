"""
get_active_compounds.py
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
from data.chembl import get_assay_results
from utils import core_utils

#TODO: add flag for verbose
#TODO: add flag to switch between getting results from molecules vs targets

parser = argparse.ArgumentParser(prog='get_assay_results.py',
    description='download the assay results for set of molecules from chembl and put data into appropriate csv format.')
parser.add_argument('-i', '--id', dest='chembl_id', help='a single chembl id - default target_chembl_id change using -m flag')
parser.add_argument('-f', '--file', dest='in_file', help='txt file containing list of chembl_ids')
parser.add_argument('-m', '--m_id', default=False, action="store_true", help='interpret ids as molecule_chembl_id. Otherwise default to target_chembl_id')
parser.add_argument('-o', '--output', dest='out_file', help='file to write csv or json results -- otherwise stdout')
parser.add_argument('-j', '--json', default=False, action="store_true", help='write output as json rather than csv')
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

id_type = "target"
if args.m_id:
    id_type = "molecule"

if id_type == "molecule":
    assay_results = get_assay_results.by_molecules(id_list)
else:
    assay_results = get_assay_results.by_targets(id_list)

if not len(assay_results):
    logger.error("Stopping. Did not return any assay results for {} id(s) of type {}".format(len(id_list),id_type))
    exit()

logger.level = logging.INFO
# write out the assay results
logger.info("found {} assay_results for {} ids of type {}".format(len(assay_results), len(id_list), id_type))
if (not args.json):
    core_utils.write_assay_results(out_file, assay_results)
else:
    dataset = get_assay_results.convert_to_standard(assay_results)
    if (out_file):
        core_utils.write_json(out_file, dataset)
    else:
        print(json.dumps(dataset, indent=2))
