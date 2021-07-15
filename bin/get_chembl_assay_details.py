"""
get_chembl_assay_details.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with chembl data
"""
import sys
import argparse
import logging
import json
from chem_utils.data import chembl_api

parser = argparse.ArgumentParser(prog='get_chembl_assay_details.py',
    description='retrieve assay details for assay_chembl_id on commandline')
parser.add_argument('assay_chembl_ids', nargs="+", help='list of assay_chembl_id(s)')
parser.add_argument('-v', '--verbose', default=False, action="store_true", help='verbose output')
parser.add_argument('-d', '--debug', default=False, action="store_true", help='very verbose output')

args = parser.parse_args()

# create logger
logger = logging.getLogger()
logger.level = logging.ERROR
if args.verbose:
    logger.level = logging.INFO
if args.debug:
    logger.level = logging.DEBUG
stream_handler = logging.StreamHandler(sys.stderr)
logger.addHandler(stream_handler)

assay_list = chembl_api.get_assay_details(args.assay_chembl_ids)
if len(assay_list) == 1:
    print(json.dumps(assay_list[0], indent=2))
else:
    print(json.dumps(assay_list, indent=2))
