"""
get_compounds.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with chembl data
"""
import argparse
import logging
import sys
import os
from rdkit import Chem
from chembl import core_utils
from chembl import compounds

parser = argparse.ArgumentParser(prog='get_compounds_by_target.py',
    description='download the molecules from chembl and put data into appropriate csv format.')
parser.add_argument('-f', '--file', dest='in_file', help='txt file containing list of target_chembl_id')
parser.add_argument('-t', '--target', dest='target', help='a single target -- target_chembl_id')
parser.add_argument('-o', '--output', dest='out_file', help='file to write results -- otherwise stdout')
parser.add_argument('-v', '--verbose', default=False, action="store_true", help='verbose output')

args = parser.parse_args()

# create logger
logger = logging.getLogger()
logger.level = logging.INFO if args.verbose else logging.ERROR
stream_handler = logging.StreamHandler(sys.stderr)
logger.addHandler(stream_handler)

# process the out_file argument
out_file = args.out_file if args.out_file else ""

# create the target list
tgt_list = list()
if (args.in_file):
    tgt_list = core_utils.read_list_from_file(args.in_file)
if (args.target):
    tgt_list.append(args.target)

if not len(tgt_list):
    print('stopping. do not have any targets defined. Use either -f or -t must be defined')
    exit()

cpd_list = compounds.get_by_targets(tgt_list)
logger.info("found " + str(len(cpd_list)) + " compounds")

content = '\n'.join(cpd_list)
if out_file:
    content = content + '\n'
    core_utils.save_text_to_file(out_file, content)
else:
    print(content)
