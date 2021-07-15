"""
create_smi.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with chembl data
"""
import sys
import argparse
import logging
import pprint
from chembl import compounds
from chembl import core_utils

#TODO: add flag for verbose
#TODO: add flag to switch between getting results from molecules vs targets

parser = argparse.ArgumentParser(prog='create_smi.py',
    description='create sdf file from list of molecules from chembl.')
parser.add_argument('-f', '--file', dest='in_file', help='txt file containing list of molecule_chembl_id')
parser.add_argument('-m', '--molecule', dest='molecule', help='a single molecule -- molecule_chembl_id')
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

# create the molecule list
mol_list = list()
if (args.in_file):
    mol_list = core_utils.read_list_from_file(args.in_file)
if (args.molecule):
    mol_list.append(args.molecule)

if not len(mol_list):
    print('stopping. do not have any molecules defined. Use either -f or -m must be defined')
    exit()

# fetch and write the assay results
smi = compounds.get_smi(mol_list)
if out_file:
    core_utils.save_text_to_file(out_file, smi)
else:
    print(smi)
