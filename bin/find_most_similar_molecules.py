"""
find_most_similar_molecules.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with chembl data
"""
import sys
import os
import argparse
import logging
import pprint
from chem_tools import compare_molecules

def read_smiles_file(in_file):
    """ utility function to read smiles from a file """
    smiles_list = list()
    if not os.path.isfile(in_file):
        logger.critical("Unable to locate defined input file \"" + in_file +"\" for content ... skipping")
        return None
    file_ref = open(in_file, "r")
    for row in file_ref:
        row = row.rstrip()
        parts = row.split(' ')
        smiles_list.append({'name':parts[1], 'smiles':parts[0]})
    logger.info('read {} smiles from file {}'.format(len(smiles_list), in_file))

    return smiles_list

#TODO: add flag for verbose
#TODO: add flag to switch between getting results from molecules vs targets

parser = argparse.ArgumentParser(prog='find_identical_molecules.py',
    description='create sdf file from list of molecules from chembl.')
parser.add_argument('-f', '--file', dest='in_file', help='input smiles to match with the reference smiles')
parser.add_argument('-r', '--reference', dest='ref_file', help='reference smiles file to match input file smiles')
parser.add_argument('-v', '--verbose', default=False, action="store_true", help='verbose output')

args = parser.parse_args()
# create logger
logger = logging.getLogger()
logger.level = logging.INFO if args.verbose else logging.ERROR
stream_handler = logging.StreamHandler(sys.stderr)
logger.addHandler(stream_handler)

if not args.in_file:
    print('stopping. Please us -f to define the smiles file')
    exit()
smiles = read_smiles_file(args.in_file)

if not args.ref_file:
    print('stopping. Please us -r to define the reference smiles file')
    exit()
smiles_ref = read_smiles_file(args.ref_file)

for s1 in smiles:
    sim = 0.0
    sim_mol = None
    for s2 in smiles_ref:
        s = compare_molecules.compute_similarity(s1['smiles'],s2['smiles'])
        #print(sim, sim_mol, s2['name'], s, s>sim)
        if (s > sim):
            sim = s
            sim_mol = s2['name']
    print('{}:{} {}'.format(s1['name'], sim_mol, sim))
