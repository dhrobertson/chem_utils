"""
sdf2smi.py
------------
Author: Daniel H Robertson
https://iwatobipen.wordpress.com/2018/01/06/simple-way-for-making-smiles-file-rdkit/

Part of the IBRI cheminformatics system

Convert SDF File to Smiles
"""
import argparse
import logging
import pprint
import sys
import os
from rdkit import rdBase
from rdkit import Chem

parser = argparse.ArgumentParser(prog='sdf2smi.py',
    description='convert sdf file to smiles using RDkit.')
parser.add_argument('-f', '--file', dest='in_file', help='input sdf file')
parser.add_argument('-o', '--output', dest='out_file', help='file to write results -- otherwise stdout')
parser.add_argument('-v', '--verbose', default=False, action="store_true", help='verbose output')

args = parser.parse_args()

# create logger
logger = logging.getLogger()
logger.level = logging.INFO if args.verbose else logging.ERROR
stream_handler = logging.StreamHandler(sys.stderr)
logger.addHandler(stream_handler)

logger.info("RDKit Version: {}".format(rdBase.rdkitVersion))

# process the out_file argument
out_file = args.out_file if args.out_file else ""

# create the target list
sdf = Chem.SDMolSupplier(args.in_file)
name = None
nskipped = 0
for mol in sdf:
    try:
        name = mol.GetProp('Catalog ID') # need to make this an option
        name = name.replace(' ', '_')
        smi = Chem.MolToSmiles(mol)
        print("{} {}".format(smi, name))
    except:
        nskipped += 1
        logger.error("Issue with structure of {} ... skipping.".format(name))

if nskipped:
    logger.error("Issue with {} structures that were skipped.".format(nskipped))
#with open('smiles.smi', 'w') as f:
#    for mol in sdf:
#        smi = Chem.MolToSmiles(mol)
#        f.write("{}\n".format(smi)
