"""
smi2sdf.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Convert SDF File to Smiles
"""
import argparse
import logging
import pprint
import sys
import os
from chem_utils.objects import Molecules
from rdkit import rdBase
from rdkit import Chem

parser = argparse.ArgumentParser(prog='smi2sdf.py',
    description='convert smi file to sdf using RDkit.')
parser.add_argument('-f', '--file', dest='in_file', help='input smi file', required=True)
parser.add_argument('-o', '--output', dest='out_file', help='file to write sdf results -- otherwise stdout', required=True)
parser.add_argument('-id', '--idfield', dest='id_field', help='id field name in sdf to store molecule name/id')
parser.add_argument('-v', '--verbose', default=False, action="store_true", help='verbose output')

args = parser.parse_args()

# create logger
logger = logging.getLogger()
logger.level = logging.INFO if args.verbose else logging.ERROR
stream_handler = logging.StreamHandler(sys.stderr)
logger.addHandler(stream_handler)

logger.info("RDKit Version: {}".format(rdBase.rdkitVersion))

id_field = args.id_field if args.id_field else "ID"

molecules = Molecules.Molecules()
molecules.addFromSmiFile(args.in_file)

# write out as sdf
SDWriter = Chem.rdmolfiles.SDWriter(args.out_file)
SDWriter.SetProps([id_field])
for i in range(0, molecules.length()):
    mol = molecules.getMol(i)
    rd_mol = mol['mol']
    rd_mol.SetProp(id_field, mol['name'])
    SDWriter.write(mol['mol'])

logger.info("Wrote {} molecules to {} file".format(molecules.length(), args.out_file))
