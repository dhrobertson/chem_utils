"""
similarities.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Compute either the intra or inter similarities
"""
import argparse
import logging
import pprint
import sys
import os
from chem_tools import Molecules

parser = argparse.ArgumentParser(prog='similarities.py',
    description='compute the similarities either intra (one file) or inter (between two files).')
parser.add_argument('-f1', '--file1', dest='file1', help='file1 - primary file for intra or inter')
parser.add_argument('-f2', '--file2', dest='file2', help='file2 - the reference file to compute inter similarities with file 1')
parser.add_argument('-o', '--output', dest='out_file', help='file to write results -- otherwise stdout')
parser.add_argument('-v', '--verbose', default=False, action="store_true", help='verbose output')

args = parser.parse_args()

# create logger
logger = logging.getLogger()
logger.level = logging.INFO if args.verbose else logging.ERROR
stream_handler = logging.StreamHandler(sys.stderr)
logger.addHandler(stream_handler)

# process the out_file argument
file1 = args.file1 if args.file1 else ""
file2 = args.file2 if args.file2 else ""

molecules_1 = Molecules.Molecules()
molecules_1.addFromSmiFile(file1)

sim = list()
if file2 == "":
    sim = molecules_1.getIntraSimilarities()
else:
    molecules_2 = Molecules.Molecules()
    molecules_2.addFromSmiFile(file2)
    sim = molecules_1.getInterSimilarities(molecules_2)

if len(sim):
    print("mol1 mol2 similarity")
    for s in sim:
        print(s)
