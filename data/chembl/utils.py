"""
utils.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with chembl data
"""
from pathlib import Path
import unittest
import os
import sys
import logging
import csv
import pprint

# Python modules used for API access...
from chembl_webresource_client.new_client import new_client
#
# define helper files
#
_g2c_map_file_ = Path(__file__).parent / 'data/gene2chembl_id.csv'

# create logger
# TODO -- fix later -- check if logger exists -- and than attach to it
logger = logging.getLogger()
logger.level = logging.INFO
stream_handler = logging.StreamHandler(sys.stderr)
logger.addHandler(stream_handler)

def read_alias_map():
    alias_map = dict()
    if not os.path.isfile(_g2c_map_file_):
        logger.critical("Unable to locate target map_file \"" + _g2c_map_file_ +"\" ... exiting")
        exit(0)
    else:
        #logger.info("Reading data from file \"" + g2c_map_file + "\"" )
        with open(_g2c_map_file_) as file_ref:
            n = 0
            csvreader = csv.reader(file_ref, delimiter=',')
            for row in csvreader:
                #print(row)
                if n:
                    alias_map[row[0]] = row[1]
                n += 1
            logger.info("Read " + str(n) + " lines from file \"" + str(_g2c_map_file_) + "\"" )
        file_ref.close()
    return alias_map

def get_chembl_ids_for_targets(genes):

    gene_list = list()
    if type(genes) == type(list()):
        gene_list = genes
    else:
        gene_list.append(genes)

    chembl_list = list()

    #print(pprint.pformat(gene_list))

    alias_map = read_alias_map()
    for gene in gene_list:
        chembl_id = ''
        if gene in alias_map:
            chembl_id = alias_map[gene]
        chembl_list.append({ 'gene' : gene, 'chembl_id' : chembl_id})

    if type(genes) == type(list()):
        return chembl_list
    else:
        return chembl_list[0]['chembl_id']

def get_bioactivities_for_molecules(chembl_ids):
    mol_list = list()
    if type(chembl_ids) == type(list()):
        mol_list = chembl_ids
    else:
        mol_list.append(chembl_ids)

    record_list = list()
    for mol in mol_list:
        records = new_client.activity.filter(molecule_chembl_id=mol)
        for record in records:
            record_list.append(record)

    return record_list

def get_bioactivities_for_targets(chembl_ids):
    tgt_list = list()
    if type(chembl_ids) == type(list()):
        tgt_list = chembl_ids
    else:
        tgt_list.append(chembl_ids)

    record_list = list()
    for tgt in tgt_list:
        records = new_client.activity.filter(target_chembl_id=tgt)
        for record in records:
            record_list.append(record)

    return record_list

def get_molecule_details(chembl_ids):
    mol_list = list()
    if type(chembl_ids) == type(list()):
        mol_list = chembl_ids
    else:
        mol_list.append(chembl_ids)

    records = new_client.molecule.get(mol_list)
    if len(records) != len(mol_list):
        logger.error("ERROR: Some ids not molecules?: returned {} molecules for {} ids".format(len(records),len(mol_list)))

    if type(chembl_ids) == type(list()):
        return records
    else:
        if len(records) == 1:
            return records[0]
        return {}

def get_assay_details(chembl_ids):
    assay_list = list()
    if type(chembl_ids) == type(list()):
        assay_list = chembl_ids
    else:
        assay_list.append(chembl_ids)

    records = new_client.assay.get(assay_list)

    if type(chembl_ids) == type(list()):
        return records
    else:
        return records[0]

def get_target_details(chembl_ids):
    tgt_list = list()
    if type(chembl_ids) == type(list()):
        tgt_list = chembl_ids
    else:
        tgt_list.append(chembl_ids)

    records = new_client.target.get(tgt_list)
    if not len(records):
        return None

    # bring up the GENE_SYMBOL FROM the target_components
    # TODO: Add HUGO Gene Name
    new_records = list()
    for tgt in records:
        symbols = list()
        if 'target_components' in tgt:
            #print('  2:', pprint.pformat(tgt['target_components']))
            for component in tgt['target_components']:
                #print('  2:', pprint.pformat(component))
                for syn in component['target_component_synonyms']:
                    #print(syn)
                    if syn['syn_type'] == 'GENE_SYMBOL':
                        symbols.append(syn['component_synonym'])
        tgt['gene_symbol'] = ':'.join(symbols)
        new_records.append(tgt)

    if type(chembl_ids) == type(list()):
        return new_records
    else:
        return new_records[0]
