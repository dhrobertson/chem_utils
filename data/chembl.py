"""
chembl.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with chembl data directly through the chembl api
"""
from .api import api_chembl
from pathlib import Path
import json
import os
import sys
import logging
import csv
#
# define helper files
#
_g2c_map_file_ = Path(__file__).parent / 'chembl_data/gene2chembl_id.csv'

# create logger and set logger level
logger = logging.getLogger()

def read_alias_map():
    # TODO Need to find a way to cache this so only read once per instantiated instance
    """ read in the alias map to be able to get the mapping of gene name to chembl_id """
    alias_map = dict()
    if not os.path.isfile(_g2c_map_file_):
        logger.critical("Unable to locate target map_file \"{}\" ... exiting".format(_g2c_map_file_))
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
            logger.debug("Read " + str(n) + " lines from file \"" + str(_g2c_map_file_) + "\"" )
        file_ref.close()
    return alias_map

def get_chembl_ids_for_targets(genes):
    """ get the chembl_ids for a gene name or alias -- from the configuration file """

    gene_list = list()
    if type(genes) == type(list()):
        gene_list = genes
    else:
        gene_list.append(genes)

    chembl_list = list()

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


def find_target_ids_by_accession(accessions):
    """ search for targets by uniprot accession number """
    # https://www.ebi.ac.uk/chembl/api/data/target?target_components__accession=Q13936
    # recommend using hugo_api to do the search of gene_symbol or synonyms
    acc_list = accessions
    if type(acc_list) != type(list()):
        acc_list = [acc_list]

    results = list()
    for acc in acc_list:
        acc_results = api_chembl.send_api_request("target", { 'target_components__accession' : acc})
        if 'targets' in acc_results:
            for result in acc_results['targets']:
                results.append(result)

    return results

def get_bioactivities_counts_for_assays(chembl_ids):
    """ get the count of the activities for an assay -- input either single item or a list """

    if type(chembl_ids) != type(list()):
        chembl_ids = [chembl_ids]

    for id in chembl_ids:

        results = api_chembl.send_api_request("activity", { 'assay_chembl_id__in' : ','.join(mol_list), 'limit' : 250}, True)

    if 'activities' not in results:
        return []

    return results['activities']

def get_bioactivities_for_molecules(chembl_ids):
    """ get the activities for a molecule -- input either single item or a list """

    mol_list = chembl_ids
    if type(mol_list) != type(list()):
        mol_list = [chembl_ids]

    results = api_chembl.send_api_request("activity", { 'molecule_chembl_id__in' : ','.join(mol_list), 'limit' : 250})

    if 'activities' not in results:
        return []

    return results['activities']

def get_assays_for_targets(chembl_ids):
    """ get the assays for a target -- input either single item or a list """

    not_list = False
    if type(chembl_ids) != type(list()):
        not_list = True
        chembl_ids = [chembl_ids]

    api_results = api_chembl.send_api_request("assay", { 'target_chembl_id__in' : ','.join(chembl_ids), 'limit' : 250})
    if 'assays' not in api_results:
        return []

    mapping = dict()
    results = api_results['assays']
    for result in results:
        logger.debug(result)
        tgt = result['target_chembl_id']
        assay = result['assay_chembl_id']
        if tgt not in mapping:
            mapping[tgt] = list()
        mapping[tgt].append(assay)

    logger.debug(mapping)
    if not_list:
        return mapping[chembl_ids[0]]

    return mapping

def get_bioactivities_for_targets(chembl_ids):
    """ get the activities for a target -- input either single item or a list """

    tgt_list = chembl_ids
    if type(tgt_list) != type(list()):
        tgt_list = [chembl_ids]

    results = api_chembl.send_api_request("activity", { 'target_chembl_id__in' : ','.join(tgt_list), 'limit' : 250})
    if 'activities' not in results:
        return []

    return results['activities']

def get_assay_details(ids, chunk_size=100):
    """ get the details for an assay -- input either single item or a list """

    not_list = False
    chembl_ids = list()
    if type(ids) != type(list()):
        not_list = True
        chembl_ids.append(ids)
    else:
        for id in ids:
            chembl_ids.append(id)

    total_ids = len(chembl_ids)
    results = list()
    while len(chembl_ids):
        run_list = list()
        for i in range(min(len(chembl_ids), chunk_size)):
            run_list.append(chembl_ids.pop())
        logger.info("running {} remaining {}".format(len(run_list), len(chembl_ids)))
        api_results = api_chembl.send_api_request("assay", { 'assay_chembl_id__in' : ','.join(run_list)})
        for result in api_results['assays']:
            results.append(result)

    if len(results) != total_ids:
        logger.error("ERROR: Some ids not assays?: returned {} assays for {} ids".format(
            len(results), total_ids))

    # add in the number of results
    for i in range(0, len(results)):
        assay_chembl_id = results[i]['assay_chembl_id']
        results[i]['results_count'] = api_chembl.send_api_request("activity", { 'assay_chembl_id__in' : assay_chembl_id, 'limit' : 10 }, return_count=True)

    if not_list:
        return results[0]

    return results

def get_molecule_details(chembl_ids):
    """ get the details for a molecule -- input either single item or a list """

    # list of molecule_chembl_ids
    if type(chembl_ids) == type(list()):
        records = list()
        results = api_chembl.send_api_request("molecule", { 'molecule_chembl_id__in' : ','.join(chembl_ids)})
        if len(results['molecules']) != len(chembl_ids):
            logger.error("ERROR: Some ids not molecules?: returned {} molecules for {} ids".format(
                len(results['molecules']),len(chembl_ids)))
        return results['molecules']

    # single molecule request
    return api_chembl.send_api_request("molecule/" + chembl_ids, {})


def get_target_details(chembl_ids):
    """ get the details for a target -- input either single item or a list """

    records = list()
    if type(chembl_ids) == type(list()):
        results = api_chembl.send_api_request("target", { 'target_chembl_id__in' : ','.join(chembl_ids)})
        if len(results['targets']) != len(chembl_ids):
            logger.error("ERROR: Some ids not targets?: returned {} targets for {} ids".format(
                len(results),len(chembl_ids)))
        records = results['targets']
    else:
        records = api_chembl.send_api_request("target/" + chembl_ids, {})

    if not len(records):
        return []

    # bring up the GENE_SYMBOL FROM the target_components
    # TODO: Add HUGO Gene Name
    new_records = list()
    for tgt in records:
        logger.debug("Processing tgt {}".format(tgt))
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

    return new_records
