"""
utils.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with chembl data directly through the chembl api
"""
from pathlib import Path
import requests
import requests_cache
import json
import os
import sys
import logging
import csv
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

def get_chembl_client():
    return new_client

def read_alias_map():
    # TODO Need to find a way to cache this so only read once per instantiated instance
    """ read in the alias map to be able to get the mapping of gene name to chembl_id """
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

def get_bioactivities_for_molecules(chembl_ids):
    """ get the activities for a molecule -- input either single item or a list """

    mol_list = chembl_ids
    if type(mol_list) != type(list()):
        mol_list = [chembl_ids]

    results = send_api_request("activity", { 'molecule_chembl_id__in' : ','.join(mol_list)})

    return results['activities']

def get_bioactivities_for_targets(chembl_ids):
    """ get the activities for a target -- input either single item or a list """

    tgt_list = chembl_ids
    if type(tgt_list) != type(list()):
        tgt_list = [chembl_ids]

    results = send_api_request("activity", { 'target_chembl_id__in' : ','.join(tgt_list)})

    return results['activities']

def get_assay_details(chembl_ids):
    """ get the details for an assay -- input either single item or a list """

    if type(chembl_ids) == type(list()):
        results = send_api_request("assay", { 'assay_chembl_id__in' : ','.join(chembl_ids)})
        if len(results['assays']) != len(chembl_ids):
            logger.error("ERROR: Some ids not assays?: returned {} assays for {} ids".format(
                len(results['assays']),len(chembl_ids)))
        return results['assays']

    # single molecule request
    return send_api_request("assay/" + chembl_ids, {})
    return results

def get_molecule_details(chembl_ids):
    """ get the details for a molecule -- input either single item or a list """

    # list of molecule_chembl_ids
    if type(chembl_ids) == type(list()):
        records = list()
        results = send_api_request("molecule", { 'molecule_chembl_id__in' : ','.join(chembl_ids)})
        if len(results['molecules']) != len(chembl_ids):
            logger.error("ERROR: Some ids not molecules?: returned {} molecules for {} ids".format(
                len(results['molecules']),len(chembl_ids)))
        return results['molecules']

    # single molecule request
    return send_api_request("molecule/" + chembl_ids, {})
    return results


def get_target_details(chembl_ids):
    """ get the details for a target -- input either single item or a list """

    records = list()
    if type(chembl_ids) == type(list()):
        results = send_api_request("target", { 'target_chembl_id__in' : ','.join(chembl_ids)})
        if len(results['targets']) != len(chembl_ids):
            logger.error("ERROR: Some ids not targets?: returned {} targets for {} ids".format(
                len(results['targets']),len(chembl_ids)))
        records = results['targets']
    else:
        results = send_api_request("target/" + chembl_ids, {})
        records.append(results)

    if not len(records):
        return []

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

### key method to do the api request for chembl
### Reference: https://www.ebi.ac.uk/chembl/api/data/docs

def send_api_request(type, params, cache="chembl_api"):
    """
        function to send the api request to ChEMBL API
        sends the query request and return json
        captures request error
        reference: https://www.ebi.ac.uk/chembl/api/data/docs
    """
    _DEBUG_ = 0
    url = "https://www.ebi.ac.uk/chembl/api/data/" + type
    params["format"] = "json"
    logger.info('Sending API request to {url}'.format(url=url))

    try:
        response = requests.get(
            url=url,
            params=params
        )
        logger.info('Response HTTP Status Code: {status_code}'.format(
            status_code=response.status_code))
        logger.debug( 'Response HTTP Response Body: {content}'.format(
            content=response.text))
        # capture 400 errors
        if response.status_code == 404:
            logger.error("Request Error: Status 404 - page not found")
            return {}

        if response.status_code == 400:
            # check for an html response
            if 'html' in response.text:
                soup = BeautifulSoup(response.text, "html.parser")
                error = soup.find('pre').text
                logger.error("Request Error: {} response error: {}".format(
                    response.status_code,
                    error
                ))
                return {}
            # assume json coded errors
            error_codes = json.loads(response.text)
            for error in error_codes['errors']:
                logger.error("Request Error: {} {} {}".format(error['extensions']['code'], error['locations'], error['message']))
            return {}

        # 200 correct response
        content = json.loads(response.text)
        return content

    except requests.exceptions.RequestException:
        logger.error('HTTP Request failed')
        return {}
