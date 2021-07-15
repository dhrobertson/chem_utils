"""
chembl_api.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with chembl data directly through the chembl api
"""
from pathlib import Path
import requests
import requests_cache
import json
from bs4 import BeautifulSoup
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
        acc_results = send_api_request("target", { 'target_components__accession' : acc})
        if 'targets' in acc_results:
            for result in acc_results['targets']:
                results.append(result)

    return results

def get_bioactivities_counts_for_assays(chembl_ids):
    """ get the count of the activities for an assay -- input either single item or a list """

    if type(chembl_ids) != type(list()):
        chembl_ids = [chembl_ids]

    for id in chembl_ids:

        results = send_api_request("activity", { 'assay_chembl_id__in' : ','.join(mol_list), 'limit' : 250}, True)

    if 'activities' not in results:
        return []

    return results['activities']

def get_bioactivities_for_molecules(chembl_ids):
    """ get the activities for a molecule -- input either single item or a list """

    mol_list = chembl_ids
    if type(mol_list) != type(list()):
        mol_list = [chembl_ids]

    results = send_api_request("activity", { 'molecule_chembl_id__in' : ','.join(mol_list), 'limit' : 250})

    if 'activities' not in results:
        return []

    return results['activities']

def get_assays_for_targets(chembl_ids):
    """ get the assays for a target -- input either single item or a list """

    not_list = False
    if type(chembl_ids) != type(list()):
        not_list = True
        chembl_ids = [chembl_ids]

    api_results = send_api_request("assay", { 'target_chembl_id__in' : ','.join(chembl_ids), 'limit' : 250})
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

    results = send_api_request("activity", { 'target_chembl_id__in' : ','.join(tgt_list), 'limit' : 250})
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
        api_results = send_api_request("assay", { 'assay_chembl_id__in' : ','.join(run_list)})
        for result in api_results['assays']:
            results.append(result)

    if len(results) != total_ids:
        logger.error("ERROR: Some ids not assays?: returned {} assays for {} ids".format(
            len(results), total_ids))

    # add in the number of results
    for i in range(0, len(results)):
        assay_chembl_id = results[i]['assay_chembl_id']
        results[i]['results_count'] = send_api_request("activity", { 'assay_chembl_id__in' : assay_chembl_id, 'limit' : 10 }, return_count=True)

    if not_list:
        return results[0]

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


def get_target_details(chembl_ids):
    """ get the details for a target -- input either single item or a list """

    records = list()
    if type(chembl_ids) == type(list()):
        results = send_api_request("target", { 'target_chembl_id__in' : ','.join(chembl_ids)})
        if len(results['targets']) != len(chembl_ids):
            logger.error("ERROR: Some ids not targets?: returned {} targets for {} ids".format(
                len(results),len(chembl_ids)))
        records = results['targets']
    else:
        records = send_api_request("target/" + chembl_ids, {})

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


## key method to do the api request for chembl
### Reference: https://www.ebi.ac.uk/chembl/api/data/docs

def send_api_request(api_request, params, level=0, _results_list=None, return_count=False):
    """
        function to send the api request to ChEMBL API
        sends the query request and return json
        captures request error
        reference: https://www.ebi.ac.uk/chembl/api/data/docs
    """
    _host_ = "www.ebi.ac.uk"

    types_map = {
        'molecule' : 'molecules',
        'activity' : 'activities',
        'target'   : 'targets',
        'assay'    : 'assays'
    }

    logger.info('ChEMBL API: api_request: {} parameters: {}'.format(api_request, params))

    # need to check if request is from recursion
    u_type = ""
    url  = ""
    if "chembl/api/data" in api_request:
        url = "http://" + _host_ + api_request
        u_type = api_request.split('/')[-1].split('.json?')[0]
        logger.debug("u_type: {}".format(u_type))
    else:
        url = "https://www.ebi.ac.uk/chembl/api/data/" + api_request
        u_type = api_request.split('/')[0]

    params["format"] = "json"

    # if no _results_list defined for recursive calls, create a empty list
    if _results_list == None:
        _results_list = list()
    logger.debug('Sending API request to {} of type: {} level: {}'.format(url, u_type, level, len(_results_list)))

    # confirm that it is a type that is understood


    if u_type not in types_map:
        logger.error('Unknown type \"{}\". Must be one of {}'.format(u_type, list(types_map.keys())))
        return []

    # send and parse the request
    try:
        response = requests.get(
            url=url,
            params=params
        )
        logger.info('ChEMBL API: Response HTTP Status Code: {status_code}'.format(
            status_code=response.status_code))
        logger.debug( 'Response HTTP Response Body: {content}'.format(
            content=response.text))
        # capture 400 errors
        if response.status_code == 404:
            logger.error("Request Error: Status 404 - page not found")
            return []

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
        if response.status_code == 200:
            content = json.loads(response.text)
            logger.debug("u_type: {} map: {}".format(u_type, types_map[u_type]))

            # confirm page_meta present - this means a simple api call
            if 'page_meta' not in content:
                # make sure returned as list
                return [content]

            # add results to the growing list
            logger.debug("u_type: {} map: {}".format(u_type, types_map[u_type]))
            total_count = None
            if 'total_count' in content['page_meta']:
                total_count = content['page_meta']['total_count']
            #logger.debug("type(content_json): {} content_json: {}".format(type(content_json), content_json))

            if return_count:
                return total_count

            for item in content[types_map[u_type]]:
                _results_list.append(item)
            logger.info("ChEMBL API: Results thus far: level: {} count: {} of {}".format(level, len(_results_list), total_count))

            # if page_meta['next'] != null -- continue recursion
            if content['page_meta']['next']:
                total = content['page_meta']['total_count']
                level = level + 1
                logger.debug("ChEMBL API: Sending new request: {}".format(content['page_meta']['next']))
                return send_api_request(content['page_meta']['next'], {}, level=level, _results_list=_results_list)

            # end of recursion return compiled results
            return { 'page_meta':{}, types_map[u_type]: _results_list }

        logger.error("Unkown HTTP response code: {}".format(response.status_code))
        return {}

    except requests.exceptions.RequestException:
        logger.error('HTTP Request failed')
        return {}
