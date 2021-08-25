"""
opentarget_json.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with opentarget data from downloaded opentarget json files

Command to download/sync files
    rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/21.04/output/etl/json/targets .

TODO: Need to define data directory and do the sync/download?

"""
import os
import json
import logging

# Global variables for connections so we don't create multiple
_data_cache = dict()

#
# define parquet file location
#
data_types = {
    'targets' : { 'dir' : '/Users/dhrobertson/data/opentargets/json/targets'}
}

# create logger
logger = logging.getLogger()
# return the rds connector for the database -- only create once per database

def _read_json_files(dir):
    data = list()
    if not os.path.isdir(dir):
        logger.info("Unable to locate json directory \"{}\"".format(dir))
        return data
    files = os.listdir(dir)
    logger.debug("Files from dir \"{}\": {}".format(dir, files))
    for file in files:
        if '.json' in file:
            logger.debug("  reading file \"{}\"".format(file))
            full_path = os.path.join(dir, file)
            with open(full_path) as fo:
                for row in fo.readlines():
                    json_data = json.loads(row)
                    data.append(json_data)
    logger.info("  read {} data elements from {} files".format(len(data), len(files)))

    return data

def _get_data_cache(data_type):
    global _data_cache
    global data_types

    # make sure inputted value is defined
    if data_type not in data_types:
        logger.error(" Input \"{}\" not in defined data_types: {}".format(data_type, data_types))
        return []

    # check to see if already in global memory cache -- if not read into memory
    if data_type not in _data_cache:
        _data_cache[data_type] = _read_json_files(data_types[data_type]['dir'])

    # return the data object (list) from memory cache
    return _data_cache[data_type]

def get_target_info(gene_symbols):
    """ get the details for a target -- input either single item or a list """

    logger.debug("parameter: {} type(parameter): {}".format(gene_symbols, type(gene_symbols)))
    if type(gene_symbols) == type(str()):
        gene_symbols = [gene_symbols]

    _targets_cache = _get_data_cache('targets')
    records = list()
    # search first by approvedSymbol
    found = list()
    for item in _targets_cache:
        # find first by approved symbol
        if item['approvedSymbol'] in gene_symbols:
            logger.debug('  found items approvedSymbol {} in gene_symbols {} ... keeping'.format(item['approvedSymbol'], gene_symbols))
            records.append(item)
            found.append(item['approvedSymbol'])

    # return if found all
    if len(found) == len(gene_symbols):
        return(records)

    logger.debug('  found {} items by approvedSymbol {} '.format(len(found), found))
    # find others
    not_found = list()
    for gene in gene_symbols:
        if gene not in found:
            not_found.append(gene)
    logger.debug(' Did not find {} targets by approvedSymbol -- looking by synonym {}'.format(len(not_found), not_found))

    for item in _targets_cache:
        # find first by approved symbol
        for symbol in item['symbolSynonyms']:
            if symbol in not_found:
                logger.debug('  found items symbolSynonyms {} in gene_symbols {} ... keeping {}'.format(symbol, gene_symbols, item['approvedSymbol']))
                if item['approvedSymbol'] not in found:
                    records.append(item)
                    found.append(symbol)
                else:
                    logger.error(' it appears {} is a synonym for {} which was already found in the list. not duplicating'.format(symbol, item['approvedSymbol']))
                    found.append(symbol)

    if len(found) != len(gene_symbols):
        logger.error("ERROR: Some gene_symbols not found?: returned {} targets for {} ids".format(
                len(records),len(gene_symbols)))

    return records


### key method to do the api request for chembl
### Reference: Reference: https://pharos-api.ncats.io/graphql
def send_api_request(request):
    """
        function to send the api request to Pharos API
        sends the query request and return json
        captures request error
    """
    # OpenTarget Target Info
    # POST https://api.genetics.opentargets.org/graphql
    encoded_data = json.dumps(request)
    logger.info('OPENTARGETS_API: request: {}'.format(request))

    try:
        response = requests.post(
            url="https://api.genetics.opentargets.org/graphql",
            headers={
                "Content-Type": "application/json; charset=utf-8",
            },
            data=encoded_data
        )
        logger.info('OPENTARGET_API: Response HTTP Status Code: {}'.format(response.status_code))

        # capture 400 error
        if response.status_code == 400:
            # check for an html response
            logger.debug('OPENTARGET_API: Response HTTP Text: {}'.format(response.text))
            if 'html' in response.text:
                soup = BeautifulSoup(response.text, "html.parser")
                error = soup.find('pre').text
                logger.error("Request Error: {} response error: {}".format(
                    response.status_code,
                    error
                ))
                return {}
            # assume json coded errors
            errors = json.loads(response.text)
            logger.debug('OPENTARGET_API: Response HTTP Text (json): {}'.format(errors))
            if 'syntaxError' in errors:
                logger.error("Request Syntax Error: {}".format(errors['syntaxError']))
                return {}
            logger.error('Unknown errors: {}'.format(errors))
            return {}
        # 200 correct response
        content = json.loads(response.content)
        return content

    except requests.exceptions.RequestException:
        logger.error('HTTP Request failed')
        return {}
