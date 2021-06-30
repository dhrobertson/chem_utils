"""
utils.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with pharos data from API
"""
from pathlib import Path
import requests
import requests_cache
import json
from bs4 import BeautifulSoup
import logging

# create logger
logger = logging.getLogger()

def create_single_query_symbol(gene_symbol):
    """ create the graphQL query for the API for a single target """

    logger.debug("parameter: {} type(parameter): {}".format(gene_symbol, type(gene_symbol)))
    query = {
        "query": "query targetDetails{target(q:{sym:\"" + gene_symbol + "\"}){name tdl fam sym description novelty}}",
        "variables": {}
    }
    return query

def get_target_info(gene_symbols):
    """ get the details for a target -- input either single item or a list """

    logger.debug("parameter: {} type(parameter): {}".format(gene_symbols, type(gene_symbols)))
    if type(gene_symbols) == type(str()):
        gene_symbols = [gene_symbols]

    records = list()
    for gene_symbol in gene_symbols:
        logger.debug("gene_symbol: {} type(gene_symbol): {}".format(gene_symbol, type(gene_symbol)))
        query = create_single_query_symbol(gene_symbol)
        results = send_api_request(query)
        logger.debug("gene_symbol: {} results: {}".format(gene_symbol, results))
        if 'data' in results and 'target' in results['data']:
            data = results['data']['target']
            if data != None:
                data['query'] = gene_symbol
                records.append(data)
            else:
                logger.error("Did not retrieve data from query for gene_symbol: {}".format(gene_symbol))
        else:
            logger.error("Did not retrieve data from query for gene_symbol: {}".format(gene_symbol))
    if len(records) != len(gene_symbols):
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
    # Pharos Target Info
    # POST https://pharos-api.ncats.io/graphql
    encoded_data = json.dumps(request)
    logger.info('PHAROS_API: request: {}'.format(request))

    try:
        response = requests.post(
            url="https://pharos-api.ncats.io/graphql",
            headers={
                "Content-Type": "application/json; charset=utf-8",
            },
            data=encoded_data
        )
        logger.info('PHAROS_API: Response HTTP Status Code: {status_code}'.format(
            status_code=response.status_code))

        # capture 400 error
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
        content = json.loads(response.content)
        return content

    except requests.exceptions.RequestException:
        logger.error('HTTP Request failed')
        return {}
