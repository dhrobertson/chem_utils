"""
HUGO_API.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with hugo data from API
"""
from pathlib import Path
import requests
import requests_cache
import json
from bs4 import BeautifulSoup
import logging

# setup requests_cache
requests_cache.install_cache()

# create logger
logger = logging.getLogger()

def get_target_info(genes, query_type=None):
    """ get the details for a target as symbols -- input either single item or a list """

    logger.debug("query_type: {} genes: {} type(genes): {}".format(query_type, genes, type(genes)))
    if type(genes) == type(str()):
        genes = [genes]

    try_alias_if_fail = False
    if query_type == None:
        try_alias_if_fail = True
        query_type = "symbol"
    logger.debug("query_type: {} genes: {} type(genes): {}".format(query_type, genes, type(genes)))

    records = list()
    for gene in genes:
        logger.debug("gene: {} type(gene) try_alias_if_fail {}: {}".format(gene, type(gene), try_alias_if_fail))
        results = send_api_request(query_type + '/' + gene)
        logger.debug("gene: {} results: {}".format(gene, results))
        if len(results) == 0 and try_alias_if_fail:
            logger.debug('Trying "alias_symbol" search for : {}'.format(gene))
            results = send_api_request("alias_symbol" + '/' + gene)
        if not results:
            logger.error('ERROR: Unable to locate information for specified gene: {}'.format(gene))
            results = [{ 'query' : gene}]
        for data in results:
            records.append(data)

    if len(records) != len(genes):
        logger.error("ERROR: Some gene_symbols not found?: returned {} targets for {} ids".format(
                len(records),len(genes)))

    return records


### key method to do the api request for Hugo
### Reference: https://www.genenames.org/help/rest/

def send_api_request(api_request):
    """
        function to send the api request to HUGO API
        sends the REST-based request and return json
        captures request error
        reference: https://www.genenames.org/help/rest/
    """
    _host_ = "rest.genenames.org"

    types = ['symbol', 'alias_symbol']
    query_type = api_request.split('/')[0]

    logger.info('HUGO_API: api_request: {} query_type: {}'.format(api_request, query_type))

    if query_type not in types:
        logger.error("api request of type \"{}\" not in allowed list: {}".format(query_type, types))
        return {}

    # need to check if request is from recursion
    url = "http://" + _host_ + "/fetch/" + api_request

    # send and parse the request
    try:
        response = requests.get(
            url=url,
            headers={
                "Accept": "application/json",
                "Content-Type": "multipart/form-data; charset=utf-8; boundary=__X_PAW_BOUNDARY__",
            },
        )
        logger.info('HUGO_API: Response HTTP Status Code: {status_code}'.format(
            status_code=response.status_code))
        logger.debug( 'Response HTTP Response Body: {content}'.format(
            content=response.text))
        # capture 400 errors
        if response.status_code == 404:
            logger.error("Request Error: Status 404 - page not found")
            return []

        if response.status_code == 400:
            logger.debug("400 Response: {}".format(response.text))
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
            logger.error("400 Request Error: {}".format(response.text))
            return {}

        # 200 correct response
        if response.status_code == 200:
            content = json.loads(response.text)
            content = content['response']
            logger.debug("query_type: {} content: {}".format(query_type, content))

            if 'numFound' not in content:
                logger.error("Appears no results found in response: {}".format(content))
                return {}

            if content['numFound'] == 1:
                return content['docs']

            logger.error("Do we need to handle multiple founds? Response: {}".format(content))
            return content['docs']

        logger.error("Unkown HTTP response code: {}".format(response.status_code))
        return {}

    except requests.exceptions.RequestException:
        logger.error('HTTP Request failed')
        return {}
