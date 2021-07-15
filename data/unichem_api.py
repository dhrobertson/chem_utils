"""
unichem_api.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with unichem data from API
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

# helper function for handling inputs
def as_list(ids):
    ids_list = list()
    is_list = True
    if type(ids) == type(list()):
        for id in ids:
            ids_list.append(id)
    else:
        is_list = False
        ids_list.append(ids)

    return (is_list, ids_list)

def get_all_src_ids():
    """ get all the available src_ids from unichem """
    ids_list = list()

    logger.debug("query_type: {}".format("src_ids"))
    api_results = send_api_request('src_ids/')
    for result in api_results:
        if 'src_id' in result:
            ids_list.append(result['src_id'])
        else:
            logger.warning("can't locate src_id key in result: {}".format(result))

    return ids_list

def get_src_details(ids):
    """ get details on the src_ids from unichem - single or list"""

    (is_list, ids_list) = as_list(ids)

    results = list()
    logger.debug("query_type: {} # ids: {} is_list: {} ids_list: {}".format("sources", len(ids_list), is_list, ids_list))
    for id in ids_list:
        api_results = send_api_request('sources/' + id)
        for result in api_results:
            results.append(result)

    if not is_list:
        return results[0]
    return results

def get_id_mappings(ids, type=1, type_to=None):
    """ get id_mapping from unichem - single or list"""

    (is_list, ids_list) = as_list(ids)
    query = "src_compound_id_all"

    results = dict()
    logger.debug("query_type: {} # ids: {} is_list: {} ids_list: {}".format(query, len(ids_list), is_list, ids_list))
    for id in ids_list:
        api_query = '/'.join([query, id, type])
        if type_to:
            api_query = '/'.join([api_query, type_to])
        results[id] = send_api_request(api_query)

    if len(ids_list) > 1:
        return results

    return results[ids_list[0]]

### key method to do the api request for Hugo
### Reference: https://www.ebi.ac.uk/unichem/rest
def send_api_request(api_request):
    """
        function to send the api request to unichem API
        sends the REST-based request and return json
        captures request error
        reference: https://www.ebi.ac.uk/unichem/info/webservices#GetSrcCpdIdsAll
    """
    _host_ = "www.ebi.ac.uk/unichem/rest/"

    valid_query_types = ['src_ids', 'sources', "src_compound_id_all"]
    query_type = api_request.split('/')[0]

    logger.info('UNICHEM_API: api_request: {} query_type: {}'.format(api_request, query_type))

    if query_type not in valid_query_types:
        logger.error("api request of type \"{}\" not in allowed list: {}".format(query_type, valid_query_types))
        return {}

    # need to check if request is from recursion
    url = "http://" + _host_ + api_request

    # send and parse the request
    try:
        response = requests.get(
            url=url,
            headers={
                "Accept": "application/json",
                "Content-Type": "multipart/form-data; charset=utf-8; boundary=__X_PAW_BOUNDARY__",
            },
        )
        logger.info('UNICHEM_API: Response HTTP Status Code: {status_code}'.format(
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
            #content = content['response']
            logger.debug("query_type: {} content: {}".format(query_type, content))
            return content

        logger.error("Unkown HTTP response code: {}".format(response.status_code))
        return {}

    except requests.exceptions.RequestException:
        logger.error('HTTP Request failed')
        return {}
