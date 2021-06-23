"""
api.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific function to access the chembl api
"""
import requests
import requests_cache
import json
import logging

# setup requests_cache
requests_cache.install_cache()

# create logger and set logger level
logger = logging.getLogger()
#logger.level = logging.INFO

### key method to do the api request for chembl
### Reference: https://www.ebi.ac.uk/chembl/api/data/docs

def send_api_request(api_request, params, level=0, _results_list=None):
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

    logger.debug('api_request: {} parameters: {}'.format(api_request, params))

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
    logger.info('Sending API request to {} of type: {} level: {}'.format(url, u_type, level, len(_results_list)))

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
        logger.info('Response HTTP Status Code: {status_code}'.format(
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

            for item in content[types_map[u_type]]:
                _results_list.append(item)
            logger.info("Results thus far: level: {} count: {} of {}".format(level, len(_results_list), total_count))

            # if page_meta['next'] != null -- continue recursion
            if content['page_meta']['next']:
                total = content['page_meta']['total_count']
                level = level + 1
                logger.debug("Sending new request: {}".format(content['page_meta']['next']))
                return send_api_request(content['page_meta']['next'], {}, level=level, _results_list=_results_list)

            # end of recursion return compiled results
            return { 'page_meta':{}, types_map[u_type]: _results_list }

        logger.error("Unkown HTTP response code: {}".format(response.status_code))
        return {}

    except requests.exceptions.RequestException:
        logger.error('HTTP Request failed')
        return {}
