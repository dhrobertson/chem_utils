"""
utils.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with chembl data
"""
from pathlib import Path
import requests
import json
from bs4 import BeautifulSoup
import os
import sys
import logging
import csv
import pprint

#
# define helper files
#
#_g2c_map_file_ = Path(__file__).parent / 'data/gene2chembl_id.csv'

# create logger
# TODO -- fix later -- check if logger exists -- and than attach to it
logger = logging.getLogger()
logger.level = logging.INFO
stream_handler = logging.StreamHandler(sys.stderr)
logger.addHandler(stream_handler)

# Install the Python Requests library:
# `pip install requests`

import requests

def send_api_request(request):
    """
        function to send the api request to Pharos API
        sends the query request and return json
        captures request error
    """
    # Pharos Target Info
    # POST https://pharos-api.ncats.io/graphql
    encoded_data = json.dumps(request)

    try:
        response = requests.post(
            url="https://pharos-api.ncats.io/graphql",
            headers={
                "Content-Type": "application/json; charset=utf-8",
            },
            data=encoded_data
        )
        logger.info('Response HTTP Status Code: {status_code}'.format(
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
