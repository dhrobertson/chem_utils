"""
test_pharos_pharos_api.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing Pharos information
"""
from data import pharos_api
import json
import logging

# create logger
logger = logging.getLogger()
logger.level = logging.DEBUG

test_query = {
        "query": "query targetDetails{target(q:{sym:\"GPR31\"}){name tdl fam sym description novelty diseaseCounts {name value} ligandCounts {name value}}}",
        "variables": {}
    }
test_query_error = {
        "query": "query targetDetails{target(q:{sym:\\\"GPR31\\\"}){name tdl fam sym description novelty}}",
        "variables": {}
    }
test_query_results = {
    "data": {
        "target": {
            "name": "12-(S)-hydroxy-5,8,10,14-eicosatetraenoic acid receptor",
            "tdl": "Tbio",
            "fam": "GPCR",
            "sym": "GPR31",
            "description": None,
            "novelty": 0.09160693,
            "diseaseCounts": [
                {
                    "name": "diabetes mellitus",
                    "value": 3667
                },
                {
                    "name": "rheumatoid arthritis",
                    "value": 1510
                }
            ],
                "ligandCounts": [
                    {
                        "name": "ligand",
                        "value": 0
                    },
                    {
                        "name": "drug",
                        "value": 0
                    }
                ]
            }
        }
    }

_test_target_id_ = "GPR31"
_test_target_ids_ = ["INPP5D", "SHIP1", "PAK1"]
DEBUG = False

def test_core_functions():
    """ test sending the core api request to Pharos GraphQL """
    results = pharos_api.send_api_request(test_query)
    #print(results)
    assert results == test_query_results
    results_error = pharos_api.send_api_request(test_query_error)
    assert results_error == {}
    results_fail = pharos_api.send_api_request("ill_formed_query")
    assert results_fail == {}

def test_query_functions():
    """ test functions to create the GraphQL querie strings """
    # query to look for symbol
    query = pharos_api.create_single_query_symbol(_test_target_id_)
    _new_query_ = json.dumps(query, sort_keys=True)
    _compare_query_ = json.dumps(test_query, sort_keys=True)
    assert _new_query_ == _compare_query_
    logger.debug("Target: {} query: {}".format(_test_target_id_, query))

def test_wrapper_functions():
    """ test the wrapper functions to return results in json for one or
        more targets """
    # check that works for single target as string
    results = pharos_api.get_target_info(_test_target_id_)
    if DEBUG:
        logger.debug("Target: {} Results: {}".format(_test_target_id_, results))
        logger.debug("Target: {} Test: {}".format(_test_target_id_, test_query_results))
    assert type(results) == type(list())
    assert len(results) == 1
    for key in test_query_results['data']['target']:
        assert key in results[0]
        assert results[0][key] == test_query_results['data']['target'][key]

    # check that works for the list
    results = pharos_api.get_target_info(_test_target_ids_)
    logger.debug("Target: {} Results: {}".format(_test_target_ids_, results))
    assert type(results) == type(list())
    assert len(results) == 2 # SHIP1 should fail
