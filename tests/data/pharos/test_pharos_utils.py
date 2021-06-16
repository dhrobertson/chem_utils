"""
test_utils.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing Pharos information
"""
from data.pharos import utils

test_query_str = "{\"query\": \"query targetDetails{target(q:{sym:\\\"GPR31\\\"}){name tdl fam sym description novelty}}\",\"variables\": {}}"
test_query = {
        "query": "query targetDetails{target(q:{sym:\"GPR31\"}){name tdl fam sym description novelty}}",
        "variables": {}
    }
test_query_error = {
        "query": "query targetDetails{target(q:{sym:\\\"GPR31\\\"}){name tdl fam sym description novelty}}",
        "variables": {}
    }
test_query_results = {
    "data": {
        "target": {
            "name" : "12-(S)-hydroxy-5,8,10,14-eicosatetraenoic acid receptor",
            "tdl" : "Tbio",
            "fam" : "GPCR",
            "sym" : "GPR31" ,
            "description" : None,
            "novelty" : 0.09160693
        }
    }
}

def test_core_functions():
    """ test sending the core api request to Pharos GraphQL """
    results = utils.send_api_request(test_query)
    assert results == test_query_results
    results_error = utils.send_api_request(test_query_error)
    assert results_error == {}
    results_fail = utils.send_api_request("ill_formed_query")
    assert results_fail == {}

    #assert True == False
