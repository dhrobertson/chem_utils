"""
test_chembl_api.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing chembl information
"""
from data.chembl_old import api
import sys
import logging

logger = logging.getLogger()
logger.level = logging.DEBUG

_test_target_id = 'CHEMBL2331064'
_test_target_id_list = ['CHEMBL1781870', 'CHEMBL2331064']

#TODO: Implement results cache

# test the interface for simple calls
def test_api_interface():
    _test_name = 'test_api_interface'
    # check wrong type
    results_empty = api.send_api_request("nonsense", { })
    assert len(results_empty) == 0

    # test valid single target API call
    results_tgt = api.send_api_request("target/" + _test_target_id, { })
    logger.debug("{} {} Results: {} Type: {}".format(_test_name, "SINGLE TARGET", type(results_tgt), results_tgt))
    assert len(results_tgt) == 1
    assert ('cross_references' in results_tgt[0])
    assert results_tgt[0]['target_chembl_id'] == _test_target_id
    assert results_tgt[0]['organism'] == 'Homo sapiens'

    # test valid target list API call
    results_tgts = api.send_api_request("target", { 'target_chembl_id__in' : ','.join(_test_target_id_list)})
    logger.debug("{} {} Results: {} Type: {}".format(_test_name, "MULTIPLE TARGETS", type(results_tgt), results_tgt))
    assert len(results_tgts['targets']) == 2
    for tgt in results_tgts['targets']:
        assert (tgt['target_chembl_id'] in _test_target_id_list)
    #assert True == False

# test the interface for when need multiple calls due to api limit
def test_api_iteration():
    _test_name = 'test_api_interface'
    results = api.send_api_request("activity", { 'target_chembl_id__in' : ','.join(_test_target_id_list)})
    logger.debug("{} {} Results: {} Type: {}".format(_test_name, "MULTIPLE TARGETS", type(results), results))
    assert len(results['activities']) == 153
