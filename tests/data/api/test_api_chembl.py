"""
test_api_chembl.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing chembl information
"""
from data.api import api_chembl
import json
import logging

# create logger
logger = logging.getLogger()
logger.level = logging.DEBUG

_test_target_id = 'CHEMBL1781870'
_test_target_id_list = ['CHEMBL1781870', 'CHEMBL2331064']
_test_molecule_id = 'CHEMBL3629617'
_test_molecule_id_list = ['CHEMBL3629617','CHEMBL1974235']
_test_assay_id = 'CHEMBL4274857'
_test_assay_id_list = ['CHEMBL4274857','CHEMBL1785363']
_test_assay_id_list_lg = ['CHEMBL904416', 'CHEMBL904417', 'CHEMBL939816', 'CHEMBL991464', 'CHEMBL991465',
                          'CHEMBL1052162', 'CHEMBL1052166', 'CHEMBL1226577', 'CHEMBL1226583', 'CHEMBL1226716',
                          'CHEMBL1249424', 'CHEMBL1249425', 'CHEMBL1249426', 'CHEMBL1249427', 'CHEMBL1249428' ]
_test_accession_id = 'Q13936'
_test_accession_id_list = ['Q13936', 'P07948']

# test the interface for simple calls
def test_api_interface():
    _test_name = 'test_api_interface'
    # check wrong type
    results_empty = api_chembl.send_api_request("nonsense", { })
    assert len(results_empty) == 0

    # test valid single target API call
    results_tgt = api_chembl.send_api_request("target/" + _test_target_id, { })
    logger.debug("{} {} Results: {} Type: {}".format(_test_name, "SINGLE TARGET", type(results_tgt), results_tgt))
    assert len(results_tgt) == 1
    assert ('cross_references' in results_tgt[0])
    assert results_tgt[0]['target_chembl_id'] == _test_target_id
    assert results_tgt[0]['organism'] == 'Homo sapiens'

    # test valid target list API call
    results_tgts = api_chembl.send_api_request("target", { 'target_chembl_id__in' : ','.join(_test_target_id_list)})
    logger.debug("{} {} Results: {} Type: {}".format(_test_name, "MULTIPLE TARGETS", type(results_tgt), results_tgt))
    assert len(results_tgts['targets']) == 2
    for tgt in results_tgts['targets']:
        assert (tgt['target_chembl_id'] in _test_target_id_list)
    #assert True == False

# test the interface for when need multiple calls due to api limit
def test_api_iteration():
    _test_name = 'test_api_interface'
    results = api_chembl.send_api_request("activity", { 'target_chembl_id__in' : ','.join(_test_target_id_list)})
    logger.debug("{} {} Results: {} Type: {}".format(_test_name, "MULTIPLE TARGETS", type(results), results))
    assert len(results['activities']) == 153

# test the interface for when need multiple calls due to api limit
def test_api_count():
    _test_name = 'test_api_count'
    results = api_chembl.send_api_request("activity", { 'target_chembl_id__in' : ','.join(_test_target_id_list)}, return_count=True)
    logger.debug("{} {} Results: {} Type: {}".format(_test_name, "MULTIPLE TARGETS -COUNTS", type(results), results))
    assert results == 153
