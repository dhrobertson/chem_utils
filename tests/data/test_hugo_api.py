"""
test_hugo_api.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing Pharos information
"""
from data import hugo_api
import json
import logging

# create logger
logger = logging.getLogger()
logger.level = logging.DEBUG

_test_target_id_results_ = {
    'hgnc_id' : 'HGNC:4486',
    'symbol' : 'GPR31',
    'status' : 'Approved'
}

_test_target_alias_results_ = {
    'hgnc_id' : 'HGNC:6079',
    'symbol' : 'INPP5D',
    'status' : 'Approved'
}

_test_target_id_ = "GPR31"
_test_target_alias_ = "SHIP1"
_test_target_ids_ = ["INPP5D", "SHIP1", "PAK1"]

def test_wrapper_functions():
    _test_name_ = 'test_wrapper_functions'

    # test single symbol
    results = hugo_api.get_target_info(_test_target_id_, "symbol")
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "SYMBOL", type(results), results))
    assert type(results) == type(list())
    assert len(results) == 1
    for key in _test_target_id_results_:
        assert key in results[0]
        assert results[0][key] == _test_target_id_results_[key]

    # test single alias - expect not found
    results = hugo_api.get_target_info(_test_target_id_, "alias_symbol")
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "SYMBOL", type(results), results))
    assert type(results) == type(list())
    assert len(results) == 1
    assert 'query' in results[0]
    assert results[0]['query'] == _test_target_id_

    # test single alias
    results = hugo_api.get_target_info(_test_target_alias_, "alias_symbol")
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "SYMBOL", type(results), results))
    assert type(results) == type(list())
    assert len(results) == 1
    for key in _test_target_alias_results_:
        assert key in results[0]
        assert results[0][key] == _test_target_alias_results_[key]

    # test single alias -- search both
    results = hugo_api.get_target_info(_test_target_alias_)
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "ALIAS", type(results), results))
    assert type(results) == type(list())
    assert len(results) == 1
    for key in _test_target_alias_results_:
        assert key in results[0]
        assert results[0][key] == _test_target_alias_results_[key]

# test the interface for simple calls
def test_api_interface():
    _test_name_ = 'test_api_interface'
    # check wrong type
    results_empty = hugo_api.send_api_request("nonsense")
    assert len(results_empty) == 0

    # test valid single target API call for symbol
    results_tgt = hugo_api.send_api_request("symbol/" + _test_target_id_)
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "SYMBOL", type(results_tgt), results_tgt))
    assert type(results_tgt) == type(list())
    assert len(results_tgt) == 1
    for key in _test_target_id_results_:
        assert key in results_tgt[0]
        assert results_tgt[0][key] == _test_target_id_results_[key]

    # test valid single target API call for symbol_alias
    results_tgt = hugo_api.send_api_request("alias_symbol/" + _test_target_alias_)
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "ALIAS", type(results_tgt), results_tgt))
    assert type(results_tgt) == type(list())
    assert len(results_tgt) == 1
    for key in _test_target_alias_results_:
        assert key in results_tgt[0]
        assert results_tgt[0][key] == _test_target_alias_results_[key]
