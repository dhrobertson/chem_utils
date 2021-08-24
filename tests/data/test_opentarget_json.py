"""
test_opentarget_parquet.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing OpenTarget information through parquet files
"""
from data import opentarget_json
import json
import logging

# create logger
logger = logging.getLogger()
logger.level = logging.DEBUG

_targets_dir_ = '/Users/dhrobertson/data/opentargets/json/targets'

_test_target_id_ = "ESR1"
_test_target_ids_ = ["INPP5D", "SHIP1", "PAK1"]

def test_core_functions():
    """ test sending the core api request to opentarget GraphQL """
    data_error = opentarget_json._read_json_files('bad_directory')
    assert len(data_error) == 0
    data = opentarget_json._read_json_files(_targets_dir_)
    assert len(data) >= 60608
    data_cache = opentarget_json._get_data_cache('bad_type')
    assert len(data_cache) == 0
    data_cache = opentarget_json._get_data_cache('targets')
    assert len(data_cache) >= 60608

def test_wrapper_functions():
    """ test the wrapper functions to return results in json for one or
        more targets """
    # check that works for single target as string
    results_single = opentarget_json.get_target_info(_test_target_id_)
    #logger.debug("Target: {} Results: {}".format(_test_target_id_, results_single))
    assert type(results_single) == type(list())
    assert len(results_single) == 1
    assert results_single[0]['id'] == 'ENSG00000091831'
    assert results_single[0]['approvedSymbol'] == _test_target_id_
    # do list
    results_list = opentarget_json.get_target_info(_test_target_ids_)
#    logger.debug("Target: {} Results: {}".format(_test_target_ids_, results))
    assert type(results_list) == type(list())
    assert len(results_list) == 2 # SHIP1 should fail -- matches INPP5D
