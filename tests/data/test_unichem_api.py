"""
test_unichem_api.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing Unichem information
"""
from data import unichem_api
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

_test_src_id_ = "1"
_test_src_id_list_ = ["1", "2"]
_test_chembl_id_ = "CHEMBL12"
_test_chembl_id_list_ = ["CHEMBL12", "CHEMBL22"]


def test_util_functions():
    _test_name_ = 'test_util_functions'

    # test single value parameter
    (is_list, ids_list) = unichem_api.as_list(_test_src_id_)
    assert is_list == False
    assert len(ids_list) == 1

    # test list parameters
    (is_list, ids_list) = unichem_api.as_list(_test_src_id_list_)
    assert is_list == True
    assert len(ids_list) == 2
    assert ids_list[0] == _test_src_id_list_[0]
    assert ids_list[1] == _test_src_id_list_[1]

def test_wrapper_functions():
    _test_name_ = 'test_wrapper_functions'

    # test all srcs
    all_src_ids = unichem_api.get_all_src_ids()
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "get_all_src_ids", type(all_src_ids), all_src_ids))
    assert type(all_src_ids) == type(list())
    assert len(all_src_ids) == 48

    # test single src
    src_details = unichem_api.get_src_details(_test_src_id_)
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "get_src_details", type(src_details), src_details))
    assert type(src_details) == type(dict())
    assert src_details['name'] == 'chembl'

    # test single from one source
    src_details_list = unichem_api.get_src_details(_test_src_id_list_)
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "get_id_mappings", type(src_details_list), src_details_list))
    assert type(src_details_list) == type(list())
    assert src_details_list[0]['name'] == 'chembl'
    assert src_details_list[1]['name'] == 'drugbank'

    # test single from one source
    cpd_mappings = unichem_api.get_id_mappings(_test_chembl_id_, '1')
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "get_id_mappings", type(cpd_mappings), cpd_mappings))
    assert type(cpd_mappings) == type(list())
    assert cpd_mappings[0]['src_id'] == '1'
    assert cpd_mappings[0]['src_compound_id'] == _test_chembl_id_
    assert cpd_mappings[1]['src_id'] == '2'
    assert cpd_mappings[1]['src_compound_id'] == 'DB00829'

    # test single from one source to another
    cpd_mappings_2 = unichem_api.get_id_mappings(_test_chembl_id_, '1', '2')
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "get_id_mappings", type(cpd_mappings_2), cpd_mappings_2))
    assert type(cpd_mappings_2) == type(list())
    assert len(cpd_mappings_2) == 2
    assert cpd_mappings_2[0]['src_compound_id'] == 'DB00829'
    assert cpd_mappings_2[1]['src_compound_id'] == 'DB07699'

    # test list from one source
    cpd_mappings_list = unichem_api.get_id_mappings(_test_chembl_id_list_, '1')
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "get_id_mappings", type(cpd_mappings_list), cpd_mappings_list))
    assert type(cpd_mappings_list) == type(dict())
    assert cpd_mappings_list[_test_chembl_id_list_[0]][0]['src_id'] == '1'
    assert cpd_mappings_list[_test_chembl_id_list_[0]][0]['src_compound_id'] == _test_chembl_id_
    assert cpd_mappings_list[_test_chembl_id_list_[1]][1]['src_id'] == '2'
    assert cpd_mappings_list[_test_chembl_id_list_[1]][1]['src_compound_id'] == 'DB00440'

    # test list from one source to another
    cpd_mappings_list_2 = unichem_api.get_id_mappings(_test_chembl_id_list_, '1', '2')
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "get_id_mappings", type(cpd_mappings_list_2), cpd_mappings_list_2))
    assert type(cpd_mappings_list) == type(dict())
    assert cpd_mappings_list_2[_test_chembl_id_list_[0]][0]['src_compound_id'] == 'DB00829'
    assert cpd_mappings_list_2[_test_chembl_id_list_[1]][0]['src_compound_id'] == 'DB00440'

    #assert True == False

# test the interface for simple calls
def test_api_interface():
    _test_name_ = 'test_api_interface'
    # check wrong type
    results_empty = unichem_api.send_api_request("nonsense")
    assert len(results_empty) == 0

    # test valid single target API call for symbol
    results_all_ids = unichem_api.send_api_request("src_ids/")
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "src_ids", type(results_all_ids), results_all_ids))
    assert type(results_all_ids) == type(list())
    assert len(results_all_ids) == 48

    # test valid single target API call for symbol
    src_1_details = unichem_api.send_api_request("sources/1")
    logger.debug("{} {} Results: {} Type: {}".format(_test_name_, "sources/1", type(src_1_details), src_1_details))
    assert type(src_1_details) == type(list())
    assert len(src_1_details) == 1
    assert src_1_details[0]['name'] == 'chembl'
