"""
test_core_utils.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions in the core_utils.py file
"""
from pathlib import Path
import csv
import json
from data import core_utils

# test file definitions
_null_file_ = 'my_file.txt'
_test_read_csv_file_ = str(Path(__file__).parent / 'test_data/table.csv')
_test_read_json_file_ = str(Path(__file__).parent / 'test_data/table.json')
_test_tmp_json_file_ = str(Path(__file__).parent / 'test_data/tmp/table.json')
_test_read_item_file_ = str(Path(__file__).parent / 'test_data/items.txt')
_test_content_file_ = str(Path(__file__).parent / 'test_data/tmp/content.txt')
_test_content_assay_results_json_file_ = str(Path(__file__).parent / 'test_data/assay_results.json')
_test_tmp_assay_results_file_ = str(Path(__file__).parent / 'test_data/tmp/assay_results.csv')

_test_content_ = "This is my test content"

#TODO: Add more of the full test functions
def test_file_exists():
    """ test the utility function that a file exists """
    assert core_utils.file_exists(_null_file_) == False
    assert core_utils.file_exists(_test_read_item_file_) == True

def test_read_and_write_content():
    """ test the functions to read and write content """
    if core_utils.file_exists(_test_content_file_):
        core_utils.remove_file(_test_content_file_)
    assert core_utils.file_exists(_test_content_file_) == False
    core_utils.save_text_to_file(_test_content_file_, _test_content_)
    assert core_utils.file_exists(_test_content_file_) == True
    content = core_utils.read_text_from_file(_test_content_file_)
    assert content == _test_content_
    core_utils.remove_file(_test_content_file_)
    assert core_utils.file_exists(_test_content_file_) == False
    if core_utils.file_exists(_test_content_file_):
        core_utils.remove_file(_test_content_file_)

def test_read_and_write_json():
    """ test the reading a json file """
    json_null = core_utils.read_json(_null_file_)
    assert len(json_null) == 0
    assert core_utils.file_exists(_test_read_json_file_) == True
    json_test = core_utils.read_json(_test_read_json_file_)
    assert len(json_test) == 3
    assert json_test[0]['id'] == 1
    assert json_test[0]['value'] == 20
    assert json_test[0]['type'] == 'integer'
    assert json_test[1]['id'] == 2
    assert json_test[1]['value'] == 'text, test'
    assert json_test[1]['type'] == 'string'
    assert json_test[2]['id'] == 3
    assert json_test[2]['value'] == 2.54
    assert json_test[2]['type'] == 'float'
    if core_utils.file_exists(_test_tmp_json_file_):
        core_utils.remove_file(_test_tmp_json_file_)
    assert core_utils.file_exists(_test_tmp_json_file_) == False
    core_utils.write_json(_test_tmp_json_file_, json_test)
    assert core_utils.file_exists(_test_tmp_json_file_) == True
    json_test_saved = core_utils.read_json(_test_tmp_json_file_)
    assert json_test_saved == json_test
    if core_utils.file_exists(_test_tmp_json_file_):
        core_utils.remove_file(_test_tmp_json_file_)

def test_read_csv_as_dict():
    """ test the reading a csv file as a dictionary """
    csv_null = core_utils.read_csv_as_dict(_null_file_)
    assert len(csv_null) == 0
    assert core_utils.file_exists(_test_read_csv_file_) == True
    csv_test = core_utils.read_csv_as_dict(_test_read_csv_file_)
    assert len(csv_test) == 3
    assert csv_test[0]['id'] == '1'
    assert csv_test[0]['value'] == '20'
    assert csv_test[0]['type'] == 'integer'
    assert csv_test[1]['id'] == '2'
    assert csv_test[1]['value'] == 'text, test'
    assert csv_test[1]['type'] == 'string'
    assert csv_test[2]['id'] == '3'
    assert csv_test[2]['value'] == '2.54'
    assert csv_test[2]['type'] == 'float'

def test_read_list_from_file():
    """ test the reading a list from a file """
    items_null = core_utils.read_list_from_file(_null_file_)
    assert len(items_null) == 0
    assert core_utils.file_exists(_test_read_item_file_) == True
    items_test = core_utils.read_list_from_file(_test_read_item_file_)
    assert len(items_test) == 6
    assert items_test[0] == 'item01'
    assert items_test[1] == 'item02'
    assert items_test[2] == 'item03'
    assert items_test[3] == 'item04'
    assert items_test[4] == 'item05'
    assert items_test[5] == 'item06'

def test_get_results_fields():
    """ test the only getting required fields from assay_results """
    assert core_utils.file_exists(_test_content_assay_results_json_file_) == True
    assay_results = core_utils.read_json(_test_content_assay_results_json_file_)
    results = core_utils.get_results_fields(assay_results[0])
    assert results['molecule_chembl_id'] == assay_results[0]['molecule_chembl_id']
    assert 'my_random_field' not in results

def test_write_assay_results():
    """ test the functions to write assay_results and the write_csv function """
    assert core_utils.file_exists(_test_content_assay_results_json_file_) == True
    _test_assay_results_ = core_utils.read_json(_test_content_assay_results_json_file_)
    if core_utils.file_exists(_test_tmp_assay_results_file_):
        core_utils.remove_file(_test_tmp_assay_results_file_)
    assert core_utils.file_exists(_test_tmp_assay_results_file_) == False

    core_utils.write_assay_results(_test_tmp_assay_results_file_, _test_assay_results_)
    assert core_utils.file_exists(_test_tmp_assay_results_file_) == True
    assay_results = core_utils.read_csv_as_dict(_test_tmp_assay_results_file_)
    assert len(assay_results) == 15
    assert assay_results[0]['molecule_chembl_id'] == _test_assay_results_[0]['molecule_chembl_id']
    assert assay_results[0]['canonical_smiles'] == 'CCCc1cc(=O)[nH]c(=S)[nH]1'
    if core_utils.file_exists(_test_tmp_assay_results_file_):
        core_utils.remove_file(_test_tmp_assay_results_file_)

def test_get_csvwriter():
    """ test getting csvwriter """
    csvwriter = core_utils.get_csvwriter("")
    assert csvwriter
