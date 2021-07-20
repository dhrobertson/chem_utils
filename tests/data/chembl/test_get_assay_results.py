"""
test_get_assay_results.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing chembl information

Reminder set path export PYTHONPATH=~/projects/local/chembl_webresource_client:~/projects/chembl_utils
"""
from pathlib import Path
import pprint
import logging
import json
from utils import core_utils
from data.chembl_old import get_assay_results

# test file definitions
_test_content_assay_results_json_file_ = str(Path(__file__).parent / 'data/assay_results.json')
assert core_utils.file_exists(_test_content_assay_results_json_file_) == True
_test_assay_results_ = core_utils.read_json(_test_content_assay_results_json_file_)

def test_by_targets():
    """ test getting the active molecules from target_chembl_id """
    # 2021-02-03
    # CHEMBL1839 (TPO) - 13 bioactivities
    # CHEMBL5692 (ampC) - 4 bioactivities
    logger = logging.getLogger()
    org_level = logger.getEffectiveLevel()

    logger.level = logging.ERROR
    results_1 = get_assay_results.by_targets(['CHEMBL1839'])
    assert len(results_1) >= 13

    results_2 = get_assay_results.by_targets(['CHEMBL1839','CHEMBL5692'])
    assert len(results_2) >= 17
    logger.level = org_level

def test_by_molecules():
    """ test getting the active compounds from molecule_chembl_id """
    # 2021-02-03
    # CHEMBL207881 (GW-9508) - 90 bioactivities
    # CHEMBL203297 () - 2 bioactivities
    logger = logging.getLogger()
    org_level = logger.getEffectiveLevel()

    logger.level = logging.ERROR
    results_1 = get_assay_results.by_molecules(['CHEMBL207881'])
    assert len(results_1) >= 94

    results_2 = get_assay_results.by_molecules(['CHEMBL207881','CHEMBL203297'])
    assert len(results_2) >= 92
    logger.level = org_level

def test_get_targets():
    """ test getting a list of the targets from the bioactivities """
    targets = get_assay_results.get_targets(_test_assay_results_)
    assert type(targets) == type(list())
    assert len(targets) == 1
    assert targets[0]['id'] == 'CHEMBL1839'
    assert targets[0]['gene_symbol'] == 'TPO'
    assert targets[0]['number_assays'] == 4
    assert targets[0]['number_molecules'] == 13
    assert targets[0]['number_results'] == 15

def test_get_molecules():
    """ test getting a list of the molecules from the bioactivities """
    molecules = get_assay_results.get_molecules(_test_assay_results_)
    assert type(molecules) == type(list())
    assert len(molecules) == 13
    assert molecules[0]['id'] == 'CHEMBL1518'
    assert molecules[0]['number_targets'] == 1
    assert molecules[0]['number_assays'] == 2
    assert molecules[0]['number_results'] == 2

def test_get_assays():
    """ test getting a list of the assays from the bioactivities """
    assays = get_assay_results.get_assays(_test_assay_results_)
    assert type(assays) == type(list())
    assert len(assays) == 4
    assert assays[0]['id'] == 'CHEMBL3636368'
    assert assays[0]['target'] == 'CHEMBL1839'
    assert assays[0]['number_molecules'] == 1
    assert assays[0]['number_results'] == 1

def test_get_results():
    """ test getting a list of the results from the bioactivities """
    results = get_assay_results.get_results(_test_assay_results_)
    assert type(results) == type(list())
    assert len(results) == 15
    assert results[0]['id'] == 15771999

def test_convert_to_standard():
    """ test converting bioactivities to a standard researchDataSet """
    std_assay = get_assay_results.convert_to_standard(_test_assay_results_)
    assert type(std_assay) == type(dict())
    assert 'targets' in std_assay
    assert len(std_assay['targets']) == 1
    assert std_assay['targets'][0]['id'] == 'CHEMBL1839'
    assert std_assay['targets'][0]['number_assays'] == 4
    assert std_assay['targets'][0]['number_molecules'] == 13
    assert std_assay['targets'][0]['number_results'] == 15
    assert 'molecules' in std_assay
    assert len(std_assay['molecules']) == 13
    assert std_assay['molecules'][0]['id'] == 'CHEMBL1518'
    assert std_assay['molecules'][0]['number_targets'] == 1
    assert std_assay['molecules'][0]['number_assays'] == 2
    assert std_assay['molecules'][0]['number_results'] == 2
    assert 'assays' in std_assay
    assert len(std_assay['assays']) == 4
    assert std_assay['assays'][0]['id'] == 'CHEMBL3636368'
    assert std_assay['assays'][0]['target'] == 'CHEMBL1839'
    assert std_assay['assays'][0]['number_molecules'] == 1
    assert std_assay['assays'][0]['number_results'] == 1
    assert 'results' in std_assay
    assert len(std_assay['results']) == 15
    assert std_assay['results'][0]['id'] == 15771999
    #core_utils.write_json('dataset.json', std_assay)
