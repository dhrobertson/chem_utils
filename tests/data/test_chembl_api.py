"""
test_chembl_api.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing chembl information
"""
from data import chembl_api
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
_test_accession_id = 'Q13936'
_test_accession_id_list = ['Q13936', 'P07948']

def test_read_alias_map():
    """ test the reading of the mapping file to take gene name to target_chembl_id """
    alias_map = chembl_api.read_alias_map()
    keys = list(alias_map.keys())
    assert len(keys) == 20270

def test_get_target_details():
    """ test the functions to get the details for a target identified by target_chembl_id """
    # test for a single target
    target = chembl_api.get_target_details(_test_target_id)
    assert type(target) == type(list())
    assert target[0]['target_chembl_id'] == _test_target_id
    assert target[0]['target_type'] == 'SINGLE PROTEIN'
    assert target[0]['gene_symbol'] == 'INPP5D'

    # test for sending a list of targets
    targets = chembl_api.get_target_details(_test_target_id_list)
    assert type(targets) == type(list())
    assert len(targets) == 2
    for tgt in targets:
        assert (tgt['target_chembl_id'] in _test_target_id_list)

def test_find_target_by_accession():
    """ test the functions to get the details for a target identified by target_chembl_id """
    # test for a single accession number
    target = chembl_api.find_target_ids_by_accession(_test_accession_id)
    logger.debug('results: {}'.format(target))
    assert type(target) == type(list())
    assert len(target) == 5
    assert target[0]['target_chembl_id'] == 'CHEMBL1940'
    assert target[0]['target_type'] == 'SINGLE PROTEIN'

    # test for sending a list of accession numbers
    targets = chembl_api.find_target_ids_by_accession(_test_accession_id_list)
    assert type(targets) == type(list())
    assert len(targets) == 7

def test_get_molecule_details():
    """ test the functions to get the details for a molecular identified by molecule_chembl_id """
    # test for a single molecule
    molecule = chembl_api.get_molecule_details(_test_molecule_id)
    assert type(molecule) == type(list())
    assert len(molecule) == 1
    assert molecule[0]['molecule_chembl_id'] == _test_molecule_id
    assert molecule[0]['molecule_structures']['canonical_smiles'] == 'O=C(O)CC(O)(CC(=O)O)C(=O)O.OC(c1cc(-c2ccccc2)nc2c1ccc1ccccc12)C1CCCCN1'

    # test for list of molecules
    molecules = chembl_api.get_molecule_details(_test_molecule_id_list)
    assert type(molecules) == type(list())
    assert len(molecules) == 2
    for mol in molecules:
        assert (mol['molecule_chembl_id'] in _test_molecule_id_list)


def test_get_assay_details():
    """ test the functions to get the details for a assay identified by assay_chembl_id """
    # test for a single assay
    assay = chembl_api.get_assay_details(_test_assay_id)
    assert type(assay) == type(list())
    assert len(assay) == 1
    assert assay[0]['assay_chembl_id'] == _test_assay_id
    assert assay[0]['target_chembl_id'] == 'CHEMBL2331064'

    # test for a list of assays
    assay_list = chembl_api.get_assay_details(_test_assay_id_list)
    assert type(assay_list) == type(list())
    assert len(assay_list) == 2
    # TODO need to do this test so independent on order
    assert assay_list[1]['assay_chembl_id'] == _test_assay_id_list[0]
    assert assay_list[1]['target_chembl_id'] == 'CHEMBL2331064'
    assert assay_list[0]['assay_chembl_id'] == _test_assay_id_list[1]
    assert assay_list[0]['target_chembl_id'] == 'CHEMBL1781870'

def test_get_chembl_ids_for_targets():
    """ test getting the target_chembl_id from gene name """
    cid_1 = chembl_api.get_chembl_ids_for_targets('SHIP1')
    assert type(cid_1) == type(str())
    assert cid_1 == 'CHEMBL1781870'
    cid_list = chembl_api.get_chembl_ids_for_targets(['SHIP1', 'SHIP2', 'ABCA7', 'STAT3'])
    assert type(cid_list) == type(list())
    assert len(cid_list) == 4
    assert cid_list[1]['chembl_id'] == 'CHEMBL2331064'
    assert cid_list[2]['chembl_id'] == ''

def test_get_bioactivities_for_targets():
    """ test getting the bioactivities for target identified by target_chembl_id """
    records = chembl_api.get_bioactivities_for_targets(_test_target_id)
    assert len(records) == 16
    records_2 = chembl_api.get_bioactivities_for_targets(_test_target_id_list)
    assert len(records_2) >= 128

def test_get_bioactivities_for_molecules():
    """ test getting the bioactivities form molecule identified by molecule_chembl_id """
    m_records = chembl_api.get_bioactivities_for_molecules(_test_molecule_id)
    assert len(m_records) == 2
    m_records_2 = chembl_api.get_bioactivities_for_molecules(_test_molecule_id_list)
    assert len(m_records_2) >= 20

# test the interface for simple calls
def test_api_interface():
    _test_name = 'test_api_interface'
    # check wrong type
    results_empty = chembl_api.send_api_request("nonsense", { })
    assert len(results_empty) == 0

    # test valid single target API call
    results_tgt = chembl_api.send_api_request("target/" + _test_target_id, { })
    logger.debug("{} {} Results: {} Type: {}".format(_test_name, "SINGLE TARGET", type(results_tgt), results_tgt))
    assert len(results_tgt) == 1
    assert ('cross_references' in results_tgt[0])
    assert results_tgt[0]['target_chembl_id'] == _test_target_id
    assert results_tgt[0]['organism'] == 'Homo sapiens'

    # test valid target list API call
    results_tgts = chembl_api.send_api_request("target", { 'target_chembl_id__in' : ','.join(_test_target_id_list)})
    logger.debug("{} {} Results: {} Type: {}".format(_test_name, "MULTIPLE TARGETS", type(results_tgt), results_tgt))
    assert len(results_tgts['targets']) == 2
    for tgt in results_tgts['targets']:
        assert (tgt['target_chembl_id'] in _test_target_id_list)
    #assert True == False

# test the interface for when need multiple calls due to api limit
def test_api_iteration():
    _test_name = 'test_api_interface'
    results = chembl_api.send_api_request("activity", { 'target_chembl_id__in' : ','.join(_test_target_id_list)})
    logger.debug("{} {} Results: {} Type: {}".format(_test_name, "MULTIPLE TARGETS", type(results), results))
    assert len(results['activities']) == 153
