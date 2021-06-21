"""
test_utils.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing chembl information
"""
from data.chembl import utils

_test_target_id = 'CHEMBL1781870'
_test_target_id_list = ['CHEMBL1781870', 'CHEMBL2331064']
_test_molecule_id = 'CHEMBL3629617'
_test_molecule_id_list = ['CHEMBL3629617','CHEMBL1974235']
_test_assay_id = 'CHEMBL4274857'
_test_assay_id_list = ['CHEMBL4274857','CHEMBL1785363']

def test_read_alias_map():
    """ test the reading of the mapping file to take gene name to target_chembl_id """
    alias_map = utils.read_alias_map()
    keys = list(alias_map.keys())
    assert len(keys) == 20270

def test_api_interface():
    results_empty = utils.send_api_request("nonsense/", { })
    assert results_empty == {}
    results = utils.send_api_request("target/" + _test_target_id, { })
    assert ('cross_references' in results)
    assert results['target_chembl_id'] == _test_target_id
    assert results['organism'] == 'Homo sapiens'

def test_get_target_details():
    """ test the functions to get the details for a target identified by target_chembl_id """
    # test for a single target
    target = utils.get_target_details(_test_target_id)
    assert type(target) == type(dict())
    assert target['target_chembl_id'] == _test_target_id
    assert target['target_type'] == 'SINGLE PROTEIN'
    assert target['gene_symbol'] == 'INPP5D'

    # test for sending a list of targets
    targets = utils.get_target_details(_test_target_id_list)
    assert type(targets) == type(list())
    assert len(targets) == 2
    assert targets[0]['target_chembl_id'] == _test_target_id_list[0]
    assert targets[1]['target_chembl_id'] == _test_target_id_list[1]

def test_get_molecule_details():
    """ test the functions to get the details for a molecular identified by molecule_chembl_id """
    # test for a single molecule
    molecule = utils.get_molecule_details(_test_molecule_id)
    assert type(molecule) == type(dict())
    assert molecule['molecule_chembl_id'] == _test_molecule_id
    assert molecule['molecule_structures']['canonical_smiles'] == 'O=C(O)CC(O)(CC(=O)O)C(=O)O.OC(c1cc(-c2ccccc2)nc2c1ccc1ccccc12)C1CCCCN1'

    # test for list of molecules
    molecules = utils.get_molecule_details(_test_molecule_id_list)
    assert type(molecules) == type(list())
    assert len(molecules) == 2
    assert molecules[1]['molecule_chembl_id'] == _test_molecule_id_list[0]
    assert molecules[0]['molecule_chembl_id'] == _test_molecule_id_list[1]

def test_get_assay_details():
    """ test the functions to get the details for a assay identified by assay_chembl_id """
    # test for a single assay
    assay = utils.get_assay_details(_test_assay_id)
    assert type(assay) == type(dict())
    assert assay['assay_chembl_id'] == _test_assay_id
    assert assay['target_chembl_id'] == 'CHEMBL2331064'

    # test for a list of assays
    assay_list = utils.get_assay_details(_test_assay_id_list)
    assert type(assay_list) == type(list())
    assert len(assay_list) == 2
    # TODO need to do this test so independent on order
    assert assay_list[1]['assay_chembl_id'] == _test_assay_id_list[0]
    assert assay_list[1]['target_chembl_id'] == 'CHEMBL2331064'
    assert assay_list[0]['assay_chembl_id'] == _test_assay_id_list[1]
    assert assay_list[0]['target_chembl_id'] == 'CHEMBL1781870'

def test_get_chembl_ids_for_targets():
    """ test getting the target_chembl_id from gene name """
    cid_1 = utils.get_chembl_ids_for_targets('SHIP1')
    assert type(cid_1) == type(str())
    assert cid_1 == 'CHEMBL1781870'
    cid_list = utils.get_chembl_ids_for_targets(['SHIP1', 'SHIP2', 'ABCA7', 'STAT3'])
    assert type(cid_list) == type(list())
    assert len(cid_list) == 4
    assert cid_list[1]['chembl_id'] == 'CHEMBL2331064'
    assert cid_list[2]['chembl_id'] == ''

def test_get_bioactivities_for_targets():
    """ test getting the bioactivities for target identified by target_chembl_id """
    records = utils.get_bioactivities_for_targets(_test_target_id)
    assert len(records) == 16
    records_2 = utils.get_bioactivities_for_targets(_test_target_id_list)
    assert len(records_2) >= 128

def test_get_bioactivities_for_molecules():
    """ test getting the bioactivities form molecule identified by molecule_chembl_id """
    m_records = utils.get_bioactivities_for_molecules(_test_molecule_id)
    assert len(m_records) == 2
    m_records_2 = utils.get_bioactivities_for_molecules(_test_molecule_id_list)
    assert len(m_records_2) >= 20
