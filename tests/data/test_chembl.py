"""
test_chembl.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing chembl information
"""
from data import chembl
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

def test_read_alias_map():
    """ test the reading of the mapping file to take gene name to target_chembl_id """
    alias_map = chembl.read_alias_map()
    keys = list(alias_map.keys())
    assert len(keys) == 20270

def test_get_target_details():
    """ test the functions to get the details for a target identified by target_chembl_id """
    # test for a single target
    target = chembl.get_target_details(_test_target_id)
    assert type(target) == type(list())
    assert target[0]['target_chembl_id'] == _test_target_id
    assert target[0]['target_type'] == 'SINGLE PROTEIN'
    assert target[0]['gene_symbol'] == 'INPP5D'

    # test for sending a list of targets
    targets = chembl.get_target_details(_test_target_id_list)
    assert type(targets) == type(list())
    assert len(targets) == 2
    for tgt in targets:
        assert (tgt['target_chembl_id'] in _test_target_id_list)

def test_find_target_by_accession():
    """ test the functions to get the details for a target identified by target_chembl_id """
    # test for a single accession number
    target = chembl.find_target_ids_by_accession(_test_accession_id)
    logger.debug('results: {}'.format(target))
    assert type(target) == type(list())
    assert len(target) == 5
    assert target[0]['target_chembl_id'] == 'CHEMBL1940'
    assert target[0]['target_type'] == 'SINGLE PROTEIN'

    # test for sending a list of accession numbers
    targets = chembl.find_target_ids_by_accession(_test_accession_id_list)
    assert type(targets) == type(list())
    assert len(targets) == 7

def test_get_molecule_details():
    """ test the functions to get the details for a molecular identified by molecule_chembl_id """
    # test for a single molecule
    molecule = chembl.get_molecule_details(_test_molecule_id)
    assert type(molecule) == type(list())
    assert len(molecule) == 1
    assert molecule[0]['molecule_chembl_id'] == _test_molecule_id
    assert molecule[0]['molecule_structures']['canonical_smiles'] == 'O=C(O)CC(O)(CC(=O)O)C(=O)O.OC(c1cc(-c2ccccc2)nc2c1ccc1ccccc12)C1CCCCN1'

    # test for list of molecules
    molecules = chembl.get_molecule_details(_test_molecule_id_list)
    assert type(molecules) == type(list())
    assert len(molecules) == 2
    for mol in molecules:
        assert (mol['molecule_chembl_id'] in _test_molecule_id_list)

def test_get_assay_details():
    """ test the functions to get the details for a assay identified by assay_chembl_id """
    logger.level = logging.DEBUG
    # test for a single assay
    assay = chembl.get_assay_details(_test_assay_id)
    assert type(assay) != type(list())
    assert assay['assay_chembl_id'] == _test_assay_id
    assert assay['target_chembl_id'] == 'CHEMBL2331064'
    assert assay['results_count'] == 51

    # test for a list of assays
    assert type(_test_assay_id_list) == type(list())
    assert len(_test_assay_id_list) == 2
    assay_list = chembl.get_assay_details(_test_assay_id_list)
    assert type(_test_assay_id_list) == type(list())
    assert len(_test_assay_id_list) == 2
    assert type(assay_list) == type(list())
    assert len(assay_list) == 2

    # TODO need to do this test so independent on order
    print(assay_list)
    assert assay_list[1]['assay_chembl_id'] == _test_assay_id_list[0]
    assert assay_list[1]['target_chembl_id'] == 'CHEMBL2331064'
    assert assay_list[0]['assay_chembl_id'] == _test_assay_id_list[1]
    assert assay_list[0]['target_chembl_id'] == 'CHEMBL1781870'

    # try a large list
    # test for a list of assays
    assay_list_lg = chembl.get_assay_details(_test_assay_id_list_lg, chunk_size=5)
    assert type(assay_list_lg) == type(list())
    assert len(assay_list_lg) == 15

def test_get_chembl_ids_for_targets():
    """ test getting the target_chembl_id from gene name """
    cid_1 = chembl.get_chembl_ids_for_targets('SHIP1')
    assert type(cid_1) == type(str())
    assert cid_1 == 'CHEMBL1781870'
    cid_list = chembl.get_chembl_ids_for_targets(['SHIP1', 'SHIP2', 'ABCA7', 'STAT3'])
    assert type(cid_list) == type(list())
    assert len(cid_list) == 4
    assert cid_list[1]['chembl_id'] == 'CHEMBL2331064'
    assert cid_list[2]['chembl_id'] == ''

def test_get_assays_for_targets():
    """ test getting the bioactivities for target identified by target_chembl_id """
    records = chembl.get_assays_for_targets(_test_target_id)
    assert type(records) == type(list())
    assert len(records) == 3
    records_2 = chembl.get_assays_for_targets(_test_target_id_list)
    assert type(records_2) == type(dict())
    assert len(records_2[_test_target_id_list[0]]) == 3
    assert len(records_2[_test_target_id_list[1]]) == 18

def test_get_bioactivities_for_targets():
    """ test getting the bioactivities for target identified by target_chembl_id """
    records = chembl.get_bioactivities_for_targets(_test_target_id)
    assert len(records) == 16
    records_2 = chembl.get_bioactivities_for_targets(_test_target_id_list)
    assert len(records_2) >= 128

def test_get_bioactivities_for_molecules():
    """ test getting the bioactivities form molecule identified by molecule_chembl_id """
    m_records = chembl.get_bioactivities_for_molecules(_test_molecule_id)
    assert len(m_records) == 2
    m_records_2 = chembl.get_bioactivities_for_molecules(_test_molecule_id_list)
    assert len(m_records_2) >= 20
