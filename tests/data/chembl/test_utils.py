"""
test_utils.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing chembl information
"""
from data.chembl import utils

def test_read_alias_map():
    """ test the reading of the mapping file to take gene name to target_chembl_id """
    alias_map = utils.read_alias_map()
    keys = list(alias_map.keys())
    assert len(keys) == 20270

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
    records = utils.get_bioactivities_for_targets(['CHEMBL1781870'])
    assert len(records) == 16
    #print(pprint.pformat(records))
    records_2 = utils.get_bioactivities_for_targets(['CHEMBL1781870','CHEMBL2331064'])
    assert len(records_2) >= 128
    #print(pprint.pformat(records_2))
    m_records = utils.get_bioactivities_for_molecules(['CHEMBL3629617'])
    assert len(m_records) == 2
    #print(pprint.pformat(m_records))
    m_records_2 = utils.get_bioactivities_for_molecules(['CHEMBL1974235'])
    assert len(m_records_2) == 62
    #print(pprint.pformat(m_records_2))

def test_get_target_details():
    """ test the functions to get the details for a target identified by target_chembl_id """
    # test for a single target
    target = utils.get_target_details('CHEMBL1781870')
    assert type(target) == type(dict())
    assert target['target_chembl_id'] == 'CHEMBL1781870'
    assert target['target_type'] == 'SINGLE PROTEIN'
    assert target['gene_symbol'] == 'INPP5D'
    #print(pprint.pformat(target))

    # test for sending a list of targets
    targets = utils.get_target_details(['CHEMBL1781870','CHEMBL2331064'])
    assert type(targets) == type(list())
    assert len(targets) == 2
    assert targets[0]['target_chembl_id'] == 'CHEMBL1781870'
    assert targets[1]['target_chembl_id'] == 'CHEMBL2331064'
    #print(pprint.pformat(targets))

def test_get_molecule_details():
    """ test the functions to get the details for a molecular identified by molecule_chembl_id """
    molecule = utils.get_molecule_details('CHEMBL3629617')
    assert type(molecule) == type(dict())
    assert molecule['molecule_chembl_id'] == 'CHEMBL3629617'
    assert molecule['molecule_structures']['canonical_smiles'] == 'O=C(O)CC(O)(CC(=O)O)C(=O)O.OC(c1cc(-c2ccccc2)nc2c1ccc1ccccc12)C1CCCCN1'
    #print(pprint.pformat(molecule))

    molecules = utils.get_molecule_details(['CHEMBL3629617','CHEMBL1974235'])
    assert type(molecules) == type(list())
    assert len(molecules) == 2
    assert molecules[0]['molecule_chembl_id'] == 'CHEMBL3629617'
    assert molecules[1]['molecule_chembl_id'] == 'CHEMBL1974235'
    #print(pprint.pformat(molecules))

def test_get_assay_details():
    assay = utils.get_assay_details('CHEMBL4274857')
    assert type(assay) == type(dict())
    assert assay['assay_chembl_id'] == 'CHEMBL4274857'
    assert assay['target_chembl_id'] == 'CHEMBL2331064'
    #print(pprint.pformat(assay))

    assay_list = utils.get_assay_details(['CHEMBL4274857','CHEMBL4274857'])
    assert type(assay_list) == type(list())
    assert len(assay_list) == 2
    assert assay_list[0]['assay_chembl_id'] == 'CHEMBL4274857'
    assert assay_list[0]['target_chembl_id'] == 'CHEMBL2331064'
    assert assay_list[1]['assay_chembl_id'] == 'CHEMBL4274857'
    assert assay_list[1]['target_chembl_id'] == 'CHEMBL2331064'
    #print(pprint.pformat(assay_list))
