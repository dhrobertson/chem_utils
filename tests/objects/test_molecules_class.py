"""
test_compare_molecules.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the Molecule list class
"""
import pprint
import logging
from pathlib import Path
from objects import Molecules

_null_file_ = 'data/missing.smi' # should not exist
_test_smi_file_ = str(Path(__file__).parent / 'data/test_read.smi')
_test_smi_2_file_ = str(Path(__file__).parent / 'data/test_read_2.smi')
_test_smi_tmp_file_ = str(Path(__file__).parent / 'data/tmp/test_read.smi')

_test_smiles_1_ = 'C1=CC=CN=C1'
_test_smiles_1_tup = (_test_smiles_1_, 'mol_1')
_test_smiles_2_ = 'c1cccnc1'
_test_smiles_2_tup = (_test_smiles_2_, 'mol_2')
_test_smiles_error_ = 'n1cccc1'
_test_smiles_error_tup = (_test_smiles_error_, 'mol_error')
_test_smiles_3_ = 'n1cnccc1C'
_test_smiles_3_tup = (_test_smiles_3_, 'mol_3')

logger = logging.getLogger()
logger.level = logging.ERROR


#TODO: Add more of the full test functions

def test_new_class():
    """ test the creation and add functions """
    molecules = Molecules.Molecules()
    assert type(molecules) == type(Molecules.Molecules())
    assert molecules.length() == 0
    molecules.addSmiles(_test_smiles_1_)
    assert molecules.length() == 1
    molecules.addSmiles(_test_smiles_1_tup)
    assert molecules.length() == 2
    molecules.addSmiles(_test_smiles_error_)
    assert molecules.length() == 2
    molecules.addSmiles(_test_smiles_error_tup)
    assert molecules.length() == 2
    molecules.addSmiles([_test_smiles_2_, _test_smiles_2_])
    assert molecules.length() == 4
    molecules.addSmiles([_test_smiles_2_tup, _test_smiles_3_tup])
    assert molecules.length() == 6
    my_mols = molecules.getAllSmilesList()
    assert len(my_mols) == 6
    mol = molecules.getMol(2)
    #print(mol)
    assert mol["name"] == "mol_name_00002"
    mol = molecules.getMol(20)
    assert mol == None
    mol = molecules.getMol(-1)
    assert mol["name"] == "mol_3"
    mol = molecules.getMol(-20)
    assert mol == None
    #assert True == False
    logger.level = logging.ERROR

def test_read_write():
    """ test the read and write from files """
    molecules = Molecules.Molecules()
    assert type(molecules) == type(Molecules.Molecules())
    assert molecules.length() == 0
    molecules.addFromSmiFile(_null_file_)
    assert molecules.length() == 0
    molecules.addFromSmiFile(_test_smi_file_)
    assert molecules.length() == 10
    my_mols = molecules.getAllSmilesList()
    assert len(my_mols) == 10
    # TODO: Add the ability to write
    logger.level = logging.ERROR

def test_similarities():
    """ test the similarity functions """
    molecules = Molecules.Molecules()
    molecules.addFromSmiFile(_test_smi_file_)
    assert molecules.length() == 10
    #
    similarity_matrix = molecules.getIntraSimilarityMatrix()
    assert similarity_matrix['ADM_13086138']['ADM_13086138'] == 1.0
    assert similarity_matrix['ADM_13086138']['ADM_13142206'] == 0.494
    #
    intra_sim = molecules.getIntraSimilarities()
    assert len(intra_sim) == 10
    assert intra_sim[0] == 'ADM_13086138 ADM_13142440 0.502'
    assert intra_sim[6] == 'ADM_13142258 ADM_13142206 0.969'
    #print(pprint.pformat(intra_sim))
    molecules_2 = Molecules.Molecules()
    molecules_2.addFromSmiFile(_test_smi_2_file_)
    assert molecules.length() == 10
    inter_sim = molecules.getInterSimilarities(molecules_2)
    assert len(inter_sim) == 10
    assert inter_sim[0] == 'ADM_13086138 SYN_18547648 0.49'
    assert inter_sim[3] == 'ADM_13142223 SFA_21708021 0.449'
    #print(pprint.pformat(inter_sim))
    #assert True == False
    logger.level = logging.ERROR
