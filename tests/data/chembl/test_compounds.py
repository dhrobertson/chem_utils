"""
test_compounds.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing chembl information
"""
import pprint
import logging
from pathlib import Path
from data.chembl import compounds
from utils import core_utils

_test_sdf_read_file_ = str(Path(__file__).parent / 'data/cpds_test_lg.txt')

logger = logging.getLogger()

def test_get_by_targets():
    # 2021-02-03
    # CHEMBL1781870 (SHIP1) -  16 compounds
    # CHEMBL2331064 (SHIP2) -  61 compounds
    # SHIP1 & SHIP2 - 64
    """ test getting the active compounds from target_chembl_id """

    cpds = compounds.get_by_targets(['CHEMBL1781870'])
    assert len(cpds) == 16
    cpds_2 = compounds.get_by_targets(['CHEMBL1781870','CHEMBL2331064'])
    assert len(cpds_2) == 71

def test_get_details():
    # 2021-02-03
    # CHEMBL207881 (GW-9508) - 90 bioactivities
    # CHEMBL203297 () - 2 bioactivities
    """ test getting the active compounds from molecule_chembl_id """

    cpd_error = compounds.get_details('CHEMBL1781870')
    assert type(cpd_error) == type(list())
    assert len(cpd_error) == 0

    cpd = compounds.get_details('CHEMBL207881')
    assert type(cpd) == type(list())
    assert len(cpd) == 1
    assert 'molecule_chembl_id' in cpd[0]
    assert cpd[0]['molecule_chembl_id'] == 'CHEMBL207881'

    cpd_list = compounds.get_details(['CHEMBL207881'])
    assert type(cpd_list) == type(list())
    assert len(cpd_list) == 1
    assert 'molecule_chembl_id' in cpd_list[0]
    assert cpd_list[0]['molecule_chembl_id'] == 'CHEMBL207881'

    cpds_list = compounds.get_details(['CHEMBL207881','CHEMBL203297'])
    assert type(cpds_list) == type(list())
    assert len(cpds_list) == 2


def test_get_sdf():
    """ test the creating of an sdf from chembl_ids """

    sdf_1 = compounds.get_sdf(['CHEMBL207881'])
    assert 'CHEMBL207881' in sdf_1

    sdf_2 = compounds.get_sdf(['CHEMBL207881','CHEMBL203297'])
    assert 'CHEMBL207881' in sdf_2
    assert 'CHEMBL203297' in sdf_2
    #print(pprint.pformat(sdf_2))

    cpd_list = core_utils.read_list_from_file(_test_sdf_read_file_)
    assert len(cpd_list) == 200
    sdf_3 = compounds.get_sdf(cpd_list)
    assert 'CHEMBL207881' in sdf_3
    assert 'CHEMBL483546' in sdf_3
    #print(pprint.pformat(sdf_3))

def test_get_smi():
    """ test the creating of an sdf from chembl_ids """

    smi_1 = compounds.get_smi(['CHEMBL207881'])
    assert 'CHEMBL207881' in smi_1
    print(pprint.pformat(smi_1))

    sdf_2 = compounds.get_smi(['CHEMBL207881','CHEMBL203297'])
    assert 'CHEMBL207881' in sdf_2
    assert 'CHEMBL203297' in sdf_2
    #print(pprint.pformat(sdf_2))

    cpd_list = core_utils.read_list_from_file(_test_sdf_read_file_)
    assert len(cpd_list) == 200
    sdf_3 = compounds.get_smi(cpd_list)
    assert 'CHEMBL207881' in sdf_3
    assert 'CHEMBL483546' in sdf_3
    #print(pprint.pformat(sdf_3))
