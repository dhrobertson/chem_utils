"""
test_api_chembl.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing chembl information
"""
from data.db import db_chembl
import json
import logging

# create logger
logger = logging.getLogger()
logger.level = logging.DEBUG

# one to one
_tgt_chembl_id_to_tid = {
    'CHEMBL1781870' : 104162,
    'CHEMBL2331064' : 105660
}
# set up for reverse : one to one
_tgt_tid_to_chembl_id = dict()
for id in _tgt_chembl_id_to_tid.keys():
    _tgt_tid_to_chembl_id[_tgt_chembl_id_to_tid[id]] = id

_cpd_chembl_id_to_molregno = {
    'CHEMBL1974235' : 1295526,
    'CHEMBL3629617' : 1969582
}
# set up for reverse : one to one
_cpd_molregno_to_chembl_id = dict()
for id in _cpd_chembl_id_to_molregno.keys():
    _cpd_molregno_to_chembl_id[_cpd_chembl_id_to_molregno[id]] = id

_assay_chembl_id_to_assay_id = {
    'CHEMBL1785363' : 751726,
    'CHEMBL4274857' : 1802565
}
# set up for reverse : one to one
_assay_id_to_chembl_id = dict()
for id in _assay_chembl_id_to_assay_id.keys():
    _assay_id_to_chembl_id[_assay_chembl_id_to_assay_id[id]] = id

# one to one | reverse: many to one
_tgt_component_id_to_tids = {
    261 : [169, 104801, 105734, 117292, 117734],
    2225: [11842, 105769]
}
_test_molecule_id = 'CHEMBL3629617'
_test_molecule_id_list = ['CHEMBL3629617','CHEMBL1974235']
_test_molregno_list = [1626744, 435293]
_test_assay_id = 'CHEMBL4274857'
_test_assay_id_list = ['CHEMBL4274857','CHEMBL1785363']
_test_assay_id_list_lg = ['CHEMBL904416', 'CHEMBL904417', 'CHEMBL939816', 'CHEMBL991464', 'CHEMBL991465',
                          'CHEMBL1052162', 'CHEMBL1052166', 'CHEMBL1226577', 'CHEMBL1226583', 'CHEMBL1226716',
                          'CHEMBL1249424', 'CHEMBL1249425', 'CHEMBL1249426', 'CHEMBL1249427', 'CHEMBL1249428' ]
_test_accession_id = 'Q13936'
_test_accession_id_list = ['Q13936', 'P07948']

# test the interface for simple calls
def test_db_connection():
    _test_name = 'test_db_connection'
    # check get connector
    db_conn = db_chembl._get_db_connector()
    assert db_conn != None

    # check get same connector second time
    db_conn_2 = db_chembl._get_db_connector()
    assert db_conn == db_conn_2

    # make sure cursor is defined
    db_cur = db_chembl._get_db_cursor()
    assert db_conn != None

    # check connections and tables
    version = db_chembl._db_execute('SELECT version();')
    assert len(version) == 1
    tables = db_chembl._db_execute('SELECT * FROM information_schema.tables;')
    assert len(tables) == 273

    # check connections and tables with limit
    data = db_chembl._db_execute_structured('SELECT * from target_dictionary', 'target_dictionary', 10)
    assert 'fields' in data
    assert len(data['fields']) == 7
    assert 'values' in data
    assert len(data['values']) == 10

    # check connections and tables without limit
    data = db_chembl._db_execute_structured('SELECT * from target_dictionary', 'target_dictionary')
    assert 'fields' in data
    assert len(data['fields']) == 7
    assert 'values' in data
    assert len(data['values']) >= 14347

def test_utility_functions():
    """ test the utility functions """
    ids_str = db_chembl._ids_to_str_list(['0','1','2'])
    assert ids_str == "'0','1','2'"

    _test_structured = {
        'fields' : ['field_1', 'field_2'],
        'values' : [
            [ 0, 1 ],
            [ 2, 3 ]
        ]
    }
    my_json = db_chembl._structured_to_json(_test_structured)
    assert len(my_json) == len(_test_structured['values'])
    for i in range(0, len(_test_structured['values'])):
        for j in range(0, len(_test_structured['fields'])):
            assert my_json[i][_test_structured['fields'][j]] == _test_structured['values'][i][j]

    # Force debug listing
    #assert True == False

def test_molecule_queries():
    """ test the core queries to get molecule information from various tables """
    _test_name = 'test_molecule_queries'

    # chembl_id to molregno and vice versa

    _cpd_id_list = list(_cpd_chembl_id_to_molregno.keys())
    molregnos = db_chembl.molregnos_from_chembl_ids(_cpd_id_list[0])
    assert type(molregnos) == type(dict())
    assert len(molregnos) == 1
    for id in molregnos.keys():
        assert molregnos[id] == _cpd_chembl_id_to_molregno[id]

    molregnos = db_chembl.molregnos_from_chembl_ids(_cpd_id_list)
    assert type(molregnos) == type(dict())
    assert len(molregnos) == 2
    for id in molregnos.keys():
        assert molregnos[id] == _cpd_chembl_id_to_molregno[id]

    _cpd_id_list = list(_cpd_molregno_to_chembl_id.keys())
    chembl_ids = db_chembl.chembl_ids_from_molregnos(_cpd_id_list[0])
    assert type(chembl_ids) == type(dict())
    assert len(chembl_ids) == 1
    for id in chembl_ids.keys():
        assert chembl_ids[id] == _cpd_molregno_to_chembl_id[id]

    chembl_ids = db_chembl.chembl_ids_from_molregnos(_cpd_id_list)
    assert type(chembl_ids) == type(dict())
    assert len(chembl_ids) == 2
    for id in chembl_ids.keys():
        assert chembl_ids[id] == _cpd_molregno_to_chembl_id[id]

    all_molregnos = db_chembl.get_all_molregnos()
    assert type(all_molregnos) == type(list())
    assert len(all_molregnos) > 2105460

    # Force debug listing
    #assert True == False

def test_target_queries():
    """ test the core queries to get target information from various tables """
    _test_name = 'test_target_queries'

    # chembl_id to tid and vice versa

    _tgt_id_list = list(_tgt_chembl_id_to_tid.keys())
    tids = db_chembl.tids_from_chembl_ids(_tgt_id_list[0])
    assert type(tids) == type(dict())
    assert len(tids) == 1
    for id in tids.keys():
        assert tids[id] == _tgt_chembl_id_to_tid[id]

    tids = db_chembl.tids_from_chembl_ids(_tgt_id_list)
    assert type(tids) == type(dict())
    assert len(tids) == 2
    for id in tids.keys():
        assert tids[id] == _tgt_chembl_id_to_tid[id]

    _tgt_id_list = list(_tgt_tid_to_chembl_id.keys())
    chembl_ids = db_chembl.chembl_ids_from_tids(_tgt_id_list[0])
    assert type(chembl_ids) == type(dict())
    assert len(chembl_ids) == 1
    for id in chembl_ids.keys():
        assert chembl_ids[id] == _tgt_tid_to_chembl_id[id]

    chembl_ids = db_chembl.chembl_ids_from_tids(_tgt_id_list)
    assert type(chembl_ids) == type(dict())
    assert len(chembl_ids) == 2
    for id in chembl_ids.keys():
        assert chembl_ids[id] == _tgt_tid_to_chembl_id[id]

    all_tids = db_chembl.get_all_tids()
    assert type(all_tids) == type(list())
    assert len(all_tids) > 14550
    # assay_ids from tids

    # component_id to tid and vice versa
    #_tgt_id_list = list(_tgt_tid_to_component_id.keys())
    #tids = db_chembl.component_ids_from_tids(_tgt_id_list[0])
    #assert type(tids) == type(dict())
    #assert len(tids) == 1
    #for id in tids.keys():
    #    assert tids[id] == _tgt_tid_to_component_id[id]

    #tids = db_chembl.component_ids_from_tids(_tgt_id_list)
    #assert type(tids) == type(dict())
    #assert len(tids) == 7
    #for id in tids.keys():
    #    assert tids[id] == _tgt_tid_to_component_id[id]

    # Force debug listing
    #assert True == False

def test_assay_queries():
    """ test the core queries to get assay information from various tables """
    _test_name = 'test_assay_queries'

    # chembl_id to assay_id and vice versa

    _assay_id_list = list(_assay_chembl_id_to_assay_id.keys())
    assay_ids = db_chembl.assay_ids_from_chembl_ids(_assay_id_list[0])
    assert type(assay_ids) == type(dict())
    assert len(assay_ids) == 1
    for id in assay_ids.keys():
        assert assay_ids[id] == _assay_chembl_id_to_assay_id[id]

    assay_ids = db_chembl.assay_ids_from_chembl_ids(_assay_id_list)
    assert type(assay_ids) == type(dict())
    assert len(assay_ids) == 2
    for id in assay_ids.keys():
        assert assay_ids[id] == _assay_chembl_id_to_assay_id[id]

    _assay_id_list = list(_assay_id_to_chembl_id.keys())
    chembl_ids = db_chembl.chembl_ids_from_assay_ids(_assay_id_list[0])
    assert type(chembl_ids) == type(dict())
    assert len(chembl_ids) == 1
    for id in chembl_ids.keys():
        assert chembl_ids[id] == _assay_id_to_chembl_id[id]

    chembl_ids = db_chembl.chembl_ids_from_assay_ids(_assay_id_list)
    assert type(chembl_ids) == type(dict())
    assert len(chembl_ids) == 2
    for id in chembl_ids.keys():
        assert chembl_ids[id] == _assay_id_to_chembl_id[id]

    # assays for tids

    _tgt_id_list = list(_tgt_tid_to_chembl_id.keys())
    assay_ids = db_chembl.assay_ids_from_tids(_tgt_id_list[0])
    assert type(assay_ids) == type(list())
    assert len(assay_ids) == 3
    for id in [751726, 1524224, 1524225]:
        assert id in assay_ids

    _tgt_id_list = list(_tgt_tid_to_chembl_id.keys())
    assay_ids = db_chembl.assay_ids_from_tids(_tgt_id_list)
    assert type(assay_ids) == type(dict())
    assert len(assay_ids) == len(_tgt_id_list)
    assert _tgt_id_list[0] in assay_ids
    assert len(assay_ids[_tgt_id_list[0]]) == 3
    for id in _tgt_id_list:
        #assert id in assay_ids:
        assert type(assay_ids[id]) == type(list())
        assert len(assay_ids[id]) > 0

    # assay details
    _assay_id_list = list(_assay_id_to_chembl_id.keys())
    assay_details = db_chembl.assay_details_from_assay_ids(_assay_id_list[0])
    assert type(assay_details) == type(dict())
    assert 'assay_id' in assay_details
    assert assay_details['assay_id'] == _assay_id_list[0]
    assert 'assay_type' in assay_details
    assert assay_details['assay_type'] == 'B'

    assay_details = db_chembl.assay_details_from_assay_ids(_assay_id_list)
    assert type(assay_details) == type(list())
    assert len(assay_details) == 2

    all_assay_ids = db_chembl.get_all_assay_ids()
    assert type(all_assay_ids) == type(list())
    assert len(all_assay_ids) > 1383550
    # Force debug listing
    #assert True == False

def test_activities_queries():
    """ test the core queries to get activity information from various tables """
    _test_name = 'test_activity_queries'

    # get activities for assay_id
    _assay_id_list = list(_assay_id_to_chembl_id.keys())
    activities = db_chembl.activities_from_assay_ids(_assay_id_list[0])
    assert type(activities) == type(list())
    assert len(activities) == 3

    activities = db_chembl.activities_from_assay_ids(_assay_id_list)
    assert type(activities) == type(list())
    assert len(activities) == 54

    # get activities for target chembl_id
    # get activities for assay_id
    _assay_id_list = list(_assay_id_to_chembl_id.keys())
    activities = db_chembl.activities_from_molregnos(_test_molregno_list[0])
    assert type(activities) == type(list())
    assert len(activities) == 8

    activities = db_chembl.activities_from_molregnos(_test_molregno_list)
    assert type(activities) == type(list())
    assert len(activities) == 47

    # Force debug listing
    #assert True == False
