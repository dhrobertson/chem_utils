"""
api_chembl.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with chembl data directly through the chembl api
"""
import psycopg2
import json
import logging

# Global variables for connections so we don't create multiple
_postgresql_connector = dict()
_databases = ['chembl_29', 'chembl_28', 'oprm1_chembl_29']
_selected_database = _databases[0]

# create logger and set logger level
logger = logging.getLogger()

## key methods to do the db requests for chembl

def get_selected_database():
    return _selected_database

#TODO: Need to make this dynamic
def set_database(database):
    global _selected_database
    if database not in _databases:
        logger.error("Unable to set database to \"{}\" as it is not one of the supported databases")
        return False
    _selected_database = database
    return True

# return the rds connector for the database -- only create once per database
def _get_db_connector():
    global _postgresql_connector
    global _selected_database

    if _selected_database not in _postgresql_connector:
        _postgresql_connector[_selected_database] = dict()

    if 'conn' not in _postgresql_connector[_selected_database] or _postgresql_connector[_selected_database]['conn'] == None:
        # go ahead and try to connect
        try:
            _postgresql_connector[_selected_database]['conn'] = psycopg2.connect(dbname=_selected_database)
        except Exception as e:
            logger.error('***DB Connection Create Error for {}: {}'.format(_selected_database, e))
            _postgresql_connector[_selected_database]['conn'] = None

    return _postgresql_connector[_selected_database]['conn']

# return the rds connector for the database -- only create once per database
def _get_db_cursor():
    global _postgresql_connector
    global _selected_database

    if _selected_database not in _postgresql_connector:
        _postgresql_connector[_selected_database] = dict()

    if 'cur' not in _postgresql_connector[_selected_database] or _postgresql_connector[_selected_database]['cur'] == None:
        # go ahead and try to connect
        try:
            conn = _get_db_connector()
            _postgresql_connector[_selected_database]['cur'] = conn.cursor()
        except Exception as e:
            logger.error('***DB Connection Cursor Error: {}'.format(e))
            _postgresql_connector[_selected_database]['cur'] = None

    return _postgresql_connector[_selected_database]['cur']

def _db_execute(sql_statement):
    cur = _get_db_cursor()
    data = None
    try:
        cur.execute(sql_statement)
        data = cur.fetchall()
    except Exception as e:
        logger.error('***DB Execute Error: {}'.format(e))
        conn = _get_db_connector()
        conn.rollback()
    # post process the data
    return data

# utility functions
def _ids_to_str_list(ids):
    ids_list = ids
    if type(ids_list) != type(list()):
        ids_list = [ids]

    ids_list_str = ','.join(f"'{id}'" for id in ids_list)

    return ids_list_str

def _structured_to_json(data):
    results_json = list()
    if 'fields' not in data or 'values' not in data:
        logger.error('does not looks like structured data')
        return results_json
    fields = data['fields']
    for row in data['values']:
        if len(row) == len(fields):
            new_row = dict()
            for i in range(0, len(row)):
                new_row[fields[i]] = row[i]
            results_json.append(new_row)
        else:
            logger.error('issue with len(fields) != len(rows) ... skipping')

    logger.debug("results_json: {}".format(results_json))
    return results_json

def _db_execute_structured(sql_statement, table, limit=None):

    data = {
        'fields' : list(),
        'values' : list()
    }

    # get columns from table and add to data object
    sql = "SELECT column_name FROM information_schema.columns WHERE TABLE_NAME='{}'".format(table)
    db_results = _db_execute(sql)
    for result in db_results:
        data['fields'].append(result[0])

    # do data query and put values in place
    sql = sql_statement
    if limit:
        sql = "{} limit {}".format(sql, limit)
    sql_results = _db_execute(sql)
    if sql_results:
        for row in sql_results:
            new_row = list()
            for value in row:
                new_row.append(value)
            data['values'].append(new_row)

    return data

###############################################################################
## individual methods                                                        ##
## TODO: refactor to have wrapper functions depending on whether api or db   ##
###############################################################################

def get_all_tids():
    """ get all tids from target_dictionary"""
    sql_results = _db_execute('select tid from target_dictionary')
    logger.debug("results: {}".format(sql_results))
    ids_list = list()
    for item in sql_results:
        ids_list.append(item[0])
    logger.debug("ids_list: {}".format(ids_list))

    return ids_list


def get_all_molregnos():
    """ get all molregno from molecule_dictionary"""
    sql_results = _db_execute('select molregno from molecule_dictionary')
    logger.debug("results: {}".format(sql_results))
    ids_list = list()
    for item in sql_results:
        ids_list.append(item[0])
    logger.debug("ids_list: {}".format(ids_list))

    return ids_list

def tids_from_chembl_ids(chembl_ids):
    """ search for tids by chembl_target_id in target_dictionary"""
    return _get_id_mapping(chembl_ids, 'select chembl_id,tid from target_dictionary where chembl_id')

def chembl_ids_from_tids(tids):
    """ search for chembl_ids by tids in target_dictionary """
    return _get_id_mapping(tids, 'select tid,chembl_id from target_dictionary where tid')

def molregnos_from_chembl_ids(chembl_ids):
    """ search for molregnos by chembl_molecule_id in molecule_dictionary"""
    return _get_id_mapping(chembl_ids, 'select chembl_id,molregno from molecule_dictionary where chembl_id')

def chembl_ids_from_molregnos(molregnos):
    """ search for molregnos by chembl_molecule_id in molecule_dictionary"""
    return _get_id_mapping(molregnos, 'select molregno,chembl_id from molecule_dictionary where molregno')

def assay_ids_from_chembl_ids(chembl_ids):
    """ search for assay_id by chembl_assay_id in assays"""
    return _get_id_mapping(chembl_ids, 'select chembl_id,assay_id from assays where chembl_id')

def assay_ids_from_tids(tids):
    """ search for assay_id by tid in assays"""
    results = _get_id_mapping(tids, 'select tid,assay_id from assays where tid')
    logger.debug("results: {}".format(results))

    # handle return value type (remove from dict if single value)
    logger.debug("type: {} len: {}".format(type(tids), len(results)))
    if type(tids) != type(list()) and tids in results:
        return results[tids]
    return results

def chembl_ids_from_assay_ids(assay_ids):
    """ search for assay_ids by chembl_molecule_id in molecule_dictionary"""
    return _get_id_mapping(assay_ids, 'select assay_id,chembl_id from assays where assay_id')

def assay_details_from_assay_ids(assay_ids):
    """ get assay details assay_id from assays"""
    ids_list_str = _ids_to_str_list(assay_ids)
    sql_statement = "{} in ({})".format('select * from assays where assay_id', ids_list_str)
    sql_results = _db_execute_structured(sql_statement, 'assays')
    logger.debug("results: {}".format(sql_results))
    json_results = _structured_to_json(sql_results)

    # change return based on receiving list or single value
    logger.debug("type: {} len: {}".format(type(assay_ids), len(json_results)))
    if type(assay_ids) != type(list()) and len(json_results) == 1:
        return json_results[0]

    return json_results

def get_all_assay_ids():
    """ get all assay_ids from assays"""
    sql_results = _db_execute('select assay_id from assays')
    logger.debug("results: {}".format(sql_results))
    ids_list = list()
    for item in sql_results:
        ids_list.append(item[0])
    logger.debug("ids_list: {}".format(ids_list))

    return ids_list

def activities_from_assay_ids(assay_ids):
    """ retrieve all activities for assay_ids from activities"""
    ids_list_str = _ids_to_str_list(assay_ids)
    sql_statement = "{} in ({})".format('select * from activities where assay_id', ids_list_str)
    sql_results = _db_execute_structured(sql_statement, 'activities')
    logger.debug("results: {}".format(sql_results))
    json_results = _structured_to_json(sql_results)

    logger.debug("json results: {}".format(json_results))
    return json_results

def activities_from_molregnos(molregnos):
    """ retrieve all activities for assay_ids from activities"""
    ids_list_str = _ids_to_str_list(molregnos)
    sql_statement = "{} in ({})".format('select * from activities where molregno', ids_list_str)
    sql_results = _db_execute_structured(sql_statement, 'activities')
    logger.debug("results: {}".format(sql_results))
    json_results = _structured_to_json(sql_results)

    logger.debug("json results: {}".format(json_results))
    return json_results

def get_bioactivities_for_targets(chembl_ids):
    """ retrieve all activities for target chembl_ids from activities"""

# TODO: DHR:s Fix and Test
#def tids_from_component_ids(component_ids):
#    """ search for tids by component_id in target_components"""
#    return _get_id_mapping(component_ids, 'select component_id,tid from target_components where component_id')

# TODO: DHR: Fix and Test
#def component_ids_from_tids(tids):
#    """ search for component_ids by tids in target_components """
#    return _get_id_mapping(tids, 'select tid,component_id from target_components where tid')


# generalized function assuming two values returned generate a dictionary that has one in front of the other
def _get_id_mapping(ids, sql_statement):
    ids_list_str = _ids_to_str_list(ids)
    sql_statement = "{} in ({})".format(sql_statement, ids_list_str)
    sql_results = _db_execute(sql_statement)
    logger.debug("sql results: {}".format(sql_results))
    results = dict()
    if sql_results:
        for result in sql_results:
            (key, value) = result
            # make array if already in dictionary
            if key in results:
                if type(results[key]) != type(list()):
                    prior_value = results[key]
                    results[key] = [prior_value]
                results[key].append(value)
            else:
                results[key] = value
    logger.debug("results: {}".format(results))

    return results
