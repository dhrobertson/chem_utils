"""
probeminer_tsv.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with probeminer data from downloaded probeminer tsv files

Command to download/sync files
    rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/21.04/output/etl/json/targets .

TODO: Need to define data directory and do the sync/download?

"""
import csv
import os
import json
import logging
from . import core_utils

# Global variables for connections so we don't create multiple
_data_cache_ = {
    'targets': dict(),
    'compounds': dict()
}

#
# define parquet file location
#
# TODO: Have this set in a json config file
config = {
    'dir' : '/Users/dhrobertson/data/probeminer'
}

# create logger
logger = logging.getLogger()
# return the rds connector for the database -- only create once per database

def _read_tsv_file(read_from_cache=True):

    if read_from_cache:
        data = _get_data_cache_()
        if len(data['targets']) and len(data['compounds']):
            return data

    data = {
        'targets': dict(),
        'compounds': dict()
    }
    if not os.path.isdir(config['dir']):
        logger.error("Unable to locate tsv directory \"{}\"".format(config['dir']))
        return data
    files = os.listdir(config['dir'])
    txt_files = [i for i in files if '.txt' in i]
    if len(txt_files) > 1:
        logger.critical("Exiting ... Found more than 1 file from dir \"{}\" .. confused. Files \"{}\": {}".format(dir, txt_files))
        exit()
    file = txt_files[0]
    (stem, ext) = file.split('.')
    date = file.split('_')[2]
    # Read from File
    logger.debug("  reading file \"{}\"".format(file))
    full_path = os.path.join(config['dir'], file)
    if not core_utils.file_exists(full_path):
        logger.critical("Exiting ... Issue with file {} from dir \"{}\" not existing".format(dir, file))
        exit()
    nrows = 0
    logger.info("Reading probeminer data from {}".format(file))
    with open(full_path, newline='') as fo:
        reader = csv.DictReader(fo, delimiter='\t')
        for row in reader:
            uniprot = row['UNIPROT_ACCESSION']
            #print("uniprot: {} {}".format(uniprot, type(uniprot)))
            if uniprot not in data['targets']:
                data['targets'][uniprot] = list()
            cpd_id  = row['COMPOUND_ID']
            #print("cpd_id: {} {}".format(cpd_id, type(cpd_id)))
            if cpd_id not in data['targets']:
                data['compounds'][cpd_id] = list()
            data['targets'][uniprot].append(row)
            data['compounds'][cpd_id].append(row)
            nrows += 1
    fo.close()
    logger.info("  read {} rows from {}".format(nrows, file))

    _set_data_cache_(data)

    return data

def _set_data_cache_(data):
    global _data_cache_

    _data_cache_ = data

def _get_data_cache_():
    global _data_cache_

    if _data_cache_:
        return _data_cache_

def get_target_info(uniprot_ids):
    """ get the details for a target -- input either single item or a list """

    logger.debug("parameter: {} type(parameter): {}".format(uniprot_ids, type(uniprot_ids)))
    if type(uniprot_ids) == type(str()):
        uniprot_ids = [uniprot_ids]

    data = _read_tsv_file()
    _targets_data = data['targets']
    records = list()
    # search by uniprot
    found = list()
    for id in uniprot_ids:
        if id in _targets_data:
            records.append(_targets_data[id])
            found.append(id)

    # return if found all
    if len(found) == len(uniprot_ids):
        return(records)

    logger.debug('  found {} items by submitted uniprot_ids {} '.format(len(found), found))
    # find others
    not_found = list()
    for id in uniprot_ids:
        if id not in found:
            not_found.append(id)
    logger.debug(' Did not find {} uniprot_ids {}'.format(len(not_found), not_found))

    return records
