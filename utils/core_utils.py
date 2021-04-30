"""
core_utils.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific general utility functions
"""
from pathlib import Path
import unittest
import os
import sys
import logging
import json
import csv
import pprint

# create logger
# TODO -- fix later -- check if logger exists -- and than attach to it
logger = logging.getLogger()

_assay_headers_ = [
     'molecule_chembl_id', 'molecule_pref_name', 'canonical_smiles',
     'target_chembl_id', 'target_pref_name',
     'activity_id', 'assay_chembl_id', 'assay_type', 'assay_description',
     'document_chembl_id','document_journal','document_year',
     'record_id', 'standard_type', 'standard_units',
     'type', 'units', 'value', 'upper_value'
     ]

# does a file exist
def file_exists(my_file):
    """ utility function to see if file exists """
    if not os.path.isfile(my_file):
        return False
    return True

# does a file exist
def remove_file(my_file):
    """ utility function to remove file """
    os.remove(my_file)

# read the content from a file
def read_text_from_file(in_file):
    """ utility function to read content from a file """
    if not os.path.isfile(in_file):
        logger.critical("Unable to locate defined input file \"" + in_file +"\" for content ... skipping")
        return None
    file_ref = open(in_file, "r")
    return file_ref.read()

# save content to a file
def save_text_to_file(out_file, content):
    """ utility function to save content to a file """
    file_ref = open(out_file, "w")
    file_ref.writelines(content)

# read a list from file -- assume first column
def read_csv_as_dict(in_file):
    """ utility function to read a list from a file """
    my_data = list()

    if not os.path.isfile(in_file):
        logger.critical("Unable to locate defined input csv file \"" + in_file +"\" ... skipping")
        return my_data

    with open(in_file) as file_ref:
        csv_reader = csv.DictReader(file_ref,delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for row in csv_reader:
            my_data.append(row)

    logger.info("Read " + str(len(my_data)) + " items from \"" + in_file + "\"")
    return my_data

# read a json from file
def read_json(in_file):
    """ utility function to read a json object from a file """
    my_data = list()

    if not os.path.isfile(in_file):
        logger.critical("Unable to locate defined input json file \"" + in_file +"\" ... skipping")
        return my_data

    with open(in_file) as file_ref:
        json_data = json.load(file_ref)

    logger.info("Read " + str(len(json_data)) + " items from \"" + in_file + "\"")
    return json_data

# read a json from file
def write_json(out_file, json_data):
    """ utility function to write a json object to file """
    try:
        fp = open(out_file, 'w')
        json.dump(json_data, fp, indent=2)
        fp.close()
    except:
        logger.critical("Error saving json_data to file \"" + out_file +"\" ... skipping")


# read a list from file -- assume first column
def read_list_from_file(in_file):
    """ utility function to read a list from a file """
    my_list = list()

    if not os.path.isfile(in_file):
        logger.critical("Unable to locate defined input file \"" + in_file +"\" for items ... skipping")
        return my_list

    with open(in_file) as file_ref:
        for line in file_ref.readlines():
            item = line.rstrip()
            if ',' in item:
                item = item.split(',')[0]
            if ' ' in item:
                item = item.split(' ')[0]
            if item:
                my_list.append(item)

    logger.info("Read " + str(len(my_list)) + " items from \"" + in_file + "\"")
    return my_list

def get_results_fields(result):
    result_fields = dict()
    for f in _assay_headers_:
        value = ''
        if f in result:
            value = result[f]
        result_fields[f] = value

    return result_fields

def write_assay_results(out_file, assay_results):
    """ write out the assay results as csv file """
    # convert to the object
    results = list()
    for item in assay_results:
        res_fields = get_results_fields(item)
        results.append(res_fields)
    csvwriter = get_csvwriter(out_file)
    write_csv(out_file, results, _assay_headers_)

def write_assay_results_json(out_file, assay_results):
    """ write out the assay results as csv file """
    # convert to the object
    fp = open(out_file, 'w')
    json.dump(assay_results, fp, indent=2)

# setup output stream
def get_csvwriter(out_file):
    #
    # either use out_file if defined or stdout otherwise
    #
    out_stream = sys.stdout
    if (out_file):
        out_stream = open(out_file, 'w', newline = '')

    #print(pp.pformat(out_stream))
    csvwriter = csv.writer(out_stream, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

    return csvwriter

# function to write out the csv files
def write_csv(out_file, item_list, headers):
    no_lines = 0
    csvwriter = get_csvwriter(out_file)
    csvwriter.writerow(headers)
    for item in item_list:
        #print(pprint.pformat(tgt))
        row = list()
        for field in headers:
            if field in item:
                row.append(item[field])
            else:
                row.append('')
        csvwriter.writerow(row)
        no_lines += 1
    return no_lines
