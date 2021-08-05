"""
get_ade_data.py
------------
Author: Daniel H Robertson

Part of the IBRI informatics system

compile the data needed for the ade tool
"""
import argparse
import os
import sys
import csv
import logging

import chem_utils
from chem_utils.data.db import db_chembl

ignore = ['gene_symbol']

mol_cols = [ 'query_id', 'molecule_chembl_id', 'molregno']

activity_cols = [
    'activity_id', 'assay_id', 'record_id', 'standard_relation',
    'standard_relation', 'standard_value', 'standard_unit', 'standard_flag', 'standard_type',
    'bao_endpoint', 'src_id', 'type', 'relation', 'value', 'units'
]

#
def read_file(in_file, id_column):
    #print(pprint.pformat(results))
    cpds = list()
    if not os.path.isfile(in_file):
        logger.critical("Unable to locate defined input file \"" + in_file +"\" for genes ... skipping")
        return []

    # read file
    nread = 0
    id_column_error = False
    with open(in_file) as file_ref:
# TODO: make delimiter an option - chembl uses ;
        csv_reader = csv.DictReader(file_ref, delimiter=';')
        for row in csv_reader:
            if id_column in row:
                cpds.append(row[id_column])
            else:
                if not id_column_error: # only print once
                    logger.error("Unable to find id_column \"{}\" in keys \"{}\"".format(
                        id_column, list(row.keys())))
                    id_column_error = True
            #print(row)
            nread += 1
    file_ref.close()

    logger.info("\n\nRead {} molecules from \"{}\"\n\n".format(len(cpds), in_file))

    return cpds

# return the chembl_id for the single protein --

#@smart_cache
# read a json mapping file
def get_molecule_activity_info(in_file, id_column, out_file):

    cpds = read_file(in_file, id_column)

    if not len(cpds):
        logger.critical("Did not read any genes from input file {} ... exiting".format(in_file))
        return

    # cache results
    #cpd_cache = get_cache()
    cpd_data = list()

    cpd2molregno = db_chembl.molregnos_from_chembl_ids(cpds)
    logger.debug("cpd2molregnos: {}".format(cpd2molregno))

    count = 0
    # do the pipeline for the genes
    for molecule_chembl_id in cpds:
        logger.debug("\n\nCPD: {}".format(molecule_chembl_id))
        count += 1
        if False and molecule_chembl_id in cpd_cache:
            cpd_info = cpd_cache[molecule_chembl_id]
        else:
            cpd_info = dict()
            cpd_info['query_id'] = molecule_chembl_id
            if True:
            #if molecule_chembl_id in cpd2molregno:
                molregno = cpd2molregno[molecule_chembl_id]
                logger.debug("molregno: {}".format(molregno))
                cpd_info['molregno'] = molregno
                results = db_chembl.activities_from_molregnos(molregno)
                cpd_info['activities'] = results
                logger.info('Found {} results for molecule_chembl_id: {}'.format(len(results), molecule_chembl_id))
                logger.debug("cpd: {} dict_results: {}".format(molecule_chembl_id, results))
            else:
                logger.critical("Unable to find moleregno for molecule_chembl_id: {} .. skipping".format(molecule_chembl_id))
            #at the risk of more I/O -- save cache
            #cpd_cache[molecule_chembl_id] = cpd_info
            #save_cache(cpd_cache)
        # add to growing list and
        cpd_data.append(cpd_info)


    #
    # either use out_file if defined or stdout otherwise
    #
    out_stream = sys.stdout
    if (out_file):
        out_stream = open(out_file, 'w', newline = '')
    csvwriter = csv.writer(out_stream, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

    #
    # write out the data from gene_data as csv
    #
    logger.debug('cpddata: {}'.format(cpd_data))
    headers = list()
    for col in mol_cols:
        headers.append(col)
    for col in activity_cols:
        headers.append(col)
    csvwriter.writerow(headers)
    for cpd_info in cpd_data:
        if len(cpd_info['activities']):
            for activity in cpd_info['activities']:
                new_row = list()
                # base info
                for col in mol_cols:
                    value = ""
                    if col in cpd_info:
                        value = cpd_info[col]
                    new_row.append(value)
                #activity info
                for act_col in activity_cols:
                    value = ""
                    if act_col in activity:
                        value = activity[act_col]
                    new_row.append(value)
                csvwriter.writerow(new_row)
        else:
            new_row = list()
            for col in mol_cols:
                value = ""
                if col in cpd_info:
                    value = cpd_info[col]
                new_row.append(value)
            for act_col in activity_cols:
                new_row.append("")
            csvwriter.writerow(new_row)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='get_ade_data.py',
        description='compile the ADE target information for a set of genes specified in an file.')
    parser.add_argument('-f', '--file', dest='in_file', help='csv file containing list of gene_symbols', required=True)
    parser.add_argument('-c', '--column', dest='id_column', help='column with chembl id for molecule', required=True)
    parser.add_argument('-o', '--output', dest='out_file', help='optional output file -- otherwise write to stdout')
    parser.add_argument('-v', '--verbose', dest='verbose', action="store_true", default=False, help='turn on addition information that is sent to stderr')
    parser.add_argument('-d', '--debug', dest='debug', action="store_true", default=False, help='turn on very verbose debug information to stderr')

    args = parser.parse_args()
    # create logger
    logger = logging.getLogger()
    logger.level = logging.WARN
    if args.verbose:
        logger.level = logging.INFO
    if args.debug:
        logger.level = logging.DEBUG
    stream_handler = logging.StreamHandler(sys.stderr)
    logger.addHandler(stream_handler)

    # process the out_file argument
    out_file = args.out_file if args.out_file else ""

    get_molecule_activity_info(args.in_file, args.id_column, out_file)
