"""
get_cpd_tools_by_accession.py
------------
Author: Daniel H Robertson

Part of the IBRI informatics system

compile the data for potential tool compounds from CheMBL and probeminer
"""
import argparse
import os
import json
import sys
import csv
import pprint
import logging
import math

import chem_utils
from chem_utils.data import probeminer_tsv
from chem_utils.data import chembl

ignore = ['gene_symbol']

csv_map = [
    { 'field' : 'GENE_SYMBOL', 'label': 'GENE_SYMBOL'},
    { 'field' : 'UNIPROT_ACCESSION', 'label': 'UNIPROT_ACCESSION'},
    { 'field' : 'generic_name', 'label': 'GENERIC_NAME'},
    { 'field' : 'IBRI_REG_NO', 'label': 'IBRI_REG_NO'},
    { 'field' : 'ChEMBL', 'label': 'CHEMBL_ID'},
    { 'field' : 'target_potency', 'label': 'TARGET_POTENCY'},
    { 'field' : 'off_targets', 'label': 'OFF_TARGETS_POTENCY'},
    { 'field' : 'pubmed_ids', 'label' : 'PM_PUBMED_IDS'},
    { 'field' : 'smiles', 'label': 'SMILES'},
    { 'field' : 'scaffold', 'label': 'PM_SCAFFOLD'},
    { 'field' : 'mw', 'label': 'MW'},
    { 'field' : 'clogp', 'label': 'clogP'},
    { 'field' : 'tPSA', 'label': 'tPSA'},
    { 'field' : 'HBA', 'label': 'HBD'},
    { 'field' : 'pains', 'label' : 'PM_PAINS'},
    { 'field' : 'pains_text', 'label' : 'PM_PAINS_TEXT'},
    { 'field' : 'indications', 'label' : 'INDICATIONS'},
    { 'field' : 'COMPOUND_ID', 'label': 'PM_COMPOUND_ID'},
    { 'field' : 'no_secondary_targets', 'label' : 'PM_NUM_SEC_TARGETS'},
    { 'field' : 'selectivity_comp2', 'label' : 'PM_SELECTIVITY'},
]
# TODO -- make script to find this from HUGO
# set this mapping in an easy utility to grab
uniprot_map = {
    'Q13153' : 'PAK1',
    'Q13177' : 'PAK2',
    'O60885' : 'BRD4'
}

# cache results
cache_file = 'cache.json'
def get_cache():
    cache = dict()
    if not os.path.isfile(cache_file):
        logger.info("Unable to locate cache file \"" + cache_file +"\" for genes ... skipping")
        return cache

    with open(cache_file) as fp:
        cache = json.load(fp)
    fp.close()

    return cache

def save_cache(cache):
    with open(cache_file, "w") as fp:
        json.dump(cache, fp, indent=2)
    fp.close()

#
def read_file(in_file):
    #print(pprint.pformat(results))
    ids = list()
    extra_data = {
        'fields' : list(),
        'data'   : dict()
    }
    if not os.path.isfile(in_file):
        logger.critical("Unable to locate defined input file \"" + in_file +"\" for ids ... skipping")
        return (ids, extra_data)

    # read in as csv
    with open(in_file, newline='', encoding='utf-8') as fp:
        n = 0
        csvreader = csv.reader(fp, delimiter=',')
        for row in csvreader:
            if n == 0:
                headers = row
                for i in range(0, len(headers)):
                    extra_data['fields'].append(row[i])
            else:
                id = row[0]
                id = id.replace(" ", "")
                if id not in extra_data['data']:
                    extra_data['data'][id] = list()
                if id not in ids: # in case of duplicates
                    ids.append(id)
                data = dict()
                for i in range(0, len(headers)):
                    data[headers[i]] = row[i]
                extra_data['data'][id].append(data)
            n += 1
        logger.info("Read " + str(n) + " lines from file \"" + str(in_file) + "\"" )
    fp.close()

    logger.debug("ids: {} extra_data: {}".format(ids, extra_data))

    return (ids, extra_data)


def get_tool_info(in_file, out_file):

    (uniprot_ids, extra_data) = read_file(in_file)

    if not len(uniprot_ids):
        logger.critical("Did not read any uniprot_ids from input file {} ... exiting".format(in_file))
        return

    logger.debug("ids: {} \n extra_data: {}".format(uniprot_ids, extra_data))

    results = list()
    not_found = list()
    data_cache = get_cache()
    for id in uniprot_ids:
        if id in data_cache:
            results.append(data_cache[id])
        else:
            not_found.append(id)

    # get target results if not already in cache
    if len(not_found):
        results_new = probeminer_tsv.get_target_info(not_found)
        for i in range(0, len(results_new)):
            uniprot_id = results_new[i][0]['UNIPROT_ACCESSION']
            data_cache[uniprot_id] = results_new[i]
        save_cache(data_cache)

    if len(results) != len(uniprot_ids):
        logger.info("Did not return information for all ids: {} sent {} information returned".format(len(uniprot_ids), len(results)))

    cpd_data = list()
    for i in range(0, len(results)):
        #print(results[i][0])
        uniprot_id = results[i][0]['UNIPROT_ACCESSION']
        logger.info("For id {} returned {} records".format(uniprot_id, len(results[i])))
        for cpd in results[i]:
            cpd_row = dict()
            uniprot_id = results[i][0]['UNIPROT_ACCESSION']
            if uniprot_id in uniprot_map:
                cpd_row['GENE_SYMBOL'] = uniprot_map[uniprot_id]
            for k in cpd.keys():
                cpd_row[k] = cpd[k]

            # unwrap xrefs
            xrefs_json = json.loads(cpd['xrefs'])
            logger.debug("xrefs: {} {} {}".format(cpd['xrefs'], type(cpd['xrefs']), type(xrefs_json)))
            for ref in ['ChEMBL', 'BindingDB', 'canSAR']:
                if ref in xrefs_json:
                    cpd_row[ref] = xrefs_json[ref]
            # unwrap pubmed_ids
            pubmed_json = json.loads(cpd['pubmed_ids'])
            cpd_row['pubmed_ids'] = ' '.join(pubmed_json)

            # add row to list
            cpd_data.append(cpd_row)
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
    logger.debug('cpd_data: {}'.format(cpd_data))
    headers = list()
    fields = list()
    for item in csv_map:
        headers.append(item['label'])
        fields.append(item['field'])
    csvwriter.writerow(headers)
    for cpd in cpd_data:
        new_row = list()
        for f in fields:
            value = ""
            if f in cpd:
                value = cpd[f]
            new_row.append(value)
        csvwriter.writerow(new_row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='get_ade_data.py',
        description='compile the ADE target information for a set of genes specified in an file.')
    parser.add_argument('-f', '--file', dest='in_file', help='txt file containing list of gene_symbols', required=True)
    parser.add_argument('-o', '--output', dest='out_file', help='optional output file -- otherwise write to stdout')
    parser.add_argument('-v', '--verbose', action="store_true", default=False,help='turn on addition information that is sent to stderr')
    parser.add_argument('-d', '--debug', action="store_true", default=False,help='turn on debug information that is sent to stderr')

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

    get_tool_info(args.in_file, out_file)
