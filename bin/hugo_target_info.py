"""
get_ade_data.py
------------
Author: Daniel H Robertson

Part of the IBRI informatics system

compile the data needed for the ade tool
"""
import argparse
import os
import json
import sys
import csv
import pprint
import logging
import copy

import chem_utils
from chem_utils.data import hugo_api
from chem_utils.data import pharos_api
from chem_utils.data import opentarget_json
from chem_utils.data import chembl as chembl_api
from chem_utils.data.db import db_chembl

ignore = ['gene_symbol']

csv_map = [
    { 'field' : 'query_symbol', 'label': 'ORIG_GENE_SYMBOL'},
    { 'field' : 'hgnc_symbol', 'label': 'HGNC_GENE_SYMBOL'},
    { 'field' : 'hgnc_name', 'label': 'HGNC_NAME'},
    { 'field' : 'hgnc_entrez_id', 'label': 'HGNC_ENTREZ_ID'},
    { 'field' : 'hgnc_ensembl_gene_id', 'label': 'HGNC_ENSEMBLE_GENE_ID'},
    { 'field' : 'hgnc_uniprot_ids', 'label': 'HGNC_UNIPROT_IDS'},
    { 'field' : 'hgnc_refseq_accession', 'label': 'HGNC_REFSEQ_ACCESSION'}
]

data_map = [
    { 'source_field': 'symbol', 'target_field': 'hgnc_symbol' },
    { 'source_field': 'hgnc_id', 'target_field': 'hgnc_id' },
    { 'source_field': 'name', 'target_field': 'hgnc_name' },
    { 'source_field': 'entrez_id', 'target_field': 'hgnc_entrez_id' },
    { 'source_field': 'ensembl_gene_id', 'target_field': 'hgnc_ensembl_gene_id' },
    { 'source_field': 'uniprot_ids', 'target_field': 'hgnc_uniprot_ids' },
    { 'source_field': 'refseq_accession', 'target_field': 'hgnc_refseq_accession'}
]

#
def read_file(in_file):
    #print(pprint.pformat(results))
    genes = list()
    if not os.path.isfile(in_file):
        logger.critical("Unable to locate defined input file \"" + in_file +"\" for genes ... skipping")
        return []

    # read file
    nread = 0
    with open(in_file) as file_ref:
        for line in file_ref.readlines():
            gene_symbol = line.rstrip()
            if gene_symbol not in ignore:
                genes.append(gene_symbol)
                nread += 1
    file_ref.close()

    logger.info("\n\nRead " + str(nread) + " genes from \"" + in_file + "\"\n\n")

    return genes

# get the single protein assession from list
def get_single_protein(targets):
    id_list = list()
    for target in targets:
        if 'target_components' in target and len(target['target_components']) == 1:
            id_list.append(target['target_chembl_id'])
    if not len(id_list):
        logger.critical("ChEMBL: Unable to get target_chembl_id from response")
        return ''
    if len(id_list) > 1:
        logger.error("ChEMBL: Found multiple single component -- need to check -- returning only first".format(id_list))
    return id_list[0]

def pharos_process_counts(gene_info, data):
    #logger.debug("data: {}".format(data))
    if 'ligandCounts' in data:
        for item in data['ligandCounts']:
            logger.debug("item: {}".format(item))
            field = 'pharos_ligand_drug'
            if item['name'] == 'ligand':
                field = 'pharos_ligand_active'
            gene_info[field] = item['value']
    else:
        for name in ('pharos_ligand_active', 'pharos_ligand_drug'):
            gene_info[name] = ''
    if 'diseaseCounts' in data:
        gene_info['pharos_disease_all'] = len(data['diseaseCounts'])
        n_alz = 0
        for item in data['diseaseCounts']:
            logger.debug("item: {}".format(item))
            if 'Alzheimer' in item['name']:
                n_alz += 1
        gene_info['pharos_disease_alz'] = n_alz
    else:
        for name in ('pharos_disease_all', 'pharos_disease_alz'):
            gene_info[name] = ''

def extract_value(field, data):
    value = ''
    if ':' in field:
        new_data = copy.deepcopy(data)
        depths = field.split(':')
        logger.info("extract_field: {}".format(depths))
        while len(depths):
            field = depths.pop(0)
            if field in new_data:
                new_data = new_data[field]
            else:
                new_data = dict() # not found so zero out
            logger.info("extract_field: {} {}".format(field, new_data))
        if type(new_data) != type(dict()):
            value = new_data
    else:
        if field in data:
            value = data[field]
    # if list, concatenate
    if type(value) == type(list()):
        if len(value) == 1:
            value = value[0]
        else:
            value = ":".join(value)
    return value

def get_hugo_info(in_file, out_file):

    genes = read_file(in_file)

    if not len(genes):
        logger.critical("Did not read any genes from input file {} ... exiting".format(in_file))
        return

    count = 0
    gene_data = list()

    # do the pipeline for the genes
    for gene in genes:
        count += 1
        gene_info = dict()
        gene_info['query_symbol'] = gene
        logger.debug("\n\nGENE: {}\n".format(gene))
        # HUGO data
        h_data = hugo_api.get_target_info(gene)
        logger.debug("\nHUGO: data_map: {}\n\n results: {}".format(data_map['HUGO'], h_data))
        if len(h_data) > 1:
            logger.critical("HUGO array of results > 1 for {} -- don't know how to handle".format(gene))
            break
        for mapping in data_map['HUGO']:
            gene_info[mapping['target_field']]  = extract_value(mapping['source_field'], h_data[0])
        logger.debug("\n after HUGO gene_info: {}\n".format(gene_info))

        # add to growing list and
        gene_data.append(gene_info)


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
    logger.debug('gene_data: {}'.format(gene_data))
    headers = list()
    for col in csv_map:
        headers.append(col['label'])
    csvwriter.writerow(headers)
    for g_info in gene_data:
        new_row = list()
        for col in csv_map:
            value = ""
            if col['field'] in g_info:
                value = g_info[col['field']]
            new_row.append(value)
        csvwriter.writerow(new_row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='hugo_target_info.py',
        description='find the official target information for a set of genes specified in an file.')
    parser.add_argument('-f', '--file', dest='in_file', help='txt file containing list of gene_symbols', required=True)
    parser.add_argument('-o', '--output', dest='out_file', help='optional output file -- otherwise write to stdout')
    parser.add_argument('-v', '--verbose', dest='verbose', action="store_true", default=False,help='turn on addition information that is sent to stderr')
    parser.add_argument('-d', '--debug', dest='debug', action="store_true", default=False,help='turn on debug information (overides versbose) that is sent to stderr')

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

    get_ade_info(args.in_file, out_file)
