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
    { 'field' : 'hgnc_refseq_accession', 'label': 'HGNC_REFSEQ_ACCESSION'},
    { 'field' : 'pharos_fam', 'label': 'PHAROS_GENE_FAMILY'},
    { 'field' : 'pharos_tdl', 'label': 'PHAROS_TDL'},
    { 'field' : 'pharos_novelty', 'label': 'PHAROS_NOVELTY'},
    { 'field' : 'pharos_ligand_drug', 'label': 'PHAROS_COUNT_DRUGS'},
    { 'field' : 'pharos_ligand_active', 'label': 'PHAROS_COUNT_ACTIVES'},
    { 'field' : 'pharos_disease_all', 'label': 'PHAROS_COUNT_DISEASES'},
    { 'field' : 'pharos_disease_alz', 'label': 'PHAROS_COUNT_ALZHEIMERS'},
    { 'field' : 'ot_ensembl_id', 'label': 'OT_ENSEMBL_ID'},
    { 'field' : 'ot_tract_sm_hqc', 'label': 'OT_TRACT_SM_HIGHQUALITYCOMPOUNDS'},
    { 'field' : 'ot_tract_sm_tc', 'label': 'OT_TRACT_SM_TOPCATEGORY'},
    { 'field' : 'ot_tract_ab_tc', 'label': 'OT_TRACT_AB_TOPCATEGORY'},
    { 'field' : 'chembl_target_id', 'label': 'CHEMBL_TARGET_ID'},
    { 'field' : 'chembl_assay_count', 'label': 'CHEMBL_ASSAY_COUNT'},
    { 'field' : 'chembl_results_count', 'label': 'CHEMBL_RESULTS_COUNT'},
    { 'field' : 'chembl_molecule_count', 'label': 'CHEMBL_MOLECULE_COUNT'},
    { 'field' : 'chembl_pdb_count', 'label': 'CHEMBL_PDB_COUNT'}
]

data_map = {
    'HUGO' : [
        { 'source_field': 'symbol', 'target_field': 'hgnc_symbol' },
        { 'source_field': 'hgnc_id', 'target_field': 'hgnc_id' },
        { 'source_field': 'name', 'target_field': 'hgnc_name' },
        { 'source_field': 'entrez_id', 'target_field': 'hgnc_entrez_id' },
        { 'source_field': 'ensembl_gene_id', 'target_field': 'hgnc_ensembl_gene_id' },
        { 'source_field': 'uniprot_ids', 'target_field': 'hgnc_uniprot_ids' },
        { 'source_field': 'refseq_accession', 'target_field': 'hgnc_refseq_accession'}
    ],
    'PHAROS' : [
        { 'source_field': 'tdl', 'target_field': 'pharos_tdl' },
        { 'source_field': 'fam', 'target_field': 'pharos_fam' },
        { 'source_field': 'novelty', 'target_field': 'pharos_novelty' }
    ],
    'OPENTARGET' : [
        { 'source_field': 'id', 'target_field': 'ot_ensembl_id' },
        { 'source_field': 'tractability:smallmolecule:high_quality_compounds', 'target_field': 'ot_tract_sm_hqc' },
        { 'source_field': 'tractability:smallmolecule:top_category', 'target_field': 'ot_tract_sm_tc' },
        { 'source_field': 'tractability:antibody:top_category', 'target_field': 'ot_tract_ab_tc' }
    ],
    'CHEMBL' : [
        { 'source_field': 'chembl_pdb_count', 'target_field': 'chembl_pdb_count' },
        { 'source_field': 'chembl_pdb_score', 'target_field': 'chembl_pdb_score' }
    ]
}

# cache results
cache_file = 'gene_cache.json'
def get_cache():
    gene_cache = dict()
    if not os.path.isfile(cache_file):
        logger.info("Unable to locate cache file \"" + cache_file +"\" for genes ... skipping")
        return gene_cache

    with open(cache_file) as fp:
        gene_cache = json.load(fp)
    fp.close()

    return gene_cache

def save_cache(gene_cache):
    with open(cache_file, "w") as fp:
        json.dump(gene_cache, fp, indent=2)
    fp.close()

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

# add assay data from ChEMBL

def chembl_add_assay_info(target_chembl_id, gene_info):
    mol_field = 'molecule_chembl_id'
    asy_field = 'assay_chembl_id'
    if False: # api
        records = chembl_api.get_bioactivities_for_targets(target_chembl_id)
    else: # db
        mol_field = 'molregno'
        asy_field = 'assay_id'
        tids_dict = db_chembl.tids_from_chembl_ids(target_chembl_id)
        tids = list(tids_dict.values())
        if len(tids) == 1:
            tids = tids[0]
        assays = db_chembl.assay_ids_from_tids(tids)
        #print(tids_dict, tids, assays)
        #exit()
        records = list()
        if assays:
            records = db_chembl.activities_from_assay_ids(assays)
        #print(tids, assays)
        #exit()
    gene_info['chembl_results_count'] = len(records)
    counts_asy = dict()
    counts_mol = dict()
    for record in records:
        #logger.debug("    ChEMBL record {}".format(record))
        mol_id = record[mol_field]
        asy_id = record[asy_field]
        if mol_id not in counts_mol:
            counts_mol[mol_id] = 0
        counts_mol[mol_id] += 1
        if asy_id not in counts_asy:
            counts_asy[asy_id] = 0
        counts_asy[asy_id] += 1
    gene_info['chembl_assay_count'] = len(counts_asy)
    gene_info['chembl_molecule_count'] = len(counts_mol)
    gene_info['chembl_ligand_score'] = 0.25*float(min(len(records), 100))/100.0 + 0.25*float(min(len(counts_asy), 3))/3.0 + 0.50*float(min(len(counts_mol), 20))/20.0

# extract the data from chembl
def process_chembl_data(g_symbol, c_data):
    logger.debug("\n*** ChEMBL data: {}\n".format(len(c_data)))

    extracted_data = dict()
    # PDB from xref
    extracted_data['chembl_pdb_count'] = 0
    extracted_data['chembl_pdb_score'] = ""
    if "target_components" in c_data[0]:
        subset = c_data[0]["target_components"]
        logger.debug(" targets_components_count: {}".format(len(subset)))
        if len(subset) > 1:
            logger.critical("CHEMBL ERROR: GENE_SYMBOL {} Do not know how to handle more than one -- so taking first -- check".format(g_symbol))
        if "target_component_xrefs" in subset[0]:
            xrefs = subset[0]["target_component_xrefs"]
            pdb_count = 0
            for xref in xrefs:
                logger.debug("    xref: {}".format(xref))
                if xref['xref_src_db'] == 'PDBe':
                    pdb_count += 1
            extracted_data['chembl_pdb_count'] = pdb_count
            if pdb_count > 0:
                extracted_data['chembl_pdb_score'] = "LOW"
            if pdb_count > 10:
                extracted_data['chembl_pdb_score'] = "MED"
            if pdb_count > 30:
                extracted_data['chembl_pdb_score'] = "HIGH"

    return extracted_data

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

#@smart_cache
# read a json mapping file
def get_ade_info(in_file, out_file):

    genes = read_file(in_file)

    if not len(genes):
        logger.critical("Did not read any genes from input file {} ... exiting".format(in_file))
        return

    # cache results
    gene_cache = get_cache()
    gene_data = list()

    count = 0
    # do the pipeline for the genes
    for gene in genes:
        count += 1
        if gene in gene_cache:
            gene_info = gene_cache[gene]
        else:
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

            #PHAROS data - use HUGO gene symbol
            if True and gene_info['hgnc_symbol']:
                p_data = pharos_api.get_target_info(gene_info['hgnc_symbol'])
                if len(p_data) > 1:
                    logger.critical("PHAROS array of results > 1 for {} -- don't know how to handle".format(gene))
                    break
                if len(p_data) == 1:
                    logger.debug("\nPHAROS: data_map: {} results: {}\n".format(data_map['PHAROS'], p_data))
                    for mapping in data_map['PHAROS']:
                        gene_info[mapping['target_field']] = extract_value(mapping['source_field'],p_data[0])
                    pharos_process_counts(gene_info, p_data[0])
                else:
                    logger.critical("PHAROS array of results == 0 for {} -- skipping".format(gene))

                logger.debug("\n after PHAROS gene_info: {}\n".format(gene_info))

            #OPENTARGET data - use HUGO gene symbol
            if True and gene_info['hgnc_symbol']:
                ot_data = opentarget_json.get_target_info(gene_info['hgnc_symbol'])
                logger.debug("\nOPENTARGET: ot_data: {} \n".format(ot_data))
                if len(ot_data) > 1:
                    logger.critical("OPENTARGET array of results > 1 for {} -- don't know how to handle".format(gene))
                    break
                if len(ot_data) == 1:
                    logger.debug("\nOPENTARGET: data_map: {} results: {}\n".format(data_map['OPENTARGET'], ot_data))
                    for mapping in data_map['OPENTARGET']:
                        gene_info[mapping['target_field']] = extract_value(mapping['source_field'], ot_data[0])
                else:
                    logger.critical("OPENTARGET array of results == 0 for {} -- skipping".format(gene))

                logger.debug("\n after OPENTARGET gene_info: {}\n".format(gene_info))

            #CHEMBL data - use HUGO gene symbol
            if True and gene_info['hgnc_symbol'] and gene_info['hgnc_uniprot_ids']:
                g_symbol = gene_info['hgnc_symbol']
                uniprot_id = gene_info['hgnc_uniprot_ids']
                if ':' in uniprot_id:   # multiple uniprot_ids
                    logger.critical("\nChEMBL multiple uniprot_ids for gene {} uniprot_ids {} .. skipping".format(g_symbol, uniprot_id))
                else:
                    targets = chembl_api.find_target_ids_by_accession(uniprot_id)
                    if len(targets):
                        target_chembl_id = get_single_protein(targets)
                        logger.debug("\nChEMBL {} => {}".format(gene, target_chembl_id))
                        if target_chembl_id:
                            # extract PDB from target_details
                            gene_info['chembl_target_id'] = target_chembl_id
                            c_data = chembl_api.get_target_details(target_chembl_id)
                            c_processed_data = process_chembl_data(g_symbol, c_data)
                            for mapping in data_map['CHEMBL']:
                                gene_info[mapping['target_field']] = extract_value(mapping['source_field'], c_processed_data)
                            # get assay info
                            chembl_add_assay_info(target_chembl_id, gene_info)
                        else:
                            logger.critical("ChEMBL unable to find chembl_id for {}".format(g_symbol))
                    else:
                        logger.critical("ChEMBL unable to find access \"hgnc_uniprot_ids for {} ... skipping ".format(g_symbol))
            logger.debug("\n after CHEMBL gene_info: {}\n".format(gene_info))
            if len(genes)>10:
                logger.warning("\r Gene {:3d} of {}".format(count, len(genes)))

            #at the risk of more I/O -- save cache
            gene_cache[gene] = gene_info
            save_cache(gene_cache)
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
    parser = argparse.ArgumentParser(prog='get_ade_data.py',
        description='compile the ADE target information for a set of genes specified in an file.')
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
