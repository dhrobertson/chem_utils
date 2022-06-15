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
import math

import chem_utils
from chem_utils.data import hugo_api
from chem_utils.data import pharos_api
from chem_utils.data import pharos_api
from chem_utils.data import chembl

ignore = ['gene_symbol']

csv_map = [
    { 'field' : 'query_cpd', 'label': 'QUERY_COMPOUND_ID'},
    { 'field' : 'chembl_results_count', 'label': 'CHEMBL_RESULTS_COUNT'},
    { 'field' : 'chembl_assays_count', 'label': 'CHEMBL_ASSAYS_COUNT'},
    { 'field' : 'chembl_targets_count', 'label': 'CHEMBL_TARGETS_COUNT'},
    { 'field' : 'chembl_active_targets_count', 'label': 'CHEMBL_ACTIVE_TARGETS_COUNT'},
    { 'field' : 'chembl_most_active', 'label' : 'CHEMBL_MOST_ACTIVE'},
    { 'field' : 'chembl_similar_actives_count', 'label' : 'CHEMBL_NON_SELECTIVE_TARGETS_COUNT'},
    { 'field' : 'chembl_similar_actives', 'label' : 'CHEMBL_NON_SELECTIVE_TARGETS'}
]
# TODO -- make script to find this from HUGO
# set this mapping in an easy utility to grab
tgt_map = {
    'CHEMBL2599' : 'SYK',
    'CHEMBL3905' : 'LYN',
    'CHEMBL3234' : 'HCK',
    'CHEMBL5469' : 'PTK2B',
    'CHEMBL4005' : 'PIK3CA',
    'CHEMBL3145' : 'PIK3CB',
    'CHEMBL5554' : 'PIK3C2B',
    'CHEMBL3130' : 'PIK3CD',
    'CHEMBL3267' : 'PIK3CG',
    'CHEMBL2111367' : 'PIK3R1',
    'CHEMBL1075165' : 'PIK3C3',
    'CHEMBL2111432' : 'PIK3CD-2',
    'CHEMBL3880' : 'HSP90AA1',
    'CHEMBL5251' : 'BTK',
    'CHEMBL4247' : 'ALK',
    'CHEMBL2842' : 'mTOR',
    'CHEMBL267' : 'SRC',
    'CHEMBL203' : 'EGFR',
    'CHEMBL1862' : 'ABL1',
    'CHEMBL1824' : 'ERBB2',
    'CHEMBL258' : 'LCK',
    'CHEMBL3038510' : 'PI3KB',
    'CHEMBL1275223' : 'HSC70',
    'CHEMBL2095165' : 'HSP90AA1',
    'CHEMBL1075323' : 'HSP90B1',
    'CHEMBL1981' : 'INSR',
    'CHEMBL2096618' : 'ABL1-II',
    'CHEMBL4852' : 'MAP4K5',
    'CHEMBL3009' : 'ERBB4',
    'CHEMBL1841' : 'FYN',
    'CHEMBL2073' : 'YES',
    'CHEMBL279' : 'FLK1',
    'CHEMBL2695' : 'FAK1',
    'CHEMBL2147' : 'PIM1',
    'CHEMBL5719' : 'CDK8',
    'CHEMBL2276' : 'JNK1',
    'CHEMBL1957' : 'IGF1R',
    'CHEMBL5568' : 'ROS1',
    'CHEMBL1974' : 'FLT3'
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
    compounds = list()
    extra_data = {
        'fields' : list(),
        'data'   : dict()
    }
    if not os.path.isfile(in_file):
        logger.critical("Unable to locate defined input file \"" + in_file +"\" for compounds ... skipping")
        return (compounds, extra_data)

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
                cpd_id = row[0]
                cpd_id = cpd_id.replace(" ", "")
                if cpd_id not in extra_data['data']:
                    extra_data['data'][cpd_id] = list()
                if cpd_id not in compounds:
                    compounds.append(cpd_id)
                data = dict()
                for i in range(0, len(headers)):
                    data[headers[i]] = row[i]
                extra_data['data'][cpd_id].append(data)
            n += 1
        logger.info("Read " + str(n) + " lines from file \"" + str(in_file) + "\"" )
    fp.close()

    logger.debug("compounds: {} extra_data: {}".format(compounds, extra_data))

    return (compounds, extra_data)

# clean up the activity records
def clean_activity(activity):
    tgt = activity['target_chembl_id']
    assay = activity['assay_chembl_id']
    active = 'UNKNOWN'
    pIC50 = None
    original = "\"{}\" {} {}".format(activity['type'], activity['value'], activity['units'])
    #logger.info("type: {} regular: {} {} standard: {} {}".format(
    #    activity['type'],
    #    activity['value'], activity['units'],
    #    activity['standard_value'], activity['standard_units']
    #))
    try:
        if activity['standard_units'] and activity['standard_value']:
            if activity['type'] in ['pIC50']:
                pIC50 = activity['standard_value']
            elif activity['type'] in ['IC50', 'EC50'] and activity['standard_units'] == 'nM':
                pIC50 = float(int(100*(9.0 - math.log10(float(activity['standard_value']))))/100.0)
            elif activity['type'] not in ['INH', 'Inhibition', 'Potency', 'Activity', 'Ratio IC50', 'Ratio',
                                          'pKi', 'K', 'Kd', 'KD', 'kon', 'koff', 'Ki', 'KD', 'Kd apparent', 'Ki(app)',
                                          'GI50', 'Kinact', 'DC50', 't1/2',
                                          'deltaTm', 'Time', 'Thermal melting change', 'Residual Activity']:
                logger.error("Unable to process type \"{}\"".format(original))
        else:
            logger.debug("Issue with value or units -- skipping \"{}\"".format(original))
    except:
        logger.debug("Unable to process bioactivity \"{}\"".format(original))
    #set active flag < 10uM -> pIC50 5
    if pIC50:
        try:
            if pIC50 < 5.5:
                active = 'INACTIVE'
            else:
                active = 'ACTIVE'
        except:
            logger.debug('Issue with pIC50 {} type {}'.format(pIC50, type(pIC50)))
    #print("tgt: {} assay: {} active: {} pIC50: {} original: \"{}\"".format(tgt, assay, active, pIC50, original))

    return(tgt, assay, active, pIC50, original)

# get the information through API for the selected compound
def get_cpd_info(cpd):
    cpd_info = dict()
    logger.debug("\n\nCPD: {}\n".format(cpd))

    #CHEMBL data
    logger.level = logging.ERROR
    activities = chembl.get_bioactivities_for_molecules(cpd)
    logger.level = logging.INFO
    cpd_info['chembl_results_count'] = len(activities)
    actives = dict()
    counts = {
        'targets': dict(),
        'assays': dict()
    }
    if len(activities):
        logger.info("CPD: {} found {} results".format(cpd, len(activities)))
        for activity in activities:
            if activity['assay_type'] == 'B' and activity['target_organism'] == 'Homo sapiens':
                (target, assay, active, pIC50, original) = clean_activity(activity)
                if target not in counts['targets']:
                    counts['targets'][target] = 0
                if assay not in counts['assays']:
                    counts['assays'][assay] = 0
                counts['targets'][target] += 1
                counts['assays'][assay] += 1
                if active == 'ACTIVE':
                    gene_symbol = None
                    if target in tgt_map:
                        gene_symbol = tgt_map[target]
                    else:
                        gene_symbol = target
                        #print("unknown mapping target: {} target_map: {}".format(target, gene_symbol))
                    if gene_symbol not in actives:
                        actives[gene_symbol] = list()
                    actives[gene_symbol].append(pIC50)
            else:
                logger.debug('ignoring \"{}\" \"{}\"'.format(activity['assay_type'],activity['type']))
    else:
        logger.critical("ChEMBL unable to find activities for {} ... skipping ".format(cpd))

    # work through the actives to fill out the final columns
    max_pIC50 = dict()
    targets = list(actives.keys())
    for target in targets:
        max = 1.0
        for pIC50 in actives[target]:
            if pIC50 > max:
                max = pIC50
        max_pIC50[max] = target
    #print(max_pIC50)
    pIC50s = list(max_pIC50.keys())
    if len(pIC50s):
        #print(pIC50s)
        pIC50s_sorted = sorted(pIC50s, reverse=True)
        #print(pIC50s_sorted)
        cpd_info['chembl_most_active'] = "{} (pIC50: {})".format(max_pIC50[pIC50s_sorted[0]], pIC50s_sorted[0])
        non_selective = list()
        for pIC50 in pIC50s_sorted:
            if pIC50 != pIC50s_sorted[0] and pIC50 - pIC50s_sorted[0] < 1.01:
                non_selective.append("{} (pIC50: {})".format(max_pIC50[pIC50], pIC50))
        cpd_info['chembl_similar_actives'] = ' : '.join(non_selective)
        cpd_info['chembl_similar_actives_count'] = len(non_selective)
    cpd_info['chembl_targets_count'] = len(counts['targets'])
    cpd_info['chembl_active_targets_count'] = len(actives)
    cpd_info['chembl_assays_count'] = len(counts['assays'])
    logger.level = logging.DEBUG

    logger.debug("CPD: {} cpd_info: {}".format(cpd, cpd_info))

    return cpd_info

# read a json mapping file
def get_ade_info(in_file, out_file):

    (compounds, extra_data) = read_file(in_file)

    if not len(compounds):
        logger.critical("Did not read any compounds from input file {} ... exiting".format(in_file))
        return

    # cache results
    cpd_cache = get_cache()
    cpd_data = dict()

    count = 0
    # do the pipeline for the genes
    for cpd in compounds:
        count += 1
        cpd_info = dict()
        if cpd in cpd_cache:
            cpd_info = cpd_cache[cpd]
        else:
            cpd_info = get_cpd_info(cpd)
            cpd_info['query_cpd'] = cpd
            #at the risk of more I/O -- save cache
            cpd_cache[cpd] = cpd_info
            #save_cache(cpd_cache)
        # add to growing list and
        cpd_data[cpd] = cpd_info

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
    for field in extra_data['fields']:
        headers.append(field)
    for col in csv_map:
        headers.append(col['label'])
    csvwriter.writerow(headers)
    for cpd in compounds:
        for row in extra_data['data'][cpd]:
            new_row=list()
            for field in extra_data['fields']:
                value = ""
                if field in row:
                    value = row[field]
                new_row.append(value)
            for col in csv_map:
                value = ""
                if col['field'] in cpd_data[cpd]:
                    value = cpd_data[cpd][col['field']]
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

    get_ade_info(args.in_file, out_file)
