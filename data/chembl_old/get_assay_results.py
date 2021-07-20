"""
get_assay_results.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utility functions for returning assay results
"""
import pprint
from .. import chembl

def by_molecules(molecule_chembl_ids):
    """ retrieve the assay results by molecule_chembl_ids """

    return chembl.get_bioactivities_for_molecules(molecule_chembl_ids)

def by_targets(target_chembl_ids):
    """ retrieve the assay results by target_chembl_ids """

    return chembl.get_bioactivities_for_targets(target_chembl_ids)

def get_targets(chembl_assay_results):
    """ retrieve the unique targets in chembl_assay_results """
    targets = dict()
    for a in chembl_assay_results:
        if 'target_chembl_id' in a:
            tgt = a['target_chembl_id']
            cpd = a['molecule_chembl_id']
            asy = a['assay_chembl_id']
            if tgt not in targets:
                targets[tgt] = {
                    'id' : tgt,
                    'id_type' : 'ChEMBL',
                    'name' : a['target_pref_name'],
                    'symbol' : None,
                    'organism' : a['target_organism'],
                    'molecules' : list(),
                    'assays' : list(),
                    "number_results" : 0
                }
            if cpd not in targets[tgt]['molecules']:
                targets[tgt]['molecules'].append(cpd)
            if asy not in targets[tgt]['assays']:
                targets[tgt]['assays'].append(asy)
            targets[tgt]['number_results'] += 1

    # add in the unique number for assays, molecules, results
    # TODO: use a utility function with cache so do not need to always recompute
    all_targets = list(targets.keys())
    target_details = chembl.get_target_details(all_targets)
    target_info = dict()
    for tgt in target_details:
        target_info[tgt['target_chembl_id']] = tgt

    targets_list = list()
    for tgt in targets.values():
        target_chembl_id = tgt['id']
        if target_chembl_id in target_info:
            #print(pprint.pformat(target_info[target_chembl_id]))
            tgt['gene_symbol'] = target_info[target_chembl_id]['gene_symbol']
        tgt['number_molecules'] = len(tgt['molecules'])
        tgt['number_assays'] = len(tgt['assays'])
        targets_list.append(tgt)

    return targets_list

def get_molecules(chembl_assay_results):
    """ retrieve the unique molecules in chembl_assay_results """
    molecules = dict()
    for a in chembl_assay_results:
        if 'molecule_chembl_id' in a:
            tgt = a['target_chembl_id']
            cpd = a['molecule_chembl_id']
            asy = a['assay_chembl_id']
            if cpd not in molecules:
                molecules[cpd] = {
                    'id' : cpd,
                    'id_type' : 'ChEMBL',
                    'name' : a['molecule_pref_name'],
                    'smiles' : a['canonical_smiles'],
                    'targets' : list(),
                    'assays' : list(),
                    'number_results' : 0
                }
            if tgt not in molecules[cpd]['targets']:
                molecules[cpd]['targets'].append(tgt)
            if asy not in molecules[cpd]['assays']:
                molecules[cpd]['assays'].append(asy)
            molecules[cpd]['number_results'] += 1

    # add in the unique number for assays, molecules, results
    # TODO: use a utility function with cache so do not need to always recompute
    molecules_list = list()
    for mol in molecules.values():
        mol['number_targets'] = len(mol['targets'])
        mol['number_assays'] = len(mol['assays'])
        molecules_list.append(mol)

    return molecules_list

def get_assays(chembl_assay_results):
    """ retrieve the unique assays in chembl_assay_results """
    assays = dict()
    for a in chembl_assay_results:
        if 'assay_chembl_id' in a:
            tgt = a['target_chembl_id']
            cpd = a['molecule_chembl_id']
            asy = a['assay_chembl_id']
            if asy not in assays:
                assays[asy] = {
                    'id' : asy,
                    'id_type' : 'ChEMBL',
                    'target' : tgt,
                    'type' : a['assay_type'],
                    'desc' : a['assay_description'],
                    'std_units' : a['standard_units'],
                    'std_type'  : a['standard_type'],
                    'molecules' : list(),
                    'number_results' : 0
                }
            if cpd not in assays[asy]['molecules']:
                assays[asy]['molecules'].append(cpd)
            assays[asy]['number_results'] += 1

    # add in the unique number for assays, molecules, results
    # TODO: use a utility function with cache so do not need to always recompute
    assays_list = list()
    for asy in assays.values():
        asy['number_molecules'] = len(asy['molecules'])
        assays_list.append(asy)

    return assays_list

def get_results(chembl_assay_results):
    """ retrieve the results in chembl_assay_results """
    results = list()
    results_details = list()
    for a in chembl_assay_results:
        if 'activity_id' in a:
            r_id = a['activity_id']
            if r_id not in results:
                results.append(r_id)
                result = {
                    'id' : r_id,
                    'id_type' : 'ChEMBL',
                    'target_id' : a['target_chembl_id'],
                    'target_id_type' : 'ChEMBL',
                    'molecule_id' : a['molecule_chembl_id'],
                    'molecule_id_type' : 'ChEMBL',
                    'assay_id' : a['assay_chembl_id'],
                    'assay_id_type' : 'ChEMBL',
                    'document_id' : a['document_chembl_id'],
                    'document_id_type' : 'ChEMBL',
                    'type' : a['type'],
                    'value' : a['value'],
                    'units' : a['units'],
                    'reference' : '{} {}'.format(a['document_journal'], a['document_year'])
                }
                results_details.append(result)
    return results_details

def convert_to_standard(chembl_assay_results):
    researchDataSet = dict()
    researchDataSet['targets'] = get_targets(chembl_assay_results)
    researchDataSet['molecules'] = get_molecules(chembl_assay_results)
    researchDataSet['assays'] = get_assays(chembl_assay_results)
    researchDataSet['results'] = get_results(chembl_assay_results)
    return researchDataSet
