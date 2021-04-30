"""
compounds.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Specific utilities for working with chembl data
"""
import logging
from rdkit import Chem
from data.chembl import utils

logger = logging.getLogger()

def get_by_targets(target_chembl_ids):
    """ retrieve the molecule_compound_ids for targets """
    record_list = list()
    molecules = dict()
    records = utils.get_bioactivities_for_targets(target_chembl_ids)
    for record in records:
        mol_id = record['molecule_chembl_id']
        molecules[mol_id] = '1'

    molecules_list = list(molecules.keys())

    return molecules_list

def get_details(molecule_chembl_ids):
    """ get the molecule details for the molecule_chembl_ids """

    # single as string -- return as element not list
    if type(molecule_chembl_ids) == type(str()):
        return utils.get_molecule_details(molecule_chembl_ids)

    # small list -- do as single chunk
    if len(molecule_chembl_ids) < 50:
        return utils.get_molecule_details(molecule_chembl_ids)

    # perform in chunks
    molecules = list()
    chunk_length = 50
    n_ids = len(molecule_chembl_ids)
    while len(molecule_chembl_ids):
        logger.info('{} total: {} left'.format(n_ids, len(molecule_chembl_ids)))
        subset = list()
        while len(molecule_chembl_ids) and len(subset) < chunk_length:
            subset.append(molecule_chembl_ids.pop(0))
        subset_molecules = utils.get_molecule_details(subset)
        for molecule in subset_molecules:
            molecules.append(molecule)
    return molecules

def get_sdf(molecule_chembl_ids):
    """ retrieve the details for the molecule_compound_ids and create SDF """
    sdf = ""
    molecules = get_details(molecule_chembl_ids)
    for molecule in molecules:
        m = Chem.MolFromSmiles(molecule['molecule_structures']['canonical_smiles'])
        m.SetProp("_Name", molecule['molecule_chembl_id'])
        # set cross references
        for item in molecule['cross_references']:
            if item['xref_id'] and item['xref_id']:
                m.SetProp(item['xref_src'], item['xref_id'])
        sdf = sdf + Chem.MolToMolBlock(m) + '$$$$\n'

    return sdf

def get_smi(molecule_chembl_ids):
    """ retrieve the details for the molecule_compound_ids and create SDF """
    smi = ""
    molecules = get_details(molecule_chembl_ids)
    for molecule in molecules:
        smi = smi + molecule['molecule_structures']['canonical_smiles'] + " " + molecule['molecule_chembl_id'] + "\n"

    return smi
