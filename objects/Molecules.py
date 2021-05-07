"""
Molecules.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

Class for handling sets of molecules -- wrapper class for RDKit
"""
import pprint
import logging
from chem_utils.utils import core_utils
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
logger = logging.getLogger()

class Molecules:
    """Molecules List Molecules"""

    def __init__(self):
        self.molecule_list = list()
        self.count = 0

    def length(self):
        return len(self.molecule_list)

    # adds as either a single value or list
    # it can be a tuple of (smi, name) or just
    def addSmiles(self, smiles):
        rdDepictor.SetPreferCoordGen(True)
        if type(smiles) != type(list()):
            smiles = [smiles]
        errors = list()
        for smi in smiles:
            # check if tuple
            name = None
            smiles = smi
            if type(smi) == type(tuple()):
                (smiles, name) = (smi[0], smi[1])
            if name == None:
                name = "mol_name_{:05d}".format(self.count)
                self.count += 1
            logger.debug('Molecules.addSmiles: smiles: {} name: {} '.format(smiles, name))
            # convert to Canonical Smiles - this also checks for valid structure
            try:
                c_smiles = Chem.CanonSmiles(smiles)
                cs_mol   = Chem.MolFromSmiles(c_smiles)
                fp       = FingerprintMols.FingerprintMol(cs_mol)
                d2d      = rdMolDraw2D.MolDraw2DSVG(600,600)
                d2d.DrawMolecule(cs_mol)
                d2d.FinishDrawing()
                mol = {
                    "name" : name,
                    "smiles" : c_smiles,
                    "mol" : cs_mol,
                    "fp"  : fp,
                    "svg" : d2d.GetDrawingText()
                }
                self.molecule_list.append(mol)
            except Exception as e:
                if not name:
                    name = "Undefined"
                logger.error("Issue adding ({}, {}) for reason \"{}\"..skipping".format(smiles, name, e))
                errors.append((smiles, name))

        return errors

    # return as tuple list
    def getAllSmilesList(self):
        all = list()
        for mol in self.molecule_list:
            smiles = (mol["smiles"], mol["name"])
            all.append(smiles)
        return all

    # read from smiles file
    def addFromSmiFile(self, in_file):
        if not core_utils.file_exists(in_file):
            logger.error("Unable to locate file \"{}\" .. skipping".format(in_file))
            return None
        initial = self.length()
        with open(in_file) as file_ref:
            for line in file_ref.readlines():
                parts  = line.rstrip('\n').split(' ')
                smiles = parts.pop(0)
                name = None
                if len(parts):
                    name = parts.pop(0)
                self.addSmiles((smiles, name))

        number_added = self.length() - initial
        logger.info("Added {} molecules from file \"{}\"".format(number_added, in_file))
        file_ref.close()

    # get moleule by index -- surrogate for iterator TODO: Generate Iterator function
    def getMol(self, index):
        try:
            mol = self.molecule_list[index]
            return mol
        except:
            logger.error("Molecules.getMol: Unable to get molecule index \"{}\" - valid indices are -{} to {}".format(index, self.length()-1, self.length()-1))
            return None

    def getIntraSimilarityMatrix(self):
        # set up lists and dictionaries
        logger.debug("Molecules.getSimilarityMatrix: len: {}".format(self.length()))
        name_list = list()
        similarity_matrix = dict()
        for mol in self.molecule_list:
            similarity_matrix[mol["name"]] = dict()
            name_list.append(mol["name"])
        # compute similarities -- lower triangle and reflect
        for i in range(0, len(self.molecule_list)):
            mol_i = self.molecule_list[i]
            similarity_matrix[mol_i["name"]][mol_i["name"]] = 1.0
            for j in range(i+1, len(self.molecule_list)):
                mol_j = self.molecule_list[j]
                #sim = compare_molecules.compute_similarity(mol_i["smiles"], mol_j["smiles"])
                sim = DataStructs.FingerprintSimilarity(mol_i["fp"], mol_j["fp"])
                sim = float(int(sim*1000)/1000) # truncate to 3 decimals
                similarity_matrix[mol_i["name"]][mol_j["name"]] = sim
                similarity_matrix[mol_j["name"]][mol_i["name"]] = sim
                logger.debug("Molecules.getIntraSimilarities: i j name_i name_j sim: {} {} {} {} {}".format(i,j,mol_i["name"],mol_j["name"],sim))

        logger.debug("Molecules.getIntraSimilarities: matrix: {}".format(pprint.pformat(similarity_matrix)))
        return similarity_matrix


    # compute intra molecule list similarities
    def getIntraSimilarities(self):
        # set up lists and dictionaries
        logger.debug("Molecules.getIntraSimilarities: len: {}".format(self.length()))
        name_list = list()
        similarity_matrix = self.getIntraSimilarityMatrix()
        for mol in self.molecule_list:
            name_list.append(mol["name"])

        # generate into a list of mol_i mol_j sim
        similarities = list()
        for name_1 in name_list:
            max_name = None
            max_sim = None
            for name_2 in similarity_matrix[name_1].keys():
                if name_1 != name_2:
                    if max_sim == None or similarity_matrix[name_1][name_2] > max_sim:
                        max_name = name_2
                        max_sim = similarity_matrix[name_1][name_2]
                logger.debug("Molecules.getIntraSimilarities: name_1 name_2 sim max_name max_sim: {} {} {} {} {}".format(name_1, name_2, similarity_matrix[name_1][name_2], max_name, max_sim))
            logger.debug("Molecules.getIntraSimilarities: name_1 max_name max_sim: {} {} {} ".format(name_1, max_name, max_sim))
            similarities.append("{} {} {}".format(name_1, max_name, max_sim))

        return similarities

    # compute inter molecule list similarities
    def getInterSimilarities(self, mol_ref):
        # set up lists and dictionaries
        logger.debug("Molecules.getInterSimilarities: len: {} compared_to: {}".format(self.length(), mol_ref.length()))
        name_list = list()
        similarity_matrix = dict()
        for mol in self.molecule_list:
            similarity_matrix[mol["name"]] = dict()
            name_list.append(mol["name"])
        # compute similarities -- lower triangle and reflect
        for i in range(0, len(self.molecule_list)):
            mol_i = self.molecule_list[i]
            for j in range(0, mol_ref.length()):
                mol_j = mol_ref.getMol(j)
#                sim = compare_molecules.compute_similarity(mol_i["smiles"], mol_j["smiles"])
                sim = DataStructs.FingerprintSimilarity(mol_i["fp"], mol_j["fp"])
                sim = float(int(sim*1000)/1000) # truncate to 3 decimals
                similarity_matrix[mol_i["name"]][mol_j["name"]] = sim
                logger.debug("Molecules.getInterSimilarities: i j name_i name_j sim: {} {} {} {} {}".format(i,j,mol_i["name"],mol_j["name"],sim))

        logger.debug("Molecules.getInterSimilarities: matrix: {}".format(pprint.pformat(similarity_matrix)))
        # generate into a list of mol_i mol_j sim
        similarities = list()
        for name_1 in name_list:
            max_name = None
            max_sim = None
            for name_2 in similarity_matrix[name_1].keys():
                if max_sim == None or similarity_matrix[name_1][name_2] > max_sim:
                    max_name = name_2
                    max_sim = similarity_matrix[name_1][name_2]
                logger.debug("Molecules.getIntraSimilarities: name_1 name_2 sim max_name max_sim: {} {} {} {} {}".format(name_1, name_2, similarity_matrix[name_1][name_2], max_name, max_sim))
            logger.debug("Molecules.getInterSimilarities: name_1 max_name max_sim: {} {} {} ".format(name_1, max_name, max_sim))
            similarities.append("{} {} {}".format(name_1, max_name, max_sim))

        return similarities
