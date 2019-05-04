import numpy as np
import re
import pickle

stoich_matrix = np.genfromtxt(
    '../new_kegg/data_files/pure_stoich_matrix.txt')
display_lookup = pickle.load(
    open('../new_kegg/data_files/rxn_string_dict.dat', 'rb'))
rxn_map = rxn_kegg_to_id = pickle.load(
    (open('../new_kegg/data_files/pure_kegg_to_self_rxn_map.dat', 'rb')))
inv_rxn_map = rxn_id_to_kegg = {value: key for key, value in rxn_map.items()}

def drxn(rxn):
    return new_disp_reaction(rxn, stoich_matrix, inv_rxn_map, display_lookup)


def new_disp_reaction(reaction, stoich_matrix, inv_rxn_map, display_lookup):
    reverse = False
    rxn_ind = np.where(np.all(reaction == stoich_matrix, axis=1))[0][0]
    if rxn_ind >= (np.shape(stoich_matrix)[0] / float(2)):
        rxn_ind -= (np.shape(stoich_matrix)[0] / float(2))
        reverse = True
    if rxn_ind >= (np.shape(stoich_matrix)[0] / float(2)):
        rxn_ind = 0
        reverse = True
    kegg_rxn_id = inv_rxn_map[rxn_ind]
    if not reverse:
        reactant_string, product_string = re.split(
            r'\s*&lt;=&gt;\s*', display_lookup[kegg_rxn_id])
    else:
        product_string, reactant_string = re.split(
            r'\s*&lt;=&gt;\s*', display_lookup[kegg_rxn_id])
    return reactant_string + '  ===>  ' + product_string

def give_id(reaction):
    reverse = False
    rxn_ind = np.where(np.all(reaction == stoich_matrix, axis=1))[0][0]
    if rxn_ind >= (np.shape(stoich_matrix)[0] / float(2)):
        rxn_ind -= (np.shape(stoich_matrix)[0] / float(2))
        reverse = True
    if rxn_ind >= (np.shape(stoich_matrix)[0] / float(2)):
        rxn_ind = 0
        reverse = True
    kegg_rxn_id = inv_rxn_map[rxn_ind]
    return kegg_rxn_id
