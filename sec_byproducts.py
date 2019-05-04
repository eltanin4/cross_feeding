import numpy as np
from load_data import *

def secByproducts(orgRxns, rxnMat, prodMat, Core):
    """
    Given a set of reactions in an organism, calculates which
    metabolites (products of said reactions) are not reactants
    in any other reactions in the set. Assumes that these metabolites
    will be secreted.

    RETURNS: 

    A vector of secreted metabolites, given by custom IDs.
    """

    # Calculating first the vectors of reactants and products.
    reactVec = np.logical_or.reduce( rxnMat[ orgRxns ] ) * 1
    prodVec = np.logical_or.reduce( prodMat[ orgRxns ] ) * 1

    # The byproducts are those that are products but not reactants.
    bypVec =  np.logical_and( np.logical_xor( reactVec, prodVec ), 
                              prodVec )
    bypVec[ Core + Currency + Energy ] = 0

    return bypVec

def pathwaySecByproducts(pathRxns, orgRxns, rxnMat, prodMat, Core):
    """
    Given a set of reactions in an organism, calculates which
    metabolites (products of said reactions) are not reactants
    in any other reactions in the set. Assumes that these metabolites
    will be secreted.

    RETURNS: 

    A vector of secreted metabolites, given by custom IDs.
    """

    # Calculating first the vectors of reactants and products.
    reactVec = np.logical_or.reduce( rxnMat[ orgRxns ] ) * 1
    prodVec = np.logical_or.reduce( prodMat[ pathRxns ] ) * 1

    # The byproducts are those that are products but not reactants.
    bypVec =  np.logical_and( np.logical_xor( reactVec, prodVec ), 
                              prodVec )
    bypVec[ Core + Currency + Energy ] = 0

    return bypVec
