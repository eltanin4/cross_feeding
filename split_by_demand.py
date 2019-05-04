import numpy as np
from copy import deepcopy

def isLimiting(tRct, tRxn, m, S, reactants):
    for oRct in reactants:
        if (m[tRct] * S[tRxn, oRct] / S[tRxn,tRct]) > m[oRct]:
            return True

def giveLimitingCurrency(r, tRxn):
    return np.where(r[tRxn] == max(r[tRxn][np.where(r[tRxn] < 0.0)]))[0][0]

def splitByDemand(stoich_matrix, rxnMat, prodMat, sumRxnVec, rho, pi, 
                 nutrientSet, Energy, Currency, Core, orgRxns):
    # Initializing yield counters for E and B.
    runningE, runningB = 0.0, 0.0

    # Getting a vector of which reactions are in the scope of the network.
    scopeRxns = np.zeros(len(stoich_matrix))
    scopeRxns[orgRxns] = 1
    isCoreProduced = np.zeros(len(np.transpose(stoich_matrix)))

    # Constructing a new metabolite state vector.
    metState = np.zeros(len(np.transpose(stoich_matrix)))
    metState[ Currency + nutrientSet ] = 1

    # Getting copies of the relevant matrices.
    r = (np.copy(rho).T * scopeRxns).T
    S = (np.copy(stoich_matrix).T * scopeRxns).T
    rMat = (np.copy(rxnMat).T * scopeRxns).T
    pMat = (np.copy(prodMat).T * scopeRxns).T

    # Figuring out which reactions can be performed at this step.
    procRxnVec = ((np.dot(rMat, metState != 0) - sumRxnVec) == 0) * 1

    # Continuing calculation till no more reactions can be performed
.    isChecked = np.zeros(len(rxnMat))

    # Computing which reaction gets what share of which reactants.
    mask = np.abs(np.sum(r, axis = 0)) != 0
    shareMatrix = np.zeros_like(S)
    shareMatrix[:, np.where(mask)[0]] = ((r * metState)[:, mask] / 
                                         np.abs(np.sum(r, axis = 0))[mask])
    shareMatrix[:, Currency] = -1

    while procRxnVec.any():
        # Initializing the product metabolite state vector.
        prodState = np.zeros(len(np.transpose(stoich_matrix)))

        # Updating states after all accomplishable reactions.
        for thisRxn in np.where(procRxnVec)[0]:
            # Checking if found a usable reactant.
            allowedRct = []
            isChecked[thisRxn] = 1

            # Getting the reactants and products of this reaction, except currency.
            rs, ps = np.where(rMat[thisRxn])[0], np.where(pMat[thisRxn])[0]
            reactants = [tR for tR in rs if tR not in Currency]
            products = [tP for tP in ps if tP not in Currency]

            # Checking for limiting reactants.
            for thisReactant in reactants:
                # Now assuming that all of this reactant is used up. Are others enough?
                if not isLimiting(thisReactant, thisRxn, shareMatrix[thisRxn], S, reactants):
                    allowedRct.append(thisReactant)
                    limRct = deepcopy(thisReactant)
                    break

            # If nothing is limiting, everything gets used.
            if not allowedRct:
                    limRct = giveLimitingCurrency(r, thisRxn)

            # Updating metabolite amounts post reaction.
            for thisMet in products:
                ratio = S[thisRxn, thisMet] / S[thisRxn, limRct]
                prodState[thisMet] += shareMatrix[thisRxn, limRct] * ratio

            mets = np.append(rs, ps)
            for thisMet in mets[np.where(np.in1d(mets, np.array(Core + Energy)))]:
                ratio = S[thisRxn, thisMet] / S[thisRxn, limRct]

                # Updating E and B if these metabolites are produced or consumed.
                if thisMet == 1:
                    runningE += shareMatrix[thisRxn, limRct] * ratio
                elif thisMet in [3, 4]:
                    runningE += shareMatrix[thisRxn, limRct] * ratio
                elif thisMet in Core:
                    if thisMet in ps:
                        isCoreProduced[Core] = 1
                    runningB += shareMatrix[thisRxn, limRct] * ratio

            # Updating metabolite amounts post reaction.
            for thisMet in reactants:
                ratio = S[thisRxn, thisMet] / S[thisRxn, limRct]
                shareMatrix[thisRxn, thisMet] -= shareMatrix[thisRxn, limRct] * ratio

        # Redistributing the produced metabolites among reactions by demand.
        r[np.where(isChecked)] = 0
        mask = np.abs(np.sum(r, axis = 0)) != 0
        shareMatrix[:, np.where(mask)[0]] = ((r * prodState)[:, mask] / 
                                             np.abs(np.sum(r, axis = 0))[mask])
        shareMatrix[:, Currency] = -1

        # Recalculating performable reactions.
        procRxnVec = ((np.dot(rMat, np.sum(shareMatrix, axis = 0) != 0) - sumRxnVec) == 0) * 1
        procRxnVec[np.where(isChecked)] = 0

    if isCoreProduced[Core].all():
        return runningE, runningB, True
    else:
        return -1.0, -1.0, False
