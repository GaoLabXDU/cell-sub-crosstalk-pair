'''
The codes in 'identifying_cellSubCrossTalkPairs_candidates()' and 'coupledNMF()' refers to Duren et al.
from https://github.com/SUwonglab/CoupledNMF

If you use our codes in your research, please considering citing:
    DUREN Z, CHEN X, ZAMANIGHOMI M, et al. Integrative analysis of single-cell genomics data by coupled nonnegative matrix factorizations[J]. Proceedings of the National Academy of Sciences of the United States of America, 2018, 115(30): 7723–7728.

'''

import os
from sklearn.decomposition import NMF
from numpy import linalg as LA
import itertools
import pandas as pd
from statsmodels.stats.weightstats import ttest_ind
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from sklearn.cluster import AgglomerativeClustering



def preprocessing_expression_matrix(LRTable, XL_ori, XR_ori):

    # preprocessing of expression matrix of sender cell type
    XL_ori = XL_ori[np.sum(XL_ori, axis=1) > 0]
    GL_ori = list(XL_ori.index)
    CL = list(XL_ori.columns)

    # preprocessing of expression matrix of receiver cell type
    XR_ori = XR_ori[np.sum(XR_ori, axis=1) > 0]
    GR_ori = list(XR_ori.index)
    CR = list(XR_ori.columns)

    # select ligand and receptor exist in matrix
    GL_cover1 = sorted(list(set(list(LRTable['Ligand'])) & set(GL_ori)))
    GR_cover1 = sorted(list(set(list(LRTable['Receptor'])) & set(GR_ori)))
    LRTableCover = LRTable[((LRTable['Ligand'].isin(GL_cover1)) & (LRTable['Receptor'].isin(GR_cover1)))]
    GL_cover2 = sorted(list(set(list(LRTableCover['Ligand'])) & set(GL_ori)))
    GR_cover2 = sorted(list(set(list(LRTableCover['Receptor'])) & set(GR_ori)))

    # filter gene
    GL = GL_cover2
    GR = GR_cover2
    XL = pd.DataFrame(XL_ori.loc[GL, :])
    XR = pd.DataFrame(XR_ori.loc[GR, :])

    # L1 norm for matrix
    XR = XR.div(XR.max(axis=0), axis='columns')
    XL = XL.div(XL.max(axis=0), axis='columns')
    XR = XR.fillna(0)
    XL = XL.fillna(0)
    return(XR, XL, GL, GR, CL, CR, LRTableCover)

def construct_ligand_receptor_matrix(LRTableCover, GL, GR):
    # construct ligand-receptor(A) matrix
    L2ID = dict(zip(list(GL), range(len(GL))))
    R2ID = dict(zip(list(GR), range(len(GR))))
    A = np.zeros((XR.shape[0], XL.shape[0]))
    for i in list(LRTableCover.index):
        tempTable = LRTableCover.loc[i]
        L_index = L2ID[tempTable['Ligand']]
        R_index = R2ID[tempTable['Receptor']]
        A[R_index, L_index] = 1

    # Laplacian normalization of A
    sA1 = np.diag(np.sqrt(1 / np.sum(A, axis=1)))
    sA1[sA1 == np.inf] = 0
    sA2 = np.diag(np.sqrt(1 / np.sum(A, axis=0)))
    sA2[sA2 == np.inf] = 0
    A = np.dot(np.dot(sA1, A), sA2)
    return(A)

def ttest(x,y):
    statistic, p, df = ttest_ind(x, y,alternative='larger')
    return(p)

def coupledNMF(XL, XR, A, WL0, WR0, HL0, HR0, lambda1, lambda2, mu, K):
    '''
    The codes in this function refer to Duren et al. from https://github.com/SUwonglab/CoupledNMF
    '''
    eps = 0.001
    WL = WL0
    WR = WR0
    HL = HL0
    HR = HR0

    # Initialization of coupledNMF
    maxiter = 500
    err = 1
    terms = np.zeros(maxiter)
    it = 0

    terms[it] = 0.5 * np.power(np.linalg.norm(XL - np.dot(WL, HL), ord='fro'), 2) \
                + 0.5 * lambda1 * np.power(np.linalg.norm(XR - np.dot(WR, HR), ord='fro'), 2) \
                - lambda2 * np.power(np.trace(np.dot(np.dot(np.transpose(WR), A), WL)), 2) \
                + mu * (np.power(np.linalg.norm(WL, ord='fro'), 2) + np.power(np.linalg.norm(WR, ord='fro'), 2))

    while it < maxiter - 1 and err > 1e-6:
        it = it + 1
        T1 = 0.5 * lambda2 * np.dot(np.transpose(A), WR)
        T1[T1 < 0] = 0
        WL = WL * (np.dot(XL, np.transpose(HL)) + T1) / (
                eps + np.dot(WL, np.dot(HL, np.transpose(HL))) + 2 * mu * WL)
        HL = HL * (np.dot(np.transpose(WL), XL)) / (eps + np.dot(np.dot(np.transpose(WL), WL), HL))
        T2 = 0.5 * (lambda2 / lambda1 + eps) * np.dot(A, WL)
        T2[T2 < 0] = 0
        WR = WR * (np.dot(XR, np.transpose(HR)) + T2) / (
                eps + np.dot(WR, np.dot(HR, np.transpose(HR))) + 2 * mu * WR)
        HR = HR * (np.dot(np.transpose(WR), XR) / (eps + np.dot(np.dot(np.transpose(WR), WR), HR)))
        m1 = np.zeros((K, K))
        m2 = np.zeros((K, K))
        for z in range(K):
            m1[z, z] = np.linalg.norm(HL[z, :])
            m2[z, z] = np.linalg.norm(HR[z, :])

        WR = np.dot(WR, m2)
        WL = np.dot(WL, m1)
        HL = np.dot(np.linalg.inv(m1), HL)
        HR = np.dot(np.linalg.inv(m2), HR)

        terms[it] = 0.5 * np.power(np.linalg.norm(XL - np.dot(WL, HL), ord='fro'), 2) \
                    + 0.5 * lambda1 * np.power(np.linalg.norm(XR - np.dot(WR, HR), ord='fro'), 2) \
                    - lambda2 * np.power(np.trace(np.dot(np.dot(np.transpose(WR), A), WL)), 2) \
                    + mu * (np.power(np.linalg.norm(WL, ord='fro'), 2) + np.power(np.linalg.norm(WR, ord='fro'), 2))
        err = np.abs(terms[it] - terms[it - 1]) / np.abs(terms[it - 1])
    SL = np.argmax(HL, 0)
    SR = np.argmax(HR, 0)
    return(SL,SR)

def identifying_cellSubCrossTalkPairs_candidates(XL, XR, A, K):
    '''
    The codes in this function refer to Duren et al. from https://github.com/SUwonglab/CoupledNMF
    '''
    rep = 50
    list_SL = []
    list_SR = []
    for c_i in [0,1,2,3,4]:

        # Initialization of WL and HL
        errL = np.zeros(rep)
        for i in range(0, rep):
            model = NMF(n_components=K, init='random', random_state=i + rep * c_i, solver='cd', max_iter=50)
            WL0 = model.fit_transform(XL)
            HL0 = model.components_
            errL[i] = LA.norm(XL - np.dot(WL0, HL0), ord='fro')
        model = NMF(n_components=K, init='random', random_state=np.argmin(errL), solver='cd', max_iter=1000)
        WL0 = model.fit_transform(XL)
        HL0 = model.components_
        SL0 = np.argmax(HL0, 0) # index matrix: each cell belong to which cluster

        # Initialization of WR and HR
        errR = np.zeros(rep)
        for i in range(0, rep):
            model = NMF(n_components=K, init='random', random_state=i + rep * c_i, solver='cd', max_iter=50)
            WR0 = model.fit_transform(XR)
            HR0 = model.components_
            errR[i] = LA.norm(XR - np.dot(WR0, HR0), ord='fro')
        model = NMF(n_components=K, init='random', random_state=np.argmin(errR), solver='cd', max_iter=1000)
        WR0 = model.fit_transform(XR)
        HR0 = model.components_
        SR0 = np.argmax(HR0, 0)

        # Adjust initialization
        pL = [[ttest(XL.iloc[j, SL0 == i], XL.iloc[j, SL0 != i]) for j in range(XL.shape[0])] for i in range(K)]
        pL = np.array(np.vstack(pL)).T
        pR = [[ttest(XR.iloc[j, SR0 == i], XR.iloc[j, SR0 != i]) for j in range(XR.shape[0])] for i in range(K)]
        pR = np.array(np.vstack(pR)).T

        WPL = np.zeros((WL0.shape))
        pL[np.isnan(pL)] = 1
        scores = -np.log10(pL)
        temp = int(len(GL) / 2)
        for i in range(K):
            indexs = scores[:, i].argsort()[-temp:][::-1]
            WPL[indexs, i] = 1

        WPR = np.zeros((WR0.shape))  # the 0,1 weigted matrix based on top max value (set to 1) in 'scores'
        pR[np.isnan(pR)] = 1
        scores = -np.log10(pR)
        temp = int(len(GR) / 2)
        for i in range(K):
            indexs = scores[:, i].argsort()[-temp:][::-1]
            WPR[indexs, i] = 1

        perm = list(itertools.permutations(range(K)))  # all permutation and combination of 0:K
        score = np.zeros(len(perm))
        for i in range(len(perm)):
            score[i] = np.trace(np.dot(np.dot(np.transpose(WPR[:,perm[i]]), A), WPL))

        match = np.argmax(score)
        WR0 = WR0[:, perm[match]]
        HR0 = HR0[perm[match], :]

        # select parameter
        lambda1 = pow(LA.norm(XL - np.dot(WL0, HL0), ord='fro'), 2) / \
                  pow(LA.norm(XR - np.dot(WR0, HR0), ord='fro'), 2)
        lambda2 = 1000 * 0.5 * pow(LA.norm(XL - np.dot(WL0, HL0), ord='fro'), 2) / \
                  pow(np.trace(np.dot(np.dot(np.transpose(WR0), A), WL0)), 2)
        mu = 0.5 * pow(LA.norm(XL - np.dot(WL0, HL0), ord='fro'), 2) / \
             (pow(LA.norm(WL0, ord='fro'), 2) + pow(LA.norm(WR0, ord='fro'), 2))

        try:
            SL,SR = coupledNMF(XL, XR, A, WL0, WR0, HL0, HR0, lambda1, lambda2, mu, K)
            list_SL.append(SL)
            list_SR.append(SR)
        except:
            continue
    return(list_SL,list_SR)

def construct_coupled_concensus_matrix(list_SL,list_SR):
    cMat = 0
    num_c = len(list_SL)
    for c_i in range(num_c):
        SL = list_SL[c_i].tolist()
        numCell1 = len(SL)
        SR = list_SR[c_i].tolist()
        numCell2 = len(SR)
        indexMat = np.zeros((numCell1 + numCell2, numCell1 + numCell2))
        S = SL + SR
        for i in range(0, numCell1 + numCell2):
            for j in range(i, numCell1 + numCell2):
                if S[i] == S[j]:
                    indexMat[i][j] = 1
                    indexMat[j][i] = 1
        cMat = cMat + indexMat
    return(cMat)

def hierarchical_clustering(cMat,K):
    distMat = 1 - cMat / np.max(cMat)
    hiCluster = AgglomerativeClustering(n_clusters=K).fit(distMat)
    LabelsHiCluster = hiCluster.labels_
    return(LabelsHiCluster)

def create_cellSubCrossTalkPairs_label(CL, CR, cluster):
    cellSubCrossTalkPairsLabel = pd.DataFrame({'cell':CL+CR,
                                               'cellType':['sender']*len(CL) + ['receiver']*len(CR),
                                               'cellSubCrossTalkPairID':cluster})
    return(cellSubCrossTalkPairsLabel)

if __name__ == '__main__':

    path = '..//Demo_data'

    ## load LR
    LRFileName = 'PairsLigRec_simple_mouse.txt'
    LRTable = pd.read_csv(os.path.join(path, LRFileName), sep='\t')

    ## load expression matrix(XL) of sender cell type. 'X' represents matrix, 'G' represents genes.
    fName_Lcell = 'mat_F002_Endo'
    XL_ori = pd.read_csv(os.path.join(path, fName_Lcell), sep=',', index_col=0)

    ## load expression matrix(XR) of receiver cell type. 'X' represents matrix, 'G' represents genes.
    fName_Rcell = 'mat_F002_Micro'
    XR_ori = pd.read_csv(os.path.join(path, fName_Rcell), sep=',', index_col=0)

    ## preprocessing data
    # preprocessing expression matrix
    XR, XL, GL, GR, CL, CR, LRTableCover = preprocessing_expression_matrix(LRTable, XL_ori, XR_ori)
    # construct ligand-receptor(A) matrix
    A = construct_ligand_receptor_matrix(LRTableCover, GL, GR)

    ## identifying candidates of cell sub-crosstalk Pairs
    K=2
    list_SL, list_SR = identifying_cellSubCrossTalkPairs_candidates(XL, XR, A, K)

    ## merging candidates by coupled consensus clustering
    # constructing coupled concensus matrix
    cMat = construct_coupled_concensus_matrix(list_SL, list_SR)
    # hierarchical clustering
    cluster = hierarchical_clustering(cMat,K)

    ## create cellSubCrossTalkPairs label
    cellSubCrossTalkPairsLabel = create_cellSubCrossTalkPairs_label(CL, CR, cluster)










