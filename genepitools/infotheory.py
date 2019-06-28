# -*- coding: utf-8 -*-
''' Functions to extract diversity information from composition
'''

__author__ = "Christian Brinch"
__copyright__ = "Copyright 2018"
__credits__ = ["Christian Brinch"]
__license__ = "AFL 3.0"
__version__ = "0.1"
__maintainer__ = "Christian Brinch"
__email__ = "cbri@food.dtu.dk"

import numpy as np
import pandas as pd


def richness(composition):
    ''' Calculate the richness of composition '''
    return composition.astype(bool).sum(axis=0)


def shannon(composition):
    ''' Calculate the Shannon index of a composition '''
    p = composition/composition.sum(axis=0)
    h_idx = -(np.log(p+1e-30)*p).sum(axis=0)
    return h_idx


def chao(composition):
    ''' Chao1 richness estimator for integer composition '''
    s_obs = richness(composition)
    s_chao1 = []
    for column in composition:
        f_1 = np.sum([1. for i in composition[column] if i == 1])
        f_2 = np.sum([1. for i in composition[column] if i == 2])
        if f_2 > 0:
            s_chao1.append(s_obs + (f_1**2/(2.*f_2)))
        else:
            s_chao1.append(s_obs + f_1*(f_1-1.)/2.)
    return s_chao1


def rarefy(lmappings, lcounts):
    ''' rarefy composition to lowest member
        TODO: specify rarefaction level and discard samples below
    '''

    level = min(lmappings.loc['notPhiX'])

    for dba in lmappings.T:
        if dba != 'notPhiX':
            lmappings.T[dba] = lmappings.T[dba].multiply(
                level/lmappings.T['notPhiX']).astype(int)

    for sample in lcounts:
        if lmappings[sample].loc['notPhiX'] > level:
            tmp = lcounts[sample].to_frame()

            weighted_pop = [lcounts[sample].index[elem]
                            for elem, cnt in enumerate(lcounts[sample]) for i in range(cnt)]

            for repeats in range(3):
                subset = np.random.choice(
                    weighted_pop, size=lmappings[sample].loc['ResFinder'], replace=False)

                subdict = dict(np.asarray(np.unique(subset, return_counts=True)).T)
                new = pd.DataFrame.from_dict(subdict, orient='index', columns=['r'+str(repeats)])

                tmp = tmp.merge(new, how='outer', left_index=True, right_index=True).fillna(0)

            tmp = tmp.drop(labels=sample, axis=1)
            tmp['median'] = tmp.median(axis=1, skipna=False).astype(int)

            lcounts[sample] = tmp['r0'].astype(float)

    return lmappings, lcounts
