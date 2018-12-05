# -*- coding: utf-8 -*-
''' Functions to extract diversity information from composition
'''

__author__ = "Christian Brinch"
__copyright__ = "Copyright 2018"
__credits__ = ["Christian Brinch"]
__license__ = "AFL 3.0"
__version__ = "0.1"
__maintainer__ = "Christian Brinch"
__email__ = "cbri@gfood.dtu.dk"

import numpy as np


def richness(composition):
    ''' Calculate the richness of composition '''
    rich = []
    for column in composition:
        rich.append(np.sum([1. for i in composition[column] if i > 0.]))
    return rich


def shannon(composition):
    ''' Calculate the Shannon index of a composition '''
    h_idx = []
    for column in composition:
        h_idx.append(-np.sum([i*np.log10(i)
                              for i in composition[column]/np.sum(composition[column])]))
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
