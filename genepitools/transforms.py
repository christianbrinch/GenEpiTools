# -*- coding: utf-8 -*-
""" Functions for log-ratio transforming data

"""

__author__ = "Christian Brinch"
__copyright__ = "Copyright 2018"
__credits__ = ["Christian Brinch"]
__license__ = "AFL 3.0"
__version__ = "0.1"
__maintainer__ = "Christian Brinch"
__email__ = "cbri@food.dtu.dk"

import numpy as np
import pandas as pd
import scipy.stats as ss


def transform(dataframe, transformation, n_samples=1000):
    ''' Wrapper function for alr/clr transformation of counts '''
    new_dataframe = pd.DataFrame(index=dataframe.index)
    new_errorframe = pd.DataFrame(index=dataframe.index)
    if transformation == 'alr':
        for column in dataframe.columns:
            new_dataframe[column], new_errorframe[column] = alr(list(dataframe[column]), n_samples)
        new_dataframe.drop(new_dataframe.index[-1], inplace=True)
        new_errorframe.drop(new_errorframe.index[-1], inplace=True)
    elif transformation == 'clr':
        for column in dataframe.columns:
            new_dataframe[column], new_errorframe[column] = clr(list(dataframe[column]), n_samples)
    else:
        print "Unknown transformation type"
        return None
    return new_dataframe, new_errorframe


def alr(data, n_samples):
    ''' Additive log ratio transformation
        By convention, data is normalized by last entry in data
    '''
    p_matrix = ss.dirichlet.rvs(np.array(data)+0.5, n_samples)
    a_matrix = [np.log(i/i[-1]) for i in p_matrix]
    values = [np.mean(i) for i in zip(*a_matrix)]
    errors = [np.std(i) for i in zip(*a_matrix)]
    return values, errors


def clr(data, n_samples):
    ''' Centered log ratio transformation
        data is normalized by geometric mean
    '''
    p_matrix = ss.dirichlet.rvs(np.array(data)+0.5, n_samples)
    c_matrix = [np.log(i) - np.mean(np.log(i)) for i in p_matrix]
    values = [np.mean(i) for i in zip(*c_matrix)]
    errors = [np.std(i) for i in zip(*c_matrix)]
    return values, errors


def sum_amr_and_bac(counts, mappings):
    ''' Return total AMR and total bacteria count
        TODO: bacteria count is no good - use bacteria genes instead
    '''
    counts = counts.agg(['sum'], axis=0)
    counts.loc['bac'] = 0.
    for index, _ in mappings.iterrows():
        counts.loc['bac'][index] = mappings.loc[index]['Bacteria']
    return counts
