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
import scipy.stats as ss


def alr_clr(vector, kind):
    ''' Log-transform vector based on kind '''
    if kind == 'alr':
        return [np.log10(i/i[-1]) for i in vector]
    elif kind == 'clr':
        return [np.log10(i) - np.mean(np.log10(i)) for i in vector]
    return None


def transform(frame, normfields=None, bayesian=True, n_samples=1000):
    ''' This function will perform a log-transform depending on the
        normalization. If bayesian is true, bayesian sampling will be applied
        to predict the finite probability for zero entries in frameself.
    '''
    prob = frame.copy(deep=True)
    error = frame.copy(deep=True)
    if normfields is not None:
        prob.loc['norm'] = 0.
        error.loc['norm'] = 0.
        for column in normfields:
            prob.loc['norm'][column] = np.sum(normfields[column])
        kind = 'alr'
    else:
        kind = 'clr'

    for column in frame:
        if bayesian is True:
            p_matrix = ss.dirichlet.rvs(np.array(prob[column])+0.5, n_samples)
        else:
            p_matrix = [np.array(prob[column])]

        c_matrix = alr_clr(p_matrix, kind)
        prob[column] = [np.mean(i) for i in zip(*c_matrix)]
        error[column] = [np.std(i) for i in zip(*c_matrix)]
    if normfields is not None:
        prob = prob.drop('norm')
        error = error.drop('norm')

    return prob, error


def geom_mean(composition, n_samples=1000):
    ''' Calculate the geometric mean of a composition '''
    for column in composition:
        p_matrix = ss.dirichlet.rvs(np.array(composition[column])+0.5,
                                    n_samples)
        composition[column] = np.mean(p_matrix, axis=0)

    return [np.exp(np.mean(np.log(x))) for x in np.array(composition)]


def totvar(composition, mean):
    ''' Calculate the total variation of a composition
        TODO: This is very slow at the moment. Try to optimize
    '''
    variance = 0.
    dim = np.shape(composition)
    for vector in np.array(composition.T):
        dist_a = np.sqrt(1./(2.*dim[0])
                         * np.sum([(np.log(vector[i]/vector[j])
                                    - np.log(mean[i]/mean[j]))**2
                                   for i in range(dim[0])
                                   for j in range(dim[0])]))
        variance += 1./dim[1] * dist_a**2
    return variance
