# -*- coding: utf-8 -*-
'''
    Python implementation of ALDEx2
'''

__author__ = "Christian Brinch"
__copyright__ = "Copyright 2019"
__credits__ = ["Christian Brinch"]
__license__ = "AFL 3.0"
__version__ = "0.1"
__maintainer__ = "Christian Brinch"
__email__ = "cbri@food.dtu.dk"

import random
import numpy as np
import pandas as pd
import scipy.stats as ss
from statsmodels.stats.multitest import multipletests


def anova_denom(data, denom='all'):
    ''' Subset the feature list according to denom criteria '''
    if denom == 'all':
        return np.arange(len(data))
    if denom == 'iqlr':
        clr = np.log2(data+0.5)-np.mean(np.log2(data+0.5))
        var = clr.var(axis=1)
        qnt = var.quantile([0.25, 0.75])
        return np.where((var > qnt.iloc[0]) & (var < qnt.iloc[1]))[0]
    if denom == 'lvha':
        clr = np.log2(data+0.5)-np.mean(np.log2(data+0.5))
        med = clr.median(axis=1)
        var = clr.var(axis=1)
        idx = np.where(med > 0.)
        return np.where((med > 0.) & (var.iloc[idx] < var.iloc[idx].quantile(0.25)))
    print("Unknown denominator. Falling back to all")
    return np.arange(len(data))


def anova_clr(data, nsamples, denom):
    '''
    Calculate the CLR transform and store all Monte Carlo samples. Pre-allocate
    clr array for faster access.
    '''
    clr = [[] for _ in range(np.shape(data)[1])]
    for idx, column in enumerate(data):
        p_matrix = ss.dirichlet.rvs(data[column]+0.5, nsamples)
        clr[idx] = np.log2(p_matrix) - np.mean(np.log2(p_matrix[:, denom]))
    return np.array(clr)


def anova_ttest(clr, subsets, nsamples):
    '''
    Perform Welchs T-test for each Monte Carlo sample and calculate relative
    abundances and median values
    '''
    p_val = [ss.ttest_ind(clr[:len(subsets[0]), i], clr[len(subsets[0]):, i], equal_var=False)[1]
             for i in range(nsamples)]
    rab = np.mean(p_val, axis=0)
    rab = multipletests(rab, alpha=0.003, method='fdr_bh')[1]
    median_all = np.median(clr, axis=(0, 1))
    return rab, median_all


def anova_within(clrsets, setlengths):
    '''
    Calculate the variation within groups and select the group
    with highest variation
    '''
    win = []
    for group in range(2):
        rand_seq = [random.sample(range(0, setlengths[group]), setlengths[group]) for _ in range(2)]
        win.append([np.abs(clrsets[group][rand_seq[0][i]] - clrsets[group][rand_seq[1][i]])
                    for i in range(setlengths[group])])
    if np.sum(win[0]) > np.sum(win[1]):
        return np.median(win[0], axis=0)
    return np.median(win[1], axis=0)


def anova_between(clrsets, setlengths):
    ''' Calculate the variation between groups '''
    shortset = setlengths.index(min(setlengths))
    rand_seq = [random.sample(range(0, setlengths[shortset]), setlengths[shortset])
                for _ in range(2)]
    btw = [clrsets[1][rand_seq[0][i]]-clrsets[0][rand_seq[1][i]]
           for i in range(setlengths[shortset])]
    return np.median(btw, axis=0)


def comp_anova(data, subsets, nsamples, denom):
    ''' CLR-transform, T-test, and variation calculations '''
    features = anova_denom(data, denom)
    clr = anova_clr(data, nsamples, features)
    rab, median_all = anova_ttest(clr, subsets, nsamples)
    setlengths = [len(i)*nsamples for i in subsets]
    clrsets = [clr[:len(subsets[0])].reshape((setlengths[0], len(data))),
               clr[len(subsets[0]):].reshape((setlengths[1], len(data)))]
    diff_win = anova_within(clrsets, setlengths)
    diff_btw = anova_between(clrsets, setlengths)
    return median_all, diff_win, diff_btw, rab


def anova_test():
    ''' Test against the original ALDeX2 R implementation '''
    data = pd.read_csv('selex.txt', index_col=0, encoding='utf-8', delim_whitespace=True)
    data = data[:400]
    subsets = [['1_ANS', '1_BNS', '1_CNS', '1_DNS', '2_ANS', '2_CNS', '2_DNS'],
               ['1_AS', '1_BS', '1_CS', '1_DS', '2_AS', '2_CS', '2_DS']]
    return comp_anova(data, subsets, 500, 'all')
