# -*- coding: utf-8 -*-
''' Collection of pre-formatted plot functions

'''

__author__ = "Christian Brinch"
__copyright__ = "Copyright 2018"
__credits__ = ["Christian Brinch"]
__license__ = "AFL 3.0"
__version__ = "0.1"
__maintainer__ = "Christian Brinch"
__email__ = "cbri@gfood.dtu.dk"

from operator import add, sub
from pylab import percentile
import scipy.stats as ss
import numpy as np


def histogram(axis, data, color='black', density=True, bins=None):
    ''' Return histogram using Freedman-Diaconis rule for bin size '''
    if bins is None:
        bin_width = 2.*ss.iqr(data)/pow(len(data), 1./3.)
        nbins = int(np.ceil((max(data)-min(data))/bin_width))
        axis.hist(data, nbins, density=density, facecolor=color, alpha=0.2, label='')
        return axis.hist(data, nbins, density=density, edgecolor=color,
                         facecolor='none', alpha=0.8, label='')
    axis.hist(data, bins, density=density, facecolor=color, alpha=0.2, label='')
    return axis.hist(data, bins, density=density, edgecolor=color, facecolor='none',
                     alpha=0.8, label='')


def plot_fits(func, axis, xvar, popt, pcov, color='grey'):
    ''' Plot the function f with options popt '''
    axis.plot(xvar, func(xvar, *popt), color='grey')
    xsample = np.random.multivariate_normal(popt, pcov, 10000)
    ysample = np.asarray([func(xvar, *pi) for pi in xsample])
    lower = percentile(ysample, 1., axis=0)
    upper = percentile(ysample, 99., axis=0)
    axis.plot(xvar, lower, color=color, alpha=0.7)
    axis.plot(xvar, upper, color=color, alpha=0.7)
    axis.fill_between(xvar, lower, upper, color=color, alpha=0.5)
