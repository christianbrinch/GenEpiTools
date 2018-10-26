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


from pylab import percentile
import scipy.stats as ss
import numpy as np


def histogram(axis, data, color='black', density=True, bins=None):
    ''' Return histogram using Freedman-Diaconis rule for bin size '''
    if bins is None:
        bin_width = 2.*ss.iqr(data)/pow(len(data), 1./3.)
        nbins = int(np.ceil((max(data)-min(data))/bin_width))
        axis.hist(data, nbins, density=density, facecolor=color, alpha=0.3, label='')
        return axis.hist(data, nbins, density=density, edgecolor='white',
                         facecolor='none', alpha=1.0, label='', lw=1.0)
    axis.hist(data, bins, density=density, facecolor=color, alpha=0.3, label='')
    return axis.hist(data, bins, density=density, edgecolor='white', facecolor='none',
                     alpha=1.0, lw=1.0, label='')


def plot_fits(func, axis, xvar, popt, pcov, color='grey'):
    ''' Plot the function f with options popt '''
    axis.plot(xvar, func(xvar, *popt), color=color, lw=2, alpha=0.8)
    xsample = np.random.multivariate_normal(popt, pcov, 50000)
    ysample = np.asarray([func(xvar, *pi) for pi in xsample])
    lower = percentile(ysample, 0.3, axis=0)
    upper = percentile(ysample, 99.7, axis=0)
    axis.plot(xvar, lower, color=color, alpha=0.4, lw=1.0)
    axis.plot(xvar, upper, color=color, alpha=0.4, lw=1.0)
    axis.fill_between(xvar, lower, upper, color=color, alpha=0.2)


def errorbar(axis, xvar, yvar, yerr, color='grey', fmt='o'):
    ''' Plot errorbars '''
    _, caps, bars = axis.errorbar(xvar, yvar, yerr=yerr, fmt=fmt, color=color,
                                  alpha=0.7, capsize=4., capthick=2.0,
                                  markersize=6., markeredgewidth=0.5, markeredgecolor='white')
    _ = [abar.set_alpha(0.6) for abar in bars]
    _ = [acap.set_alpha(0.6) for acap in caps]
