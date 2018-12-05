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
from matplotlib.patches import Ellipse
import scipy.stats as ss
import numpy as np


def histogram(axis, data, color='black', density=True, bins=None):
    ''' Return histogram using Freedman-Diaconis rule for bin size '''
    if bins is None:
        bin_width = 2.*ss.iqr(data)/pow(len(data), 1./3.)
        nbins = int(np.ceil((max(data)-min(data))/bin_width))
        axis.hist(data, nbins, density=density, facecolor=color, alpha=0.3, label='')
        entries, edges, _ = axis.hist(data, nbins, density=density, edgecolor='white',
                                      facecolor='none', alpha=1.0, label='', lw=1.0)
        return entries, 0.5*(edges[1:] + edges[:-1])

    axis.hist(data, bins, density=density, facecolor=color, alpha=0.3, label='')
    entries, edges, _ = axis.hist(data, bins, density=density, edgecolor='white', facecolor='none',
                                  alpha=1.0, lw=1.0, label='')

    return entries, 0.5*(edges[1:] + edges[:-1])


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


def scatter_with_err_ellipse(axis, frame, loc, color='black', name=None):
    ''' Make scatter plot with covariance ellipse. Good for biplot '''
    axis.scatter(-frame.loc[loc[0]],
                 frame.loc[loc[1]],
                 marker='.',
                 color=color,
                 alpha=0.3,
                 s=30,
                 label=None)
    cov = np.cov(-frame.loc[loc[0]],
                 frame.loc[loc[1]])
    lambda_, angle = np.linalg.eig(cov)
    lambda_ = np.sqrt(lambda_)
    ell = Ellipse(xy=(np.mean(-frame.loc[loc[0]]),
                      np.mean(frame.loc[loc[1]])),
                  width=lambda_[0]*1.*2,
                  height=lambda_[1]*1.*2,
                  angle=np.rad2deg(np.arccos(angle[0, 0])),
                  alpha=0.6,
                  edgecolor=color,
                  fill=False,
                  lw=3)
    axis.add_artist(ell)
    axis.scatter(np.mean(-frame.loc[loc[0]]),
                 np.mean(frame.loc[loc[1]]),
                 marker='x',
                 s=80,
                 color=color,
                 label=name)
