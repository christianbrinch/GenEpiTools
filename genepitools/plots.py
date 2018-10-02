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
