# -*- coding: utf-8 -*-
''' Collection of mathematical functions

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
import matplotlib.pyplot as plot


def histogram(ax, data, color='black', density=True):
    fd = 2.*ss.iqr(data)/pow(len(data), 1./3.)
    nbins = int(np.ceil((max(data)-min(data))/fd))
    ax.hist(data, nbins, density=density, facecolor=color, alpha=0.4)
    return ax.hist(data, nbins, density=density, edgecolor=color, facecolor='none', alpha=0.8)
