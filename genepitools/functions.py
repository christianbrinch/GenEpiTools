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

import numpy as np


def gauss_func(xvar, mean, sigma, amp=1., offset=0.):
    ''' A 1D Gaussian function '''
    return amp / np.sqrt(2.*np.pi*sigma**2) * np.exp(-pow(xvar-mean, 2)/(2*sigma**2)) + offset


def sine_wave(time, amp, omega, phase, ang=True):
    ''' A standard sine wave function '''
    if ang is False:
        omega = 2.*np.pi*omega
    return amp*np.sin(omega*time+phase)
