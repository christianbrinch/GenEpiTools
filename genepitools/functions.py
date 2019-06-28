# -*- coding: utf-8 -*-
''' Collection of mathematical functions

'''

__author__ = "Christian Brinch"
__copyright__ = "Copyright 2018"
__credits__ = ["Christian Brinch"]
__license__ = "AFL 3.0"
__version__ = "0.1"
__maintainer__ = "Christian Brinch"
__email__ = "cbri@food.dtu.dk"

import numpy as np
from scipy.special import factorial


def gauss_func(xvar, mean, sigma, amp=1., offset=0.):
    ''' A 1D Gaussian function '''
    # return amp / np.sqrt(2.*np.pi*sigma**2) * np.exp(-pow(xvar-mean, 2)
    #    /(2*sigma**2)) + offset
    return amp * np.exp(-pow(xvar-mean, 2)/(2*sigma**2)) + offset


def poisson(k, lamb):
    ''' A Poisson distribution '''
    return (lamb**k/factorial(k)) * np.exp(-lamb)


def erlang(xvar, k, lamb):
    ''' An Erlang distribution'''
    return 1./factorial(k-1) * pow(lamb, k)*pow(xvar, k-1)*np.exp(-lamb*(xvar))


def sine_wave(time, amp, omega, phase, offset=0.):
    ''' A standard sine wave function '''
    return amp*np.sin(omega*time+phase)+offset


def line(xvar, slope, intercept):
    ''' A straight line '''
    return xvar*slope+intercept


def sine_plus_line_wrapper(xvar, *opt):
    ''' A wrapper function for the sum of a sine and a line '''
    return sine_wave(xvar, opt[0], opt[1], opt[2], opt[3])+line(xvar, opt[4],
                                                                opt[5])


def powerlaw(xvar, scale, exp):
    ''' A generic powerlaw '''
    return scale*pow(xvar, exp)
