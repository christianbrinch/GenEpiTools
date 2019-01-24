# -*- coding: utf-8 -*-
""" Quality assesment for metagenomic mappings
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
import scipy.optimize as so
import scipy.stats as ss
import seaborn as sns
import scipy.interpolate as ip
import matplotlib.pyplot as plt
import genepitools.functions as gfunc
import genepitools.plots as gplots

sns.set()


def validate(composition):
    ''' Plot sorted loglog distributions of each component '''
    dim = np.shape(composition)

    composition = composition.div(composition.sum(axis=1), axis=0)

    _ = plt.figure(1)
    plt.clf()
    axis = plt.subplot(111)
    #axis.set_xlim(-0.25, 85)
    # axis.set_ylim(1e-6, 2e0)
    axis.set_ylabel('Percentage counts')
    axis.set_xlabel('Sorted sample #')
    xvar = np.linspace(1, dim[0], dim[0])

    for idx, column in enumerate(composition.columns):

        color = 'grey'
        if column in ['ResFinder', 'Bacteria', 'Bacteria_draft',
                      'HumanMicrobiome', 'Plasmid']:
            color = 'faded green'
        if column in ['Unmapped']:
            color = 'black'
        if column in ['Vertebrates_mammals', 'Vertebrates_other', 'Human',
                      'Mitochondrion']:
            color = 'denim blue'
        if column in ['Invertebrates', 'Plant', 'Virus']:
            color = 'amber'
        if column in ['MetalResistance']:
            color = 'pale red'

#        axis.loglog(xvar, composition[column].sort_values()[::-1],
#                    c=sns.xkcd_rgb[color], lw=3, alpha=0.8)
        axis.plot(xvar, np.log10(composition[column].sort_values()[::-1]+1e-10),
                  c=sns.xkcd_rgb[color], lw=3, alpha=0.8)

        axis.annotate(column, (-0.1, np.log10(composition[column].sort_values()[-1])+1e-10),
                      ha='right', fontsize=7)


def residuals(composition, color='grey', norm=7.0, overplot=False):
    dbase = ['Bacteria', 'Bacteria_draft', 'ResFinder', 'HumanMicrobiome',
             'Plasmid', 'Virus', 'Unmapped',
             'Plant', 'Invertebrates', 'Mitochondrion', 'Vertebrates_mammals',
             'Vertebrates_other', 'Human']

    for f in [0, 2]:
        fig = plt.figure(59+f)
        if not overplot:
            plt.clf()
        sns.set(style="darkgrid")
        sns.set(font_scale=1.05)
        # Vector components are [x_pos, y_pos, width, height]
        axt = [0.12, 0.84, 0.84, 0.13]
        axm = [0.12, 0.10, 0.84, 0.72]

        axis = (
            fig.add_axes(axt, frame_on=True),
            fig.add_axes(axm, frame_on=True))

        axis[0].set_yticks([])
        axis[0].set_ylim(-.5, .5)
        axis[0].set_xlim(5.6, 8.3)
        axis[1].set_xlim(5.6, 8.3)
        axis[1].set_xlabel('Log(Total reads)')
        axis[1].set_ylabel('Log(Reads)')

        dat_x = np.log10(composition.loc['notPhiX'])
        if f == 0:
            dat_y = np.log10(composition.loc['Bacteria']+composition.loc['Bacteria_draft'])
        else:
            dat_y = np.log10(composition.loc[dbase[f]])
        mean_x = np.mean(dat_x)
        mean_y = np.mean(dat_y)

        xvar = np.array([min(dat_x), max(dat_x)])

        popt, pcov = so.curve_fit(gfunc.line, dat_x, dat_y)
        gplots.plot_fits(gfunc.line, axis[1], xvar, popt, pcov, color=color)
        bin_width = 2.*ss.iqr(dat_x)/pow(len(dat_x), 1./3.) - 0.03
        nbins = int(np.ceil((max(dat_x)-min(dat_x))/bin_width))

        qua1 = []
        qua3 = []
        mi = []
        ma = []
        for i in range(nbins):
            subset = [idx for idx, j in enumerate(dat_x) if j > (
                min(dat_x)+i*bin_width) and j <= (min(dat_x)+(i+1)*bin_width)]

            mdn = np.median(dat_y[subset]-dat_x[subset]+mean_x-mean_y)

            if len(subset) > 1:
                qua1.append(np.median([k for k in dat_y[subset] -
                                       dat_x[subset]+mean_x-mean_y if k < mdn]))
                qua3.append(np.median([k for k in dat_y[subset] -
                                       dat_x[subset]+mean_x-mean_y if k > mdn]))
            else:
                qua1.append(0)
                qua3.append(0)
            mi.append(-np.min(dat_y[subset]-dat_x[subset]+mean_x-mean_y))
            ma.append(np.max(dat_y[subset]-dat_x[subset]+mean_x-mean_y))

        bins = np.linspace(min(dat_x), max(dat_x), nbins)
        xvar = np.linspace(min(dat_x), max(dat_x), 500)
        spl1 = ip.interp1d(bins, qua1, kind='cubic')
        spl2 = ip.interp1d(bins, qua3, kind='cubic')

        axis[1].scatter(dat_x, dat_y, alpha=0.7, s=15, color=color)
        axis[1].plot(xvar, (xvar-mean_x)+mean_y, '--', color='black', lw=2, alpha=0.5)

        axis[0].errorbar(bins, np.zeros(len(bins)), yerr=[mi, ma], alpha=0.7,
                         capsize=4., capthick=2.0, markersize=6., markeredgewidth=0.5,
                         fmt='')

        axis[0].plot(xvar, spl1(xvar), color=color)
        axis[0].plot(xvar, spl2(xvar), color=color)
        axis[0].fill_between(xvar, spl1(xvar), spl2(xvar), color=color, alpha=0.4)
        axis[0].plot([min(xvar), max(xvar)], [0, 0], '--', color='black')

        for i, name in enumerate(composition.columns):
            axis[1].annotate(name, (dat_x[i], dat_y[i]))

        # axis[1].annotate()
        if not overplot:
            axis[1].annotate(dbase[f], (6.1, axis[1].get_ylim()[1]
                                        - .1*(axis[1].get_ylim()[1]
                                              - axis[1].get_ylim()[0])))
