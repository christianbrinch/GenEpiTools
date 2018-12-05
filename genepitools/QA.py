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
import seaborn as sns
import matplotlib.pyplot as plt

sns.set()


def validate(composition):
    ''' Plot sorted loglog distributions of each component '''
    dim = np.shape(composition)

    composition = composition.div(composition.sum(axis=1), axis=0)

    _ = plt.figure(1)
    plt.clf()
    axis = plt.subplot(111)
    axis.set_xlim(1e-1, 1e3)
    axis.set_ylim(1e-6, 2e0)
    axis.set_ylabel('Percentage counts')
    axis.set_xlabel('Sorted sample #')
    xvar = np.linspace(1, dim[0], dim[0])

    print composition.columns
    for column in composition.columns:
        color = 'white'
        if column in ['ResFinder', 'Bacteria', 'Bacteria_draft', 'HumanMicrobiome', 'Plasmid']:
            color = 'faded green'
        if column in ['Unmapped']:
            color = 'black'
        if column in ['Vertebrates_mammals', 'Vertebrates_other', 'Human', 'Mitochondrion']:
            color = 'denim blue'
        if column in ['Invertebrates', 'Plant', 'Virus']:
            color = 'amber'
        if column in ['MetalResistance']:
            color = 'pale red'
        axis.loglog(xvar, composition[column].sort_values()[::-1],
                    c=sns.xkcd_rgb[color], lw=3, alpha=0.8)
        axis.annotate(column, (9e-1, composition[column].sort_values()[-1]),
                      ha='right', fontsize=7)
