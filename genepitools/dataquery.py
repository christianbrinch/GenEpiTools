# -*- coding: utf-8 -*-
''' Extract subset of data

'''

__author__ = "Christian Brinch"
__copyright__ = "Copyright 2018"
__credits__ = ["Christian Brinch"]
__license__ = "AFL 3.0"
__version__ = "0.1"
__maintainer__ = "Christian Brinch"
__email__ = "cbri@gfood.dtu.dk"

import numpy as np
import pandas as pd
import random as rd
import scipy.stats as ss


class DataQuery(object):
    ''' A class that a selection of data to be plottet '''

    def __init__(self, data, error, metadata, first_day=np.datetime64('1970-01-01')):
        self.data = data
        self.error = error
        self.metadata = metadata
        self.first_day = first_day
        self.length = len(self.days())

    def subset(self, column, values, excl=False):
        ''' Extract a subset of metadata based on input '''
        if not isinstance(values, list):
            values = values.split()
        if excl:
            return DataQuery(self.data, self.error,
                             self.metadata[~self.metadata[column].isin(values)],
                             self.first_day)
        return DataQuery(self.data, self.error,
                         self.metadata[self.metadata[column].isin(values)],
                         self.first_day)

    def flag(self, gene, threshold, excl=True):
        ''' Filter data based on a threshold value '''
        if excl:
            for column in self.data.loc[gene].index:
                if self.data[column].loc[gene] > threshold:
                    self.data[column].loc[gene] = 0.
        else:
            for column in self.data.loc[gene].index:
                if self.data[column].loc[gene] < threshold:
                    self.data[column].loc[gene] = 0.

        return DataQuery(self.data, self.metadata, self.first_day)

    def gene(self, gene):
        ''' Return a single gene abundance '''
        subset = self.data[self.metadata['Sample_ID']]
        return subset.loc[gene]

    def gene_err(self, gene):
        ''' Return a single gene abundance error'''
        subset = self.error[self.metadata['Sample_ID']]
        return subset.loc[gene]

    def shannon(self, n_samples):
        ''' Return the Shannon index. Downsample to 800 counts.
            Todo: This can (and should) be done better. Also, sample N times
            and average to get better fidelity.
        '''
        subset = self.data[self.metadata['Sample_ID']]
        shannon = pd.DataFrame(index=[], columns=['H', 'L'])

        for column in subset:

            p_matrix = ss.dirichlet.rvs(np.array(subset[column])+0.5, n_samples)
            a_vector = np.mean(p_matrix, axis=0)
            sample = np.random.choice(a_vector, size=1000)
            s = -np.sum([np.log10(i)*i for i in sample])

            shannon.loc[column] = (s, len(
                np.array([i for i in list(subset[column]) if i > 0.])))
            '''
            reduc = np.array([i for i in list(subset[column]) if i > 0.])
            reduc = reduc/np.sum(reduc)
            tmp = np.zeros(len(reduc))
            count = 0
            while count < n_samples:
                myint = rd.randint(0, len(reduc)-1)
                if np.random.rand() < reduc[myint]:
                    tmp[myint] += 1
                    count += 1
            tmp = [i for i in tmp if i > 0.]
            tmp = tmp/np.sum(tmp)
            '''

            #shannon.loc[column] = (- np.sum(values*np.log(values)), 60)

        return shannon['H'], shannon['L']

    def days(self):
        ''' Return a list days passed since the earliest date in data set '''
        return list((self.metadata['Sampling_date']-self.first_day).dt.days+1)
