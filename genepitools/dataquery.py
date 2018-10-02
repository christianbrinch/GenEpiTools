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
from uncertainties.umath import *


class DataQuery(object):
    ''' A class that a selection of data to be plottet '''

    def __init__(self, data, metadata, first_day=np.datetime64('1970-01-01')):
        self.data = data
        self.metadata = metadata
        self.first_day = first_day
        self.length = len(self.days())

    def subset(self, column, values, excl=False):
        ''' Extract a subset of metadata based on input '''
        if not isinstance(values, list):
            values = values.split()
        if excl:
            return DataQuery(self.data, self.metadata[~self.metadata[column].isin(values)],
                             self.first_day)
        return DataQuery(self.data, self.metadata[self.metadata[column].isin(values)],
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

    def sum(self):
        ''' Sum of the data columns '''
        subset = self.data[self.metadata['Sample_ID']]
        return [sum(subset[column]) for column in self.metadata['Sample_ID']]

    def gene(self, gene):
        ''' Return a single gene abundance '''
        subset = self.data[self.metadata['Sample_ID']]
        return subset.loc[gene]

    def days(self):
        ''' Return a list days passed since the earliest date in data set '''
        return list((self.metadata['Sampling_date']-self.first_day).dt.days+1)