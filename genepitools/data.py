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
import scipy.optimize as so
import genepitools.transforms as gtrans
import genepitools.functions as gfunc


class DataSheet(object):
    ''' A class that holds all dataframes and transformation methods '''

    def __init__(self, counts, mappings, metadata):
        self.glc_counts = self.correct_for_gene_length(counts)
        self.mappings = mappings
        #self.mappings = self.correct_for_sequencing_depth(mappings)
        self.metadata = metadata
        self.am_classes = self.aggregate_by_class(counts)
        self.firstday = min(metadata['Sampling_date'])

    @classmethod
    def correct_for_sequencing_depth(cls, mappings):
        ''' A method that corrects for the nonlinear response to variations
            in the sequencing depth.
            TODO: This should be applied to the count frames as well
        '''
        dbase = ['ResFinder', 'Plasmid', 'HumanMicrobiome', 'Bacteria', 'Bacteria_draft']
        for i in range(len(dbase)):
            dat_x = np.log10(mappings.sum(axis=0))
            dat_y = np.log10(mappings.loc[dbase[i]])
            popt, _ = so.curve_fit(gfunc.line, dat_x, dat_y)
            corr = (dat_y - (gfunc.line(dat_x, *popt)-(dat_x-np.mean(dat_x)+np.mean(dat_y))))
            mappings.loc[dbase[i]] = pow(10, corr)

        return mappings

    @classmethod
    def correct_for_gene_length(cls, counts):
        ''' A method that corrects the counts for the gene length '''
        tmp_df = counts.iloc[:, 1:].div(counts['GeneSize']/1000., axis=0)
        return tmp_df.T.drop(['GeneSize']).T

    @classmethod
    def aggregate_by_class(cls, counts):
        ''' Aggregate genes by resistance to AM classes '''
        tmp_df = counts.iloc[:, 1:].div(counts['GeneSize']/1000., axis=0)
        tmp_df = tmp_df.join(counts['Description']).groupby('Description').sum()
        return tmp_df.T.drop(['GeneSize']).T

    def days(self, subset=None):
        ''' Return a list days passed since the earliest date in data set '''
        if subset is None:
            return list((self.metadata['Sampling_date']-self.firstday).dt.days+1)
        return list((self.metadata['Sampling_date'].loc[subset]-self.firstday).dt.days+1)


class Filter(list):
    ''' A list subclass containing a subset of sample metadata '''

    def __init__(self, datasheet):
        list.__init__(self, datasheet.metadata.index.values.tolist())
        self.data = datasheet

    def days(self):
        ''' A wrapper method for outputting days '''
        return self.data.days(self)

    def isin(self, field, values):
        ''' Select samples where value is in field '''
        if not isinstance(values, list):
            values = values.split()
        tmp_frame = self.data.metadata[self.data.metadata.index.isin(self)]
        tmp_frame = tmp_frame[tmp_frame[field].isin(values)]
        super(Filter, self).__init__(tmp_frame.index.values.tolist())
        return self

    def isnotin(self, field, values):
        ''' Select samples where value is not in field '''
        if not isinstance(values, list):
            values = values.split()
        tmp_frame = self.data.metadata[self.data.metadata.index.isin(self)]
        tmp_frame = tmp_frame[~tmp_frame[field].isin(values)]
        super(Filter, self).__init__(tmp_frame.index.values.tolist())
        return self
