# -*- coding: utf-8 -*-
''' A container for metagenomic count data and GenEpi related extensions to
    pandas.
'''

__author__ = "Christian Brinch"
__copyright__ = "Copyright 2018-2019"
__credits__ = ["Christian Brinch"]
__license__ = "AFL 3.0"
__version__ = "0.2"
__maintainer__ = "Christian Brinch"
__email__ = "cbri@gfood.dtu.dk"

import pandas as pd
import numpy as np
import scipy.stats as ss


@pd.api.extensions.register_dataframe_accessor("gen_count")
class GenAccessor(object):
    ''' An extension to pandas objects countaining counts '''

    def __init__(self, pandas_obj):
        # self._validate(pandas_obj)
        self._obj = pandas_obj

    # @staticmethod
    # def _validate(obj):
    #    if 'notPhiX' not in obj.index.values:
    #        raise AttributeError("No notPhiX row in mappings frame!")

    def filter_low_abundance(self, threshold):
        ''' Remove low abundant genes '''
        return self._obj[self._obj.T.sum() > threshold]

    def correct_for_gene_length(self, refdata):
        ''' A method that corrects the counts for the gene length '''
        return self._obj.iloc[:, :].div(refdata['refLength']/1000., axis=0)

    def homology_reduce(self, refdata):
        ''' Simple homology reduction based on the gene identifier '''
        return self._obj.join(refdata['RepresentativeSeqName']).groupby(
            'RepresentativeSeqName').sum()

    def agg(self, refdata):
        ''' Aggregate genes by resistance to AM classes and sort by abundance'''
        tmp_df = self._obj.join(refdata['description']).groupby('description').sum()
        tmp_df = tmp_df.assign(sum=tmp_df.sum(axis=1)).sort_values(
            by='sum', ascending=False).iloc[:, :-1]
        return tmp_df


@pd.api.extensions.register_dataframe_accessor("gen_map")
class MapAccessor(object):
    ''' An extension to pandas objects countaining mappings '''

    def __init__(self, pandas_obj):
        # self._validate(pandas_obj)
        self._obj = pandas_obj

    def agg(self):
        ''' A method that aggregates all types of Bacteria
            NOTE: Plasmids are currently not included'''
        return (self._obj.rename({'Bacteria': 'Bacteria_agg',
                                  'Bacteria_draft': 'Bacteria_agg',
                                  'HumanMicrobiome': 'Bacteria_agg',
                                  'MetaHitAssembly': 'Bacteria_agg'})
                .reset_index().groupby('Database').sum().drop('notPhiX'))


@pd.api.extensions.register_dataframe_accessor("gen_trans")
class TransAccessor(object):
    ''' An extension to pandas objects countaining various transformations '''

    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    def amal(self, kind):
        ''' Amalgamte all rows but one. Kind determines which one '''
        for idx, _ in self._obj.iterrows():
            if idx != kind:
                self._obj.rename({idx: 'Other'}, inplace=True)
        return self._obj.reset_index().groupby(self._obj.index).sum()

    def alr(self, denom=None, bayesian=False, n_samples=2000):
        ''' A method to calculate the additive log transform '''
        if denom is None:
            raise AttributeError("Please specify valid denominator")
        if not self._obj.values.all() and not bayesian:
            raise AttributeError(
                "Dataframe contains zeros. Use Bayesian inference to estimate zeros.")

        row = self._obj.index.get_loc(denom)

        return self._transform('alr', bayesian, n_samples, row)

    def clr(self, bayesian=False, n_samples=2000):
        ''' A method to calculate the centered log transform '''
        if not self._obj.values.all() and not bayesian:
            raise AttributeError(
                "Dataframe contains zeros. Use Bayesian inference to estimate zeros.")

        return self._transform('clr', bayesian, n_samples)

    def _transform(self, kind, bayesian, n_samples, row=None):
        ''' Calculate the alr/clr transforms '''
        error = self._obj.copy(deep=True)
        for column in self._obj:
            if bayesian is True:
                p_matrix = ss.dirichlet.rvs(self._obj[column]+0.5, n_samples)
            else:
                p_matrix = [self._obj[column]]

            if kind == 'clr':
                c_matrix = [np.log10(i) - np.mean(np.log10(i)) for i in p_matrix]
            else:
                c_matrix = [np.log10(i/i[row]) for i in p_matrix]
            self._obj[column] = [np.mean(i) for i in zip(*c_matrix)]
            error[column] = [np.std(i) for i in zip(*c_matrix)]

        return (self._obj, error)

    def geom_mean(self, bayesian=False, n_samples=5000):
        ''' Calculate the geometric mean of a composition '''

        if not self._obj.values.all() and not bayesian:
            raise AttributeError(
                "Dataframe contains zeros. Use Bayesian inference to estimate zeros.")

        for column in self._obj:
            if bayesian is True:
                p_matrix = ss.dirichlet.rvs(self._obj[column]+0.5, n_samples)
            else:
                p_matrix = [self._obj[column]]

            self._obj[column] = np.mean(p_matrix, axis=0)

        return [np.exp(np.mean(np.log(x))) for x in np.array(self._obj)]

    def totvar(self, mean):
        ''' Calculate the total variation of a composition
            TODO: This is very slow at the moment. Try to optimize
        '''
        variance = 0.
        dim = np.shape(self._obj)
        for vector in np.array(self._obj.T):
            dist_a = np.sqrt(1./(2.*dim[0])
                             * np.sum([(np.log(vector[i]/vector[j])
                                        - np.log(mean[i]/mean[j]))**2
                                       for i in range(dim[0])
                                       for j in range(dim[0])]))
            variance += 1./dim[1] * dist_a**2
        return variance


class DataSheet(object):
    ''' A class that holds all dataframes and transformation wrapper methods '''

    def __init__(self, counts, refdata, mappings, metadata):
        self.raw_counts = counts[(counts.T != 0).any()]
        self.raw_mappings = mappings
        self.metadata = metadata.T
        self.refdata = refdata
        self.errors = None

    def mappings(self, agg=True):
        ''' Wrapper for the mappings  '''
        if agg:
            return self.raw_mappings.gen_map.agg()
        return self.raw_mappings

    def counts(self, corr=False, reduc=False, filter_rare=0, agg=False):
        ''' Return count matrix with corrections '''
        tmp_df = self.raw_counts
        if filter_rare > 0:
            tmp_df = tmp_df.gen_count.filter_low_abundance(filter_rare)
        if corr:
            tmp_df = tmp_df.gen_count.correct_for_gene_length(self.refdata)
        if reduc:
            tmp_df = tmp_df.gen_count.homology_reduce(self.refdata)
        if agg:
            tmp_df = tmp_df.gen_count.agg(self.refdata)

        return tmp_df

    def days(self, subset=None):
        ''' Return a list days passed since the earliest date in data set '''
        if subset is None:
            return list((self.metadata.loc['Sampling_date'] -
                         min(self.metadata.loc['Sampling_date'])).dt.days + 1)
        return list((self.metadata.loc['Sampling_date'].loc[subset]
                     - min(self.metadata.loc['Sampling_date'])).dt.days + 1)


class Filter(list):
    ''' A list class extension containing a subset of sample metadata '''

    def __init__(self, datasheet):
        # list.__init__(self, datasheet.metadata.index.values.tolist())
        list.__init__(self, datasheet.metadata.columns.tolist())
        self.data = datasheet

    def days(self):
        ''' A wrapper method for outputting days '''
        return self.data.days(self)

    def isin(self, field, values):
        ''' Select samples where value is in field '''
        if not isinstance(values, list):
            values = values.split()
        tmp_frame = self.data.metadata.loc[:, self.data.metadata.columns.isin(self)]
        tmp_frame = tmp_frame.loc[:, tmp_frame.loc[field].isin(values)]
        super(Filter, self).__init__(tmp_frame.columns.tolist())
        return self

    def isnotin(self, field, values):
        ''' Select samples where value is not in field '''
        if not isinstance(values, list):
            values = values.split()
        tmp_frame = self.data.metadata.loc[:, self.data.metadata.columns.isin(self)]
        tmp_frame = tmp_frame.loc[:, ~tmp_frame.loc[field].isin(values)]
        super(Filter, self).__init__(tmp_frame.columns.tolist())
        return self
