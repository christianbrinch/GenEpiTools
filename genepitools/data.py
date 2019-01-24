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


class DataSheet(object):
    ''' A class that holds all dataframes and transformation methods '''

    def __init__(self, counts, refdata, mappings, metadata):
        self.counts = counts
        self.glc_counts = self.correct_for_gene_length(counts, refdata)
        self.clusters = self.homology_reduce(counts, refdata)
        self.am_classes = self.aggregate_by_class(counts, refdata)
        self.mappings = mappings
        self.metadata = metadata
        self.firstday = min(metadata['Sampling_date'])

    @classmethod
    def correct_for_gene_length(cls, counts, refdata):
        ''' A method that corrects the counts for the gene length '''

        tmp_df = counts.iloc[:, :].div(refdata['refLength']/1000., axis=0)

        return tmp_df

    @classmethod
    def aggregate_by_class(cls, counts, refdata):
        ''' Aggregate genes by resistance to AM classes '''
        tmp_df = cls.correct_for_gene_length(counts, refdata)
        tmp_df = tmp_df.join(refdata['description']).groupby('description').sum()
        tmp_df = tmp_df.assign(sum=tmp_df.sum(axis=1)).sort_values(
            by='sum', ascending=False).iloc[:, :-1]
        tmp_df = tmp_df[(tmp_df.T != 0).any()]
        return tmp_df

    @classmethod
    def homology_reduce(cls, counts, refdata):
        ''' Simple homology reduction based on the gene identifier '''
        tmp_df = cls.correct_for_gene_length(counts, refdata)
        tmp_df = tmp_df.join(refdata['RepresentativeSeqName']
                             ).groupby('RepresentativeSeqName').sum()
        tmp_df = tmp_df[(tmp_df.T != 0).any()]
        return tmp_df

    def days(self, subset=None):
        ''' Return a list days passed since the earliest date in data set '''
        if subset is None:
            return list((self.metadata['Sampling_date']-self.firstday).dt.days
                        + 1)
        return list((self.metadata['Sampling_date'].loc[subset]
                     - self.firstday).dt.days + 1)


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
