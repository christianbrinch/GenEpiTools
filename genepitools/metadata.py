# -*- coding: utf-8 -*-
""" Tools for handling metadata

"""

__author__ = "Christian Brinch"
__copyright__ = "Copyright 2018"
__credits__ = ["Christian Brinch"]
__license__ = "AFL 3.0"
__version__ = "0.1"
__maintainer__ = "Christian Brinch"
__email__ = "cbri@food.dtu.dk"


class SamplingSite(object):
    ''' A class to contain all metadata for a sampling site
    '''

    def __init__(self, attr):
        for key in attr:
            self.__dict__[key] = attr[key]
