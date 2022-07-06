# -*- coding: utf-8 -*-
"""
Script for loading contours (a list of boundaries in rows and columns) from
different formats. NOT FINISHED. ONLY WORKS FOR FITS FILES. NEED TO ADD CODE
THAT LOADS THE ID'S FOR THE CONTOURS OR CONVERTS DICTIONARY-LIKE FORMATS INTO
CONTOUR_ID'S FOR THE KEYS AND CONTOURS FOR THE VALUES.
"""
import numpy as np
import json
import regions
import pickle


def main():
    from . import contoursDS9 as ds
   # a1()


# from . import contoursDS9 as ds

"""
Script for loading in data saved in a few different formats
"""

import json
import numpy as np
import pickle
from astropy.io import fits

# files that can be loaded in io
file_types = {
    'fits':'from_fits_file',
    'pkl' :'from_pickle',
    'npy' :'from_npy',
    'json':'from_json',
    'reg' :'from_ds9_regions'
    }

class Loader:
    """
    Class for loading in (astro) data saved in different formats given in
    the dictonary `file_types`. The lazy loader will always output
    the data and a `header`. If no `header` was given or found, the `header`
    will be returned as a `None`

    """
    def __init__(self, filename, filename_ids=None, **file_args):

        self.filename = filename
        self.filename_ids = filename.split('.')[0] + '.txt' if filename_ids is None else filename_ids
        self.header = file_args.pop('header', None)
        self.file_args = file_args
        self.hdu_i = self.file_args.pop('hdu_i', 0)
        file_type = self.filename.split('.')[-1]

        try:
            self.method = file_types[file_type]
        except KeyError as error:

            print('%s is not an accepted format. Loaded file must be one of the following:' % error, file_types.keys())
            raise

    def load_from_fits_file(self, **kwargs):
        """
        Load data from a fits file. If loading a fits file
        and you need to specify the extention, add the extension to
        the variable `hdu_i` when initalising the class

        Parameters
        ----------
        **kwargs : dict
            Optional parameter

        Returns
        -------
        data : ndarray
            The data
        header : astropy.io.Header.header
            Header associated with the data. Acts like a dictionary

        """
        with fits.open(self.filename) as hdulist:
            data = hdulist[self.hdu_i].data
            self.header = hdulist[self.hdu_i].header

        return data, self.header

    def load_from_pickle(self, **kwargs):
        """
        Loads data that has been pickled

        Parameters
        ----------
        **kwargs : dict
            Optional parameters for pickle load.

        Returns
        -------
        data : ndaarray or object or dict
            Data saved within the pickle.

        """
        with open(self.filename, 'rb') as get:
            data = pickle.load(get, **kwargs)

        return data

    def load_from_json(self, **kwargs):
        """
        Load from json file

        Parameters
        ----------
        **kwargs : dict
            Optional parameters for json load.

        Returns
        -------
        data : ndaarray or object or dict
            Data saved within the json.

        """

        with open(self.filename + '.json') as get:
            data = json.load(get, **kwargs)

        return data

    def load_from_npy(self, **kwargs):
        """
        Load from a numpy object. Numpy has its own file type when using
        `numpy.save`. Can save data in this format and load it back in

        Parameters
        ----------
        **kwargs : dict
            Optional parameters for numpy load.

        Returns
        -------
        data : ndaarray or object or dict
            Data saved within the json.

        """
        allow_pickle = kwargs.pop('allow_pickle', True)
        data = np.load(self.filename, allow_pickle=allow_pickle, **kwargs)
        return data

    def load_from_ds9_regions(self,  **kwargs):
        """
        Loads wcs contours from ds9 region file. This is very slow
        """
        ds9 = ds.ds9_regions_contours(self.filename + '.reg')
        data = ds9.get_contours_wcs_fk5(**kwargs)
        return data

    def lazy_load(self):
        """
        Can load your data using this function. If loading a fits file
        and you need to specify the extention, add the extension to
        the variable `hdu_i` when initalising the . The default is 0

        Returns
        -------
        data : ndaarray or object or dict
            Data saved within the file.
        header : astropy.io.Header.header or None
            Header associated with the data. Acts like a dictionary. If
            `header` is not given when initalising and loaded data was
            not from fits file, `header`=None

        """

        data = getattr(self, 'load_' + self.method)(**self.file_args)

        if self.method == 'from_fits_file' :
            data = data[0]

        return data, self.header

# class load_contours():
#     """
#     Load contours saved in different formats
#     """
#     def __init__(self, contours, contour_ids, filename, **kwargs):

#         self.contours = contours
#         self.contour_ids = contour_ids
#         self.filename = filename
#         self.filename_ids = filename + '_ids.txt'
#         self.kwargs = kwargs

#     def set_filename_for_contour_ids(self, filename):
#         """
#         Function to change the filename used to load the contour id's
#         """
#         self.filename_ids = filename

#     def load_contour_ids_txt(self):
#         """
#         Loads the contour id's from a txt file
#         """
#         data_ids  = np.loadtxt(self.filename_ids, self.contour_ids, dtype=int, **self.kwargs)
#         return data_ids

#     def load_contours_npy(self):
#         """
#         Load contours from a numpy object
#         """
#         data = np.load(self.filename+'.npy', allow_pickle=True, **self.kwargs)
#         return data

#     def load_contours_pickle(self):
#         """
#         Load contours from pickled object
#         """
#         with open(self.filename + '.pk', 'rb') as get:
#             data = pickle.load(get, **self.kwargs)
#         return data

#     def load_contours_json(self):
#         """
#         Load contours from json file
#         """
#         with open(self.filename + '.json') as json_file:
#             data = json.load(json_file, **self.kwargs)
#         return data



if __name__ == "__main__":
    main()