# -*- coding: utf-8 -*-
"""
Script for loading contours (a list of boundaries in rows and columns) from
different formats
"""
import numpy as np
import json
import regions
import pickle


def main():
    from . import contoursDS9 as ds
   # a1()


# from . import contoursDS9 as ds

class load_contours():
    """
    Load contours saved in different formats
    """
    def __init__(self, contours, contour_ids, filename, **kwargs):

        self.contours = contours
        self.contour_ids = contour_ids
        self.filename = filename
        self.filename_ids = filename + '_ids.txt'
        self.kwargs = kwargs

    def set_filename_for_contour_ids(self, filename):
        """
        Function to change the filename used to load the contour id's
        """
        self.filename_ids = filename

    def load_contour_ids_txt(self):
        """
        Loads the contour id's from a txt file
        """
        data_ids  = np.loadtxt(self.filename_ids, self.contour_ids, dtype=int, **self.kwargs)
        return data_ids

    def load_contours_npy(self):
        """
        Load contours from a numpy object
        """
        data = np.load(self.filename+'.npy', allow_pickle=True, **self.kwargs)
        return data

    def load_contours_pickle(self):
        """
        Load contours from pickled object
        """
        with open(self.filename + '.pk', 'rb') as get:
            data = pickle.load(get, **self.kwargs)
        return data

    def load_contours_json(self):
        """
        Load contours from json file
        """
        with open(self.filename + '.json') as json_file:
            data = json.load(json_file, **self.kwargs)
        return data

    def load_contours_wcs_ds9(self):
        """
        Loads wcs contours from ds9 region file. Recomend using the class object
        """
        ds9 = ds.ds9_regions_contours(self.filename + '.reg')
        data = ds9.get_contours_wcs_fk5(**self.kwargs)
        return data

if __name__ == "__main__":
    main()