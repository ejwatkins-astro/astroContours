# -*- coding: utf-8 -*-
"""
Script for saving contours (a list of boundaries in rows and columns)
"""
import numpy as np
import json
import regions
import pickle

from . import contoursDS9 as ds9

class save_contours():
    """
    Class for saving contours in different formats
    """
    def __init__(self, filename, contours, contour_ids, save_ids_as_txt=True, save_as_dictionary=False):

        if save_as_dictionary:
            self.contours = make_dictionary_from_contours(contours, contour_ids)
        else:
            self.contours = contours

        self.contour_ids = contour_ids
        if '.' in filename:
            print('No `.` files needed. The right extension will be added with the save')
        self.filename = filename
        self.filename_ids = filename + '_ids.txt'


    def set_filename_for_contour_ids(self, filename):
        """
        sets the filename used to save contour ids
        """
        self.filename_ids = filename

    def save_contour_ids_txt(self):
        """
        Saves the contour id's as a .txt file
        """
        np.savetxt(self.filename_ids, self.contour_ids, fmt='%d')

    def save_contours_npy(self):
        """
        numpy object is simualar to a pickle. Can save complex python
        instances. Saves as a .npy file
        """
        np.save(self.filename , self.contours)
        self.save_contour_ids_txt()


    def save_contours_pickle(self):
        """
        Save as a pickled object
        """
        with open(self.filename+'.pk', 'wb') as outfile:
            pickle.dump(self.contours, outfile, pickle.HIGHEST_PROTOCOL)
        self.save_contour_ids_txt()

    def save_contours_json(self):
        """
        Saves a dictionary like object as a json file
        """
        with open(self.filename+'.json', 'w') as outfile:
            json.dump(self.contour, outfile)
        self.save_contour_ids_txt()

    def save_contours_ds9(self, wcs_projection='fk5'):
        """
        Converts then saves contours as ds9 regions. This process is slow
        """
        ds9_regions = ds9.contours_to_ds9_regions(self.contours, self.contour_ids, wcs_projection)
        regions.write_ds9(ds9_regions.get_ds9_regions(), self.filename+'.reg')

        self.save_contour_ids_txt()


def make_dictionary_from_contours(contours, contour_ids):
    """
    Note: This will mess up the order of where holes are located. Typically
    holes will be listed directly after the contour it is within. Making
    a dictionary will shove all holes into a single key.

    This format means disconected contours will be saved under one ID key

    """
    contour_dictionary = {}
    for i in range(len(contour_ids)):
        id_key = str(int(contour_ids[i]))
        try:
            obj_props[id_key].append(contours[i])
        except KeyError:
            obj_props[id_key] = [contours[i]]
