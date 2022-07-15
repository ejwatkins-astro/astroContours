# -*- coding: utf-8 -*-
"""
Script contains classes and functions that wrap around methods that find
contours (a list of boundaries in rows and columns) from a labled-maksed image
"""

import matplotlib.pyplot as plt
import numpy as np
from skimage import measure
from .. import fitsTools as mf #wrappers for astropy.io
from .. import maskContoursUtility as m2c

from astropy.io import fits
import functools
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

# No to tickers sticking out!
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
# Makes mathematical text regular rather than italics
plt.rcParams['mathtext.default'] = 'regular'

#Decorator than removes extra dimentsions
def correct_dimension(func):
    def wrapper(self, *args, **kwargs):
        output = func(self, *args, **kwargs)
        output_shape = np.shape(output)
        if len(output_shape) == 2 and output_shape[0] == 1:
            return output[0]
        elif len(output_shape) == 3 and output_shape[1] == 1:
            return [out[0] for out in output]
        else:
            return output
    return wrapper

class Contours:

    """
    This class is a wrapper around the functions needed to convert either a
    2d or 3d masked astronomical data into their coresponding contours.
    It accepts labeled masked data -- where 0 is no mask, and 1, 2, 3 etc --
    are the ID's of the masked data. If you want the on sky projection,
    the 'header' needs to be given too.
    Using the method `from_fits_file` allows you to just enter the path
    to the masked fits image instead

    Most contour finding algorthms do not work well if the masked touch
    (i.e., masks for object 1 and 2 have no 0's between them). I have
    written a work around but it is much slower. So if you know that all
    your masks are separate, change `touching_masks` to be False.

    'contour_dilation' sets how the boundary is cast (it is the contour level).

    But in most contour finding algoithms, the level=1 dilates the contour
    half a pixel larger than the input mask. In gerneral I found:

    0 contracts the contour boundary by half a pixel.
    0.5 approximates the exact boundary
    1 dilates the contour boundary by half a pixel.
    So leave this setting as `0.5`unless you want the behaviour

    """

    def __init__(self, labeled_mask_data, header=None, contour_dilation=0.5, touching_masks=True, add_to_mask=0):

        data_shape = labeled_mask_data.shape
        if len(data_shape) == 2:
            self.labeled_mask_data = np.array([labeled_mask_data])
        else:
            self.labeled_mask_data = labeled_mask_data + add_to_mask

        if header is not None:
            header = mf.remove_axis_from_fits_header(header, 4)
            header = mf.remove_axis_from_fits_header(header, 3)
        self.header = header

        self.contour_dilation = contour_dilation
        self.touching_masks = touching_masks

    @classmethod
    def from_fits_file(cls, filename, hdu_i=0, contour_dilation=0.5, touching_masks=True, add_to_mask=0):
        with fits.open(filename) as hdulist:
            data = hdulist[hdu_i].data + add_to_mask
            header = hdulist[hdu_i].header
        return cls(data, header, contour_dilation, touching_masks)

    @functools.lru_cache # caching result to speed up any potential following uses of this
    def _get_contours_yx(self):
        """
        Given a masked, label image, this method will find the contours
        in pixel coordinates for each mask

        """
        contours_yx, self.contour_ids = map(list, zip(*[m2c.get_contours(labeled_image=mask,
                                                      contour_dilation=self.contour_dilation,
                                                      get_contour_id=True,
                                                      touching_masks=self.touching_masks) for mask in self.labeled_mask_data]))

        return contours_yx, self.contour_ids

    def set_header(self, header):
        self.header = header

    def _get_contours_wcs(self):

        contours_yx, contour_ids = self._get_contours_yx()
        contours_wcs = [get_contours_wcs_from_contours_yx(contours, self.header) for contours in contours_yx]

        return contours_wcs, contour_ids

    @correct_dimension
    def get_contours_yx(self):
        """
        Given a masked, label image, this method will find the contours
        in pixel coordinates for each mask

        """
        contours_yx, contour_ids = self._get_contours_yx()

        return contours_yx,  contour_ids


    @correct_dimension
    def get_contours_wcs(self):
        """
        Given a masked image, this method will find the contours
        in wcs coordinates for each mask
        """
        contours_wcs, contour_ids = self._get_contours_wcs()

        return contours_wcs, contour_ids

    @correct_dimension
    def convert_contours_yx_between_frames(self, new_header):
        """
        Convenience function. This method will output the
        contours in the pixel coordinate of the header given.

        """
        new_header = mf.remove_axis_from_fits_header(new_header, 4)
        new_header = mf.remove_axis_from_fits_header(new_header, 3)
        contours_wcs, contour_ids = self._get_contours_wcs()
        contours_yx_new = [get_contours_yx_from_contours_wcs(contours, new_header) for contours in contours_wcs]

        return contours_yx_new, contour_ids

    @correct_dimension
    def convert_contours_wcs_to_new_world(self, new_world_type):
        contours_wcs, contour_ids = self._get_contours_wcs()
        current_world_type = get_world_type(self.header)

        contours_wcs_new = [get_new_world_contours_wcs(contours, current_world_type, new_world_type) for contours in contours_wcs]

        return contours_wcs_new, contour_ids


class Levels_to_contours(Contours):

    def __init__(self, data, levels, header=None):

    # data_shape = data.shape
    # if len(data_shape) == 2:
    #     self.data = data[None,:,:]#np.array([data])
    # else:
        self.data = data
        self.levels = levels

        if header is not None:
            header = mf.remove_axis_from_fits_header(header, 4)
            header = mf.remove_axis_from_fits_header(header, 3)
        self.header = header

    @classmethod
    def from_fits_file(cls, filename, levels, hdu_i=0):
        with fits.open(filename) as hdulist:
            data = hdulist[hdu_i].data
            header = hdulist[hdu_i].header
        return cls(data=data, levels=levels, header=header)

    @functools.lru_cache # caching result to speed up any potential following uses of this
    def _get_contours_yx(self):
        """
        Given a masked, label image, this method will find the contours
        in pixel coordinates for each mask

        """
        contours_yx = [measure.find_contours(self.data, level=lvl) for lvl in self.levels]
        hole_locations = self.find_holes_using_orientation(contours_yx)

        return contours_yx, hole_locations

    def find_holes_using_orientation(self, contours_yx):
        hole_locations = [find_all_holes(contours) for contours in contours_yx]
        return hole_locations

world_types = {
    'RA':'fk5',
    'DEC':'fk5',
    'GLON':'galactic',
    'GLAT':'galactic'
    }
#Maybe move this function to untilities?
def get_world_type(header):
    try:
        current_world = header['RADESYS'].lower()
    except KeyError: # will guess from the CTYPE using `world_types` dictionary
        ctype = header['CTYPE1']
        for key in world_types.keys():
            if key in ctype:
                current_world = world_types[key]
                break
            else:
                raise KeyError('World system could not be identified')
    return current_world

def find_all_holes(contours_yx):
    """
    1 is a contour, 0 is a hole

    Parameters
    ----------
    contours_yx : 2d array like
        list of y,x vertices of a contour

    Returns
    -------
    hole_locations: 1d array like

    """

    hole_locations = [m2c.find_hole_contours(cont) for cont in contours_yx]

    return hole_locations

#These are wrappers around the functions astropy wcs funtions
def get_new_world_contours_wcs(contours_wcs, current_wcs, new_wcs):
    """
    This method converts contours to different world types
    """
    #contours are found with vertices y,x. wcs requires x,y. [:,::-1]
    #switches from y,x to x,y
    new_contours_wcs = [m2c.convert_between_worlds(cont, current_wcs, new_wcs) for cont in contours_wcs]

    return new_contours_wcs


def get_contours_wcs_from_contours_yx(contours_yx, header):
    """
    This method converts contours and to WCS coordinates

    """
    #contours are found with vertices y,x. wcs requires x,y. [:,::-1]
    #switches from y,x to x,y
    contours_wcs = [m2c.convert_pixel2world(cont[:,::-1], header) for cont in contours_yx]

    return contours_wcs

def get_contours_yx_from_contours_wcs(contours_wcs, new_header):
    """
    Given contours in wcs coordinates this method converts the contours to
    the given projection in the header.

    """
    contours_yx = [m2c.convert_world2pixel(cont, new_header)[:,::-1] for cont in contours_wcs]

    return contours_yx

#quick plotting function if you want to plot the contours
def plot_contours(contour_yx_list, ax=None, reverse_yx=False, **kwargs):
    """
    Given a list of contours, this function will plot them. Kwargs are for
    `plt.plot`
    """
    label = kwargs.pop('label', None)
    if ax is None:
        ax = plt.gca()

    for cont in contour_yx_list[1:]:
        x, y = cont[:,1], cont[:,0]
        if reverse_yx:
            x, y  = y, x

        ax.plot(x, y, **kwargs)

    x, y = cont[0,1], cont[0,0]
    if reverse_yx:
        x, y  = y, x

    ax.plot(x, y, label=label, **kwargs)

def get_single_contour_from_3d_mask(mask3d, ids):
    """
    This collapses a 3d mask catalogue into a single 2d projected contours
    (i.e., the CPROPS GMC catalogue)

    """
    single_mask = np.zeros(mask3d.shape[1:])

    remove_ind = []
    for i in range(len(mask3d)):
        if np.all(mask3d[i] == 0):
            remove_ind.append(i)
    mask3d = np.delete(mask3d, remove_ind, axis = 0)

    contours = []
    contour_ids = []
    for i in range(1, len(ids)+1):
        print(100*i/len(ids))
        single_mask = np.zeros(mask3d.shape[1:])
        z, y, x = np.where(mask3d == i)
        single_mask[y, x] = i
        c = Contours(labeled_mask_data=single_mask, header=None, touching_masks=False)
        cont, cont_id = c.get_contours_yx()
        contours.append(cont[0])
        contour_ids.append(cont_id[0])

    return contours, contour_ids, c

# def convert_point_sources_yx_to_contour(point_sources_yx):
#     """
#     This is a vectorised operation for converting an entire
#     catalogue of point sources into small contours around thier
#     pixel boundary.
#     Note: If you have a point source catalogue, can just use
#     """
#     point_source_yx = np.array(point_sources_yx)

#     x_values = np.array([-0.5, -0.5, 0.5, 0.5, -0.5])
#     y_values = np.array([-0.5, 0.5, 0.5, -0.5, -0.5])

#     yx_values = np.column_stack((y_values, x_values))

#     #the [:,None] adds the extra dimension needed for the vectorised operation
#     pixel_contours = point_sources_yx[:,None] + yx_values

#     return pixel_contours
