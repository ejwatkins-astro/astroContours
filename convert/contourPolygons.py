# -*- coding: utf-8 -*-
"""
Script contains classes and functions that wrap around methods that convert
contours (a list of boundaries in rows and columns) into polygon object that
can be used in shapely and geopanda data frames.
"""


import shapely.geometry as sg
import descartes
from matplotlib.collections import PatchCollection
import numpy as np
import copy as copy
import functools

from .. import load
from . import maskContours as mC

import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)



def correct_dimension(func):
    def wrapper(self):
        output = func(self)
        output_shape = np.shape(output)
        if len(output_shape) == 2 and output_shape[0] == 1:
            return output[0]
        elif len(output_shape) == 3 and output_shape[1] == 1:
            return [out[0] for out in output]
        else:
            return output
    return wrapper

class ContoursAsPolygons():

    """
    Wrapper to converts a set of contours from a 2d or 3d image into
    polygon objects, while trying to correctly deal with 'holes' from
    the origonal object.

    Requires a list of contours for each 2d slice (so if the image is only
    2d, you will need to add a list around it i.e. this > [[cont1, cont2]]
    not this > [cont1, cont2]

    Requires the corresponding ID's, strucutres simulary as the contours

    Optional boolean for if the imputted coordinates of the contours are
    ordered x first then y. Default is False (since python normally
    sets things to y then x.)

    Note: Each contour is a 2d numpy array of separate coordinates i.e.,
    cont1 = np.array([[y1, y2, y3], [x1, x2, x3]]).

        """

    def __init__(self, contours, contour_ids, are_contours_xy=False, fix_invalid_geometries=False):

        if len(np.shape(contour_ids)) == 1:
            self.contour_ids = [contour_ids]
        if len(np.shape(contours)) == 1:
            self.contours = [contours]
        else:
            self.contours = contours


        self.num_of_3d = len(self.contours)
        self.fix_invalid_geometries = fix_invalid_geometries

        self.r = -1 # contours are usully yx, but polygon needs xy
        if are_contours_xy:
            self.r = 1

    @classmethod
    def from_file(cls, filename, convert_wcs=False,
                  current_header=None, new_header=None,
                  hdu_i=0, contour_dilation=0.5,
                  touching_masks=True, fix_invalid_geometries=False):

        data, current_header = load.Loader(filename, header=current_header, hdu_i=hdu_i).lazy_load()

        cont_obj = mC.Contours(data, header=current_header,
                               contour_dilation=contour_dilation,
                               touching_masks=touching_masks)

        return ContoursAsPolygons.from_contour_object(cont_obj, convert_wcs, current_header, new_header)

    @classmethod
    def from_contour_object(cls, contour_object, convert_wcs=False,
                            current_header=None, new_header=None,
                            fix_invalid_geometries=False):

        are_contours_xy = False
        if convert_wcs:
            if new_header is not None:
                contours, contour_ids = contour_object.convert_contours_yx_between_frames(new_header=new_header)
            else:
                contours, contour_ids = contour_object.get_contours_wcs()
                are_contours_xy = True
        else:
            contours, contour_ids = contour_object.get_contours_yx()

        return cls(contours, contour_ids, are_contours_xy=are_contours_xy, fix_invalid_geometries=fix_invalid_geometries)

    @functools.lru_cache
    def _get_holed_polygons(self):
        """
        Turns the contours into polygons, accounting for holes
        """
        polygons, self.polygons_ids = map(list, zip(*[get_holed_polygons_from_contours(self.contours[i], self.contour_ids[i], self.r, self.fix_invalid_geometries) for i in range(self.num_of_3d)]))
        return polygons, self.polygons_ids

    @correct_dimension
    def get_polygons(self):
        """
        Turns the contours into polygons, but each hole is given a polygon
        """
        return [get_polygons_from_contours(contours_yx, self.r, self.fix_invalid_geometries) for contours_yx in self.contours]

    def get_hole_info(self):
        hole_posistions, hole_posistion_indices, number_of_holes = map(list, zip(*[get_hole_info_from_contour_ids(contour_ids) for contour_ids in self.contour_ids]))
        return hole_posistions, hole_posistion_indices, number_of_holes

    @correct_dimension
    def get_holed_polygons(self):
        """
        Turns the contours into polygons, accounting for holes
        """
        polygons, polygons_ids = self._get_holed_polygons()
        return polygons, polygons_ids

    @correct_dimension
    def get_patches(self):
        """
        Turns the contours directly into patches, accounting for holes
        """
        polygon_contours, polygons_ids = self._get_holed_polygons()
        patches = [get_patches_from_shapely_objects(polygons) for polygons in polygon_contours]
        return patches, polygons_ids


def get_polygons_from_contours(contours_yx, r=-1, fix_invalid_geometries=False):
    """
    Wrapper function for converting a list of contours into shapely geometry
    polygons. If the contours are xy, not yx, set `r` to equal 1

    Parameters
    ----------
    contours_yx : list
        List containing a set of contours i.e, [cont1, cont2, cont3]
    r : int
        Integer that indicates and changes the  order of the coordinates of
        each contour. Python normally creates objects ordered y, x, but most
        polygon codes require x, y. if y, x, the order is reversed by setting
        `r=-1`. If the order is already x, y, set `r=1`. The default for this
        function is -1.
    fix_invalid_geometries : bool
        Rejects invalid geometries, such as self intercections (lines, bowties).
        `buffer(0)` will try and make the polygon valid in anyway it can, so
        might drastically change the geometry. If no valid geometry is found
        it returns an empty polygon.

    Returns
    -------
    polygon_contours_xy : list
        List containing polygon objects i.e., [poly1, poly2 poly3]

    """
    if fix_invalid_geometries:
        polygon_contours_xy = [sg.Polygon(cont_yx[:,::r]).buffer(0) for cont_yx in contours_yx]
    else:
        polygon_contours_xy = [sg.Polygon(cont_yx[:,::r]) for cont_yx in contours_yx]
    return polygon_contours_xy

def get_points_from_coordinates(point_sources_yx, r=-1, buffer=None):

    if buffer is not None:
        shapely_points_xy = [sg.Point(*point_source[::r]).buffer(buffer) for point_source in point_sources_yx]
    else:
        shapely_points_xy = [sg.Point(*point_source[::r]) for point_source in point_sources_yx]

    return shapely_points_xy

def get_hole_info_from_contour_ids(contour_ids):
    """
    Function identify if and what contours are not objects, but are holes
    TODO: skimage.measure defines holes by making the coordinates anti-clockwise
    and they are listed directly after the contour containing the hole. Use this
    to better locate holes and correctly remove only the holes

    Parameters
    ----------
    contour_id : List of ints
        List containing the ID's corresponding to each contour.

    Returns
    -------
    hole_posistions : Boolean list
        List containing the boolean True when a hole was found

    hole_posistion_indices : 1d numpy array of ints
        indicies where the object is a hole (i.e., =0).
    number_of_holes : int
        Number of holes identified in the set of contours

    """
    if isinstance(contour_ids, list):
        contour_ids = np.array(contour_ids)
    hole_posistions = contour_ids == 0
    hole_posistion_indices = np.where(contour_ids == 0)[0]
    number_of_holes = len(hole_posistion_indices)

    return hole_posistions, hole_posistion_indices, number_of_holes

def get_holed_polygons_from_contours(contours_yx, contour_ids, r=-1, fix_invalid_geometries=False):
    """
    Main function for converting contours to correct polygons
    accounting for if the mask/contour has a hole


    Parameters
    ----------
    contours_yx : List
        List containing a set of contours i.e, [cont1, cont2, cont3]
    contour_id : List
        List containing the ID's corresponding to each contour.
    r : int
        Integer that indicates and changes the  order of the coordinates of
        each contour. Python normally creates objects ordered y, x, but most
        polygon codes require x, y. if y, x, the order is reversed by setting
        `r=-1`. If the order is already x, y, set `r=1`. The default for this
        function is -1.
    fix_invalid_geometries : bool
        Rejects invalid geometries, such as self intercections (lines, bowties)

    Returns
    -------
    polygon_contours_xy_copy : list
        List containing polygon objects with the holes defined [poly1, poly2 poly3]
    contour_id : List
        List containing the ID's corresponding to each polygon ammended to
        remove the holes that are no longer induvidual holes but are defined as
        part of the polygon is present.

    """
    #hole information. i.e if the contours have holes and how many etc
    hole_posistions, hole_posistion_indices, number_of_holes = get_hole_info_from_contour_ids(contour_ids)

    #this is a loop that turns the contours into polygon objects
    polygon_contours_xy = get_polygons_from_contours(contours_yx, r, fix_invalid_geometries)
    polygon_contours_xy_copy = copy.copy(polygon_contours_xy)

    #this loop finds if there are holes. If there are, the polygons are ammended
    if number_of_holes != 0:
        for k in range(number_of_holes):
            hole = polygon_contours_xy_copy[hole_posistion_indices[k]]

            for j in range(len(contour_ids)):
                if j in hole_posistion_indices:
                    continue

                contour = polygon_contours_xy_copy[j]
                holey_contour = contour.difference(hole)

                holey_contour = test_filled_hole(holey_contour, hole, contour)
                polygon_contours_xy_copy[j] = holey_contour

        #deleting holes now that they have been dealt with
        polygon_contours_xy_copy = np.delete(polygon_contours_xy_copy, hole_posistion_indices)
        contour_ids = np.delete(contour_ids, hole_posistion_indices)

    return polygon_contours_xy_copy, contour_ids

def get_patches_from_shapely_objects(shapely_objects):
    """
    Any shapely object, from circles, polygons, or point sources
    """
    return [descartes.PolygonPatch(obj) for obj in shapely_objects]


def get_collection_from_patches(patches, c='b', vmin=None, vmax=None, **kwargs):
    """
    Function converts patches into a collection for plotting in matplotlib
    Collections massivly speed up the plotting time.

    Parameters
    ----------
    patches : List of patch objects
        List containing patches
    c : matplotlib color or a list of values to colourise, optional
        c is how to colourise the patches. A single matplotlib colour will assign
        all patches to have that colour. To colour the patches by a parameter
        i.e., their ID number, `c` is instead a list containing the values to
        assign to each patch. Use the kwarg cmap to assign a colourmap.
        The default is 'b'.
    vmin : int, optional
        If each patch is to be colourise separatly, this defines the minimum value
        of the colourmap. The default is None.
    vmax : int, optional
        If each patch is to be colourise sepratly, this defines the maximum value
        of the colourmap. The default is None.
    **kwargs :
        extra parameters for `PatchCollection` to affect the patch properties

    Returns
    -------
    collection : obj
        A collection object

    """
    if np.isscalar(c): #isinstance(c, str)
        kwargs.setdefault('color', c)
        c = None

    if 'fc' in kwargs:
        kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs:
        kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs:
        kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs:
        kwargs.setdefault('linewidth', kwargs.pop('lw'))

    collection = PatchCollection(patches, **kwargs)

    if c is not None:
        collection.set_array(np.array(c))
        collection.set_clim(vmin, vmax)

    return collection

#utility functions for flattening your lists and checking the ID's match up
def flatten_list(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]

def test_ind_matches_id(id_values, amount_of_contours):
    test_array = np.arange(1, amount_of_contours+1)
    if not all(test_array == id_values):
        raise ValueError('Ids do not all match up with indices. Code'\
                         'will only work if the ids match the index')

def test_filled_hole(holey_contour, hole, contour):
    """
    TODO: This function doesn't work properly for complex cases. For now I've
    commented out where it should be used and will fix at another time.
    NB: Might have fixed, need to test

    # Funcion for keeping contours that fall within a 'hole' within a
    # mask i.e:
    #    [ [1,1,1,1,1],
    #      [1,0,0,0,1],
    #      [1,0,1,0,1],
    #      [1,0,0,0,1],
    #      [1,1,1,1,1],
    #    ]
    #If this happens, this function makes sure that the contour isn't removed
    # the_anti_holes = []
    """
    if holey_contour.difference(hole).area == 0:
        # print('A filled in a hole!')
        # the_anti_holes.append(j)
        return hole.difference(contour)
    else:
        return holey_contour
