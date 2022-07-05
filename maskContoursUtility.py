# -*- coding: utf-8 -*-
"""
Script contain some wrappers around functions that convert masked-labled images
into contours (a list of boundaries in rows and columns), and functions to
convert coorinates to world or pixel frames
"""

from skimage import measure
from skimage.draw import polygon
from collections import Counter
import numpy as np
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
import copy as copy


def convert_world2pixel(paired_world_coordinates_of_boarders, header):

    """
    Load the WCS information from a fits header, and use it
    to convert pixel coordinates to world coordinates.
    """
    w = wcs.WCS(header)
    paired_pixel_coordinates_of_boarders = w.all_world2pix(paired_world_coordinates_of_boarders, 0)
    return paired_pixel_coordinates_of_boarders

def convert_pixel2world(paired_pixel_coordinates_of_boarders, header):

    """
    Load the WCS information from a fits header, and use it
    to convert world coordinates to pixel coordinates.
    """
    w = wcs.WCS(header)
    paired_world_coordinates_of_boarders = w.all_pix2world(paired_pixel_coordinates_of_boarders, 0)
    return paired_world_coordinates_of_boarders

def counter_2d(data_2d, y, x):
    """
    This function counts the occurance of a number in 2d images
    """
    c = Counter(data_2d[y, x].flatten())

    ### first row is the ordered counts, second row is the occcurance.
    #   most counted is first
    most_common = np.column_stack(c.most_common())[0]

    return most_common

def get_contours_from_touching_mask(image, contour_dilation):
    """
    Finds contours for touching masks. Needs to loop through masked regions
    which is a bit slower
    """
    # finds properties of labeled masked images
    props = measure.regionprops(image.astype(int))

    all_contours = []
    for prop in props:
        blank_image = np.zeros_like(image)
        coords = prop.coords
        coords_y = coords[:,0]
        coords_x = coords[:,1]
        blank_image[coords_y, coords_x] = 1

        contours = measure.find_contours(blank_image, contour_dilation)
        for cont in contours:
            all_contours.append(cont)

    return all_contours

def get_contours(labeled_image, contour_dilation, get_contour_id=False, touching_masks=False):
    """
    Get the co-ordinates of the contour and the id number it corresponds to from
    a labeled image. Holes are always given after polygon, but because the insides
    of holes are zeros, they have an id of zero.
    """
    contour_id = []
    # this finds the contours as y, x verticies
    if touching_masks:
        contours = get_contours_from_touching_mask(image=labeled_image, contour_dilation=contour_dilation)
    else:
        contours = measure.find_contours(labeled_image, contour_dilation)

    # to get the id number for each contour, I grab the co-ords within the contour
    # and find the largest number (which will be the contour id)
    if get_contour_id:
        contour_id = [get_id_of_contour(pixel_coords, labeled_image) for pixel_coords in contours]

    return contours,  contour_id

def get_id_of_contour(pixel_coords, labeled_image):
    """
    For a contour given in x, y coordinates, and a maksed image
    this function will find id contained within the contour
    """
    y_indices_within, x_indices_within = polygon(pixel_coords[:,0], pixel_coords[:,1])

    y_indices_within, x_indices_within = correct_pixel_width_contour(pixel_coords,
                                                                     y_indices_within,
                                                                     x_indices_within)

    most_common = counter_2d(labeled_image, y_indices_within, x_indices_within)
    #ignores any holes
    most_common = most_common[most_common>0]
    try:
        num = most_common[0]
    except IndexError: # contour is only a hole
        num = 0
    return int(num)

def correct_pixel_width_contour(pixel_coords, y_indices_within, x_indices_within):
    """
    Sometimes the contour might be a single pixel accross, which `polygon`
    cannot handle. This function outputs the correct indices of the contour
    when this happens

    """
    if len(pixel_coords[:,0])>0 and len(y_indices_within) ==0:
        y_indices_within, x_indices_within = pixel_coords[:,0].astype(int), pixel_coords[:,1].astype(int)

    return y_indices_within, x_indices_within

def get_indices_within_contour(pixel_coords):
    """
    This function gets all the x,y coordinates within a single contour
    """

    y_within, x_within = polygon(pixel_coords[:,1], pixel_coords[:,0])

    #these are contours that are only 1 pixel in length and so form a line
    y_within, x_within = correct_pixel_width_contour(pixel_coords, y_within, x_within)

    return y_within, x_within

def where_within_image_boundary(y_within, x_within, image_shape):
    """
    Function removes coordinates that are out of bounds of a given image
    shape
    """
    x_p = (x_within < image_shape[1]) & (x_within >= 0) & \
          (y_within < image_shape[0]) & (y_within >= 0)

    return x_p

def get_valid_indices_within_region(pixel_coords, image):
    """
    Function return all the valid coordinates found within a contour
    """
    #pixels inside contour
    y_indices_within, x_indices_within  = get_indices_within_contour(pixel_coords)

    #ignoring coordinates that are out of bounds with the image
    im_shape = image.shape
    x_p = where_within_image_boundary(y_indices_within, x_indices_within, im_shape)
    return y_indices_within[x_p], x_indices_within[x_p]


def set_value_to_image(pixel_coords, image, value):
    """
    Given a set of pixel coordinates, this function assigns the given value to
    the given 2d image
    """
    y_indices_within, x_indices_within = get_valid_indices_within_region(pixel_coords, image)

    image[y_indices_within, x_indices_within] = value

    return image

def create_masks_from_wcs_contours(contours_WCS, contourIDs, header, image, binary_mask_out=False):
    """
    Given a set of contours in world coordinates, this function produces
    the mask of these contours for different world coordinates systems
    """
    image = np.zeros_like(image, dtype=int)

    hole_location = []
    for j in range(len(contourIDs)):
        if contourIDs[j] == 0:
            hole_location.append(j)
            continue

        pixel_coords = convert_world2pixel(contours_WCS[j], header)
        image = set_value_to_image(pixel_coords, image, contourIDs[j])

    for hole in hole_location:
        pixel_coords = convert_world2pixel(contours_WCS[hole], header)
        image = set_value_to_image(pixel_coords, image, 0)

    #if just a binary mask of 1's and 0's wanted
    if binary_mask_out:
        image[image>0] =1

    if binary_mask_out == 'both':

        return image, make_bitmap(image)
    else:
        return image

def create_masks_from_yx_contours(contours_yx, contourIDs, image, binary_mask_out=False):
    """
    Given a set of contours in world coordinates, this function produces
    the mask of these contours for different world coordinates systems
    """
    image = np.zeros_like(image, dtype=int)

    hole_location = []
    for j in range(len(contourIDs)):
        if contourIDs[j] == 0:
            hole_location.append(j)
            continue
        pixel_coords = contours_yx[j][:,::-1]

        image = set_value_to_image(pixel_coords, image, contourIDs[j])

    for hole in hole_location:
        image = set_value_to_image(pixel_coords, image, 0)

    #if just a binary mask of 1's and 0's wanted
    if binary_mask_out:
        image[image>0] =1

    if binary_mask_out == 'both':

        return image, make_bitmap(image)
    else:
        return image

def make_bitmap(image, value=1):
    """
    Function masks all masked pixels with a single value
    """

    image1 = copy.copy(image)
    y,x = np.where(image1 > 0)
    image1[y,x] = value
    return image1

def convert_between_worlds(world_coordinates, current_world, new_world):
    """
    This method converts the world co-ordinate types from equatorial to
    galactic

    """

    current_world_x, current_world_y = world_coordinates[:,0], world_coordinates[:,1]

    current_sky_coords = SkyCoord(current_world_x*u.degree, current_world_y*u.degree, frame=current_world)
    new_sky_coords = current_sky_coords.transform_to(new_world).to_string()
    new_world_x, new_world_y = map(list, zip(*[coords.split(' ') for coords in new_sky_coords]))

    new_world_coordinates = np.column_stack((new_world_x, new_world_y)).astype(float)

    return new_world_coordinates

def find_hole_contours(contour):
    """
    Contours will wind counter-clockwise (i.e. in ‘positive orientation’)
    around islands of low-valued pixels.
    A value >= 0 indicates a counter-clockwise oriented polygon.
    """
    contour = np.append(contour, [[1,1]], axis=0)
    y_vertices, x_vertices = contour[:,0], contour[:,1]

    orientation = sum(x_vertices[i] * (y_vertices[i+1] - y_vertices[i-1]) for i in range(1, len(y_vertices)-1)) / 2.0

    value = 0 if orientation >= 0 else 1

    return value