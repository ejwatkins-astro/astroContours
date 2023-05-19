# -*- coding: utf-8 -*-
"""
Script for handling contours (a list of boundaries in rows and columns) and
ds9 regions
"""

from astropy import units as u
import regions
from astropy.coordinates import SkyCoord
import numpy as np
from functools import lru_cache
import numpy as np
from astropy import wcs

from .. import fitsTools as mf
from .. import maskContoursUtility as m2c

# import mask_into_contours as cont

accepted_ds9_region_types = {
    'circle_sky': regions.shapes.circle.CircleSkyRegion,
    'circle_pix': regions.shapes.circle.CirclePixelRegion,
    'ellipse_sky':regions.shapes.ellipse.EllipseSkyRegion,
    'ellipse_pix':regions.shapes.ellipse.EllipsePixelRegion,
    'polygon_sky':regions.shapes.polygon.PolygonSkyRegion,
    'polygon_pix':regions.shapes.polygon.PolygonPixelRegion
}



def get_ds9_region_type(region):
    for key in accepted_ds9_region_types.keys():
        if accepted_ds9_region_types[key] == type(region):
            region_type = str(key)


class ds9_regions_contours():#cont.Contours):

    def __init__(self, ds9_region_file, header=None):

        convert_hst_reg_to_standard_reg(ds9_region_file)
        self.ds9_regions = regions.read_ds9(ds9_region_file)

        if header is not None:
            header = mf.remove_axis_from_fits_header(header, 4)
            header = mf.remove_axis_from_fits_header(header, 3)
            self.wcs_into = wcs.WCS(header)
        self.header = header

    def get_contours_wcs_fk5(self):
        return [get_contours_fk5_from_ds9_region(reg) for reg in self.ds9_regions]

    @lru_cache # caching result to speed up any potential following uses of this
    def get_contours_yx(self):
        return [get_contours_yx_from_ds9_region(reg, self.wcs_into) for reg in self.ds9_regions]

    def get_region_IDs_from_image(self, masked_image):
        """
        The masked image must be the same pixel scale as the header
        """
        contours_yx = self.get_contours_yx()
        return [m2c.get_id_of_contour(cont_yx[::-1], masked_image) for cont_yx in contours_yx]

    def get_region_IDs_from_meta_info(self, id_key='label'):
        return [get_ds9_meta_info(region, id_key) for region in self.ds9_regions]


class contours_to_ds9_regions():
    """
    Takes in wsc contours and saves them as ds9 regions
    """

    def __init__(self, contours, contour_ids=None, wcs_projection='fk5'):


        if contour_ids is None:
            self.contour_ids = [None] * len(contours)
        else:
            self.contour_ids = contour_ids

        self.contours = contours
        self.wcs_projection = wcs_projection


    def get_ds9_regions(self):
        return [contour_wcs_to_ds9_polygon(self.contours[i],
                                           self.contour_ids[i],
                                           frame=self.wcs_projection,
                                           unit='deg') for i in range(len(self.contours))]

def contour_wcs_to_ds9_polygon(contour, contour_id=None, **kwargs):

    # RA=0, DEC=1
    sky_coords = SkyCoord(contour[:,0], contour[:,1], **kwargs) #
    meta_data = {}
    if contour_id is not None:
        meta_data['label'] = str(contour_id)
    polygon_sky = regions.PolygonSkyRegion(vertices=sky_coords, meta=meta_data)
    # else:
    #     polygon_sky = regions.PolygonSkyRegion(vertices=sky_coords)

    return polygon_sky

def get_contours_fk5_from_ds9_region(ds9_region):

    ra_verts = ds9_region.vertices.fk5.ra.value
    dec_verts = ds9_region.vertices.fk5.dec.value

    contour_wsc = np.column_stack((ra_verts, dec_verts))
    return contour_wsc

def get_contours_yx_from_ds9_pixel_region(ds9_pixel_region):

    x_verts, y_verts = ds9_pixel_region.vertices.xy
    contour_yx = np.column_stack((y_verts, x_verts))

    return contour_yx

def get_contours_yx_from_ds9_region(ds9_region, image_wcs):

    ds9_pixel_region = ds9_region.to_pixel(image_wcs)
    x_verts, y_verts = ds9_pixel_region.vertices.xy
    contour_yx = np.column_stack((y_verts, x_verts))

    return contour_yx

def get_ds9_meta_info(region, key='label'):
    return region.meta[key]


def convert_hst_reg_to_standard_reg(hst_ds9_regions, wcs_projection='fk5', region_type='polygon'):
    """
    For the region parser I use, the coordinate system has to be specified
    """
    with open (hst_ds9_regions, "r") as myfile:
        lines = myfile.readlines()
    for i in range(len(lines)):
        if region_type in lines[i]:
            has_coordinate_system = wcs_projection in lines[i-1]
            break

    if not has_coordinate_system:
        with open (hst_ds9_regions, "w") as myfile:
            for j in range(i):
                myfile.write(lines[j])
            myfile.write(wcs_projection + '\n')

            for line in lines[i:]:
                myfile.write(line)


# DRIVE_PATH = 'C:\\Users\\Liz_J\\Documents\\'
# PROJECT_PATH = DRIVE_PATH + 'project1-postdoc\\'

# GALAXY = 'NGC1365'

# hst_regions_name = '_phangshst_associations_v_ws32pc_v1p1.reg'
# hst_catalogue_location = '\\hst_catalogue\\'

# filepath = DRIVE_PATH + GALAXY + hst_catalogue_location + GALAXY.lower() + hst_regions_name
# filepath_test = DRIVE_PATH + GALAXY + hst_catalogue_location + GALAXY.lower() + '_test' + hst_regions_name


# ds9 = ds9_regions_contours(filepath)
# contours_wcs = ds9.get_contours_wcs_fk5()


# with open (filepath, "r") as myfile:
#     data = myfile.readlines()

# with open (filepath_test, "w") as myfile:
# # Region file format: DS9 version 4.1
#     if 'fk5' not in data[1]:
#         myfile.write('# Region file format: DS9 version 4.1\n')
#         myfile.write('fk5')
#         myfile.write(data[0])

#     if ',' not in data[2]:
#         for string in data[2:]:
#             string = string[1:]
#             string = string.replace('d ', ',')
#             string = string.replace('n ','n(')
#             string = string[:-2] + ')\n'
#             myfile.write(string)

