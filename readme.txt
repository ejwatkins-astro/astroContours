To use contour codes, download all of astroContours into a folder called astroContours

Then in the same directory as the folder location, you can import astroContours

Code requirements
Python 3

Package requirements likely already met:
numpy
scipy
matplotlib
astropy
skimage

Package requirements I expect you'll need to install:

geopandas: awesome package for comparing, joining, analysing and visulising large shape catalogues (including polygons created from a mask catalogue

shapely:   package to convert contours into polygon objects for geopandas

descartes: package for converting shapley objects into matplotib patches

colorcet:  A huge list of perceptually uniform matplotlib colourmaps

regions:   astropy affiliated package for creating, opening, analysing and visualising and ds9 regions

pickle:    a way to save data

json       a way to save database data (usually dicts) in human readable format. It normally prefered over pickle for security reasons.

Notes:
Mask cataogue needs to have background values of zero, with the rest labelled using positive masks

Here is an example using this for MUSE and JWST data
```python
# -*- coding: utf-8 -*-
"""
@author: Liz_J
"""

import numpy as np

from astropy import wcs
from astropy.io import fits
from astropy.table import Table

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import colorcet as cc

#self made codes and packages
from astroContours.convert import contourPolygons as poly

#Makes tick labels go inside
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['mathtext.default'] = 'regular'


def fix_tight_layout(x_length, y_length, gridding, wspace=0, hspace=None, figscale=3, bottom=0.05, top=0.13):
    """
    For a grid of subplots with the same induvidual aspect ratios
    (so their x and y lims can different, but their x_length/y_length
     are identical), calculates the correct figure size needed to have the exact amount/correct about of white space, including no white space

    Parameters
    ----------
    x_length : float
        Length of the x axis.
    y_length : TYPE
        Length of the y axis.
    gridding : 2-array
        How the subplots are gridded. `gridding[0]` are the number
        of rows, `gridding[1]` are the the number of columns
    wspace : float, optional
        The white space between each column subplot. The default is 0.
    hspace : float or None, optional
        The white space between each row subplot. The default is None.
        If None is used, the white space is set to equal, in aspect
        ratio, the wspace given
    figscale : float, optional
        Factor by which to increase the figure size. The default is 3,
        where 3 refers to the height (in inches) of the figure

    Returns
    -------
    figsize : 2-array
        The width and height needed to get the correct white space for
        the figure.
    subplot_params : 6-array
        The adjustment parameters to use in matplotlib.pyplot.subplot_adjust
        to remove or add the correct amount of white space

    """

    aspect = x_length / y_length
    # n=number of rows, m =numberof columns
    row, col = gridding

    top = 1 - bottom
    right = 1 - left

    fig_aspect = (1 - bottom - (1 - top)) / (1 - left - (1 - right))
    if hspace is None:
        hspace = wspace/aspect

    #fix the figure height
    fig_height = figscale # inch
    fig_width = fig_height * fig_aspect * aspect * (col + (col - 1) * wspace) / (row + (row - 1) * hspace)

    if isinstance(fig_width, np.ndarray):
        fig_width = fig_width[0]

    if isinstance(fig_height, np.ndarray):
        fig_height = fig_height[0]

    figsize = [fig_width, fig_height]
    subplot_params = [left, bottom, right, top, wspace, hspace]

    return figsize, subplot_params

# plt.close('all')
#segagated colourmap:
cmap = cc.m_glasbey_light
cmap.set_bad('black')

DRIVE_PATH =
PROJECT_PATH =


GALAXY = 'NGC7496'

filename = '_IMAGE_FOV_WFI_Hasub-dr2-nati.fits'
mask_filename = '_HIIreg_mask-dr2.fits'
# hst_filename = '_'

bands = ['770']#, '1000', '1130', '2100']

for band in bands:
    jwst_root = 'jw02107-o038_t019_miri_f%sw_i2d' % band
    jwst_filename = jwst_root + '.fits'

    # for GALAXY in galaxy_names:
    #     print(GALAXY)
    mask_filepath = PROJECT_PATH + GALAXY + '\\%s' % (GALAXY + mask_filename)
    image_filepath = PROJECT_PATH + GALAXY + filename
    neb_cat_filepath = PROJECT_PATH + 'Nebulae_catalogue_v2.fits'

    # hst_filepath = DRIVE_PATH + '\\%s\\data1\\'
    jwst_filepath = DRIVE_PATH + GALAXY + '\\JWST\\' + jwst_filename


    with fits.open(image_filepath) as hdulist:
        image_data = hdulist[1].data
        header = hdulist[1].header

    t = Table.read(neb_cat_filepath, format='fits')
    dis = t['HA6562_SIGMA'][np.where(t['gal_name'] == GALAXY)[0]].value
    ids_tab = t['region_ID'][np.where(t['gal_name'] == GALAXY)[0]].value

    # --------------EXAMPLE: plotting a 2d mask catalogue onto new projection

    with fits.open(jwst_filepath) as hdulist:
        data_jwst = hdulist[1].data
        header_jwst = hdulist[1].header
    data_jwst[data_jwst==0] = np.nan

    #touching_masks=True makes this a bit slow, if you know the masks do not touch
    #set this to False and it will be much faster. `add_to_mask` needed for 
    #Muse neb catalogue since it starts at -1, with first mask=0
    p_obj = poly.ContoursAsPolygons.from_file(mask_filepath, hdu_i=0, touching_masks=True, convert_wcs=True, new_header=header_jwst, add_to_mask=1)

    #make the patches
    patch_masks, patch_id = p_obj.get_patches()

    # ----------Example for using ID colours for patch values: ----------------------

    collection_ids = poly.get_collection_from_patches(patch_masks, c=patch_id, cmap=cmap, alpha=0.5)

    # figsize, subplot_params = fix_tight_layout(header_jwst['NAXIS1'], header_jwst['NAXIS2'], gridding=[1,1], figscale=8)
    plt.figure(GALAXY), #figsize=figsize)
    plt.clf()

    pl, ph = np.nanpercentile(data_jwst, (14, 99.5))
    im = plt.imshow(data_jwst, origin='lower', cmap=plt.cm.Greys, norm=LogNorm(vmin=pl, vmax=ph))


    plt.gca().add_collection(collection_ids)
    plt.title(GALAXY +' f%s' % (band), fontsize='large')

    # plt.subplots_adjust(*subplot_params)


    # plt.savefig(DRIVE_PATH + GALAXY + '\\JWST\\' + jwst_root + '-neb_overlay-Ha.png', box_inches='tight',dpi=400)

    # ----------Example for using dispersion colours for patch values ----

    dis_order = dis[patch_id-1]
    collection_ids = poly.get_collection_from_patches(patch_masks, c=dis_order, cmap=plt.cm.inferno, alpha=0.5, norm=LogNorm(vmax=100))

    plt.figure(GALAXY + ' dis'), #figsize=figsize)
    plt.clf()

    pl, ph = np.nanpercentile(data_jwst, (14, 99.5))
    im = plt.imshow(data_jwst, origin='lower', cmap=plt.cm.Greys, norm=LogNorm(vmin=pl, vmax=ph))


    plt.gca().add_collection(collection_ids)
    plt.title(GALAXY +' f%s' % (band), fontsize='large')

    cbar = plt.colorbar(collection_ids, pad=0, label='Hα velocity dispersion [km/s]')
    cbar.set_alpha(1)
    cbar.draw_all()

    # plt.subplots_adjust(*subplot_params)


    # plt.savefig(DRIVE_PATH + GALAXY + '\\JWST\\' + jwst_root + '-neb_overlay-Ha-dis.png', box_inches='tight',dpi=400)
```


