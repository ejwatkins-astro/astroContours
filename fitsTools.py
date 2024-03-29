# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 15:43:22 2020

@author: Liz_J
"""

from astropy.wcs import WCS
from astropy.io import fits
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import LogLocator
import matplotlib.gridspec as gridspec

from matplotlib import ticker
from matplotlib.ticker import AutoMinorLocator

import matplotlib.pyplot as plt
import numpy as np
from functools import reduce


def remove_axis_from_fits_header(header, axis=3):
    """
    This function takes in a astropy fits header and removes
    dimention keys for the axis entered in the function. If no axis is
    entered, the 3rd axis is assumed the one that needs to be removed
    Returns the adjusted header
    """

    for i in range(len(header)-1,-1,-1):
        if str(axis) in list(header.keys())[i]:
            del header[i]
    header['NAXIS'] = axis-1
    try:
        del header['WCSAXES']
    except KeyError:
        pass

    return header

def get_fits_file_info(fits_file_name, *header_info) :

    """
    Parameters: String
                *arg strings \n
    Returns: 2D  Numpy array, *args mulitple type, mainly Floats
    ----------------
    This function opens up fits files into a 2D Numpy array and loads the
    header information. If only specific header information is wanted, this
    function accepts args as strings to output the corresponding header
    infomation
    """

    with fits.open('%s' % fits_file_name) as hdulist:     # importing fits file

        data = hdulist[0].data
        if header_info == ():
            header_info_array = hdulist[0].header

        else:
            header_info_array = []
            for i in range( len(header_info) ):
                header_info_array.append( hdulist[0].header[header_info[i]])

        hdulist.close()

    return data, header_info_array

def get_fits_file_info_e(fits_file_name, extention=0, *header_info) :

    """
    Parameters: String
                *arg strings \n
    Returns: 2D  Numpy array, *args mulitple type, mainly Floats
    ----------------
    This function opens up fits files into a 2D Numpy array and loads the
    header information. If only specific header information is wanted, this
    function accepts args as strings to output the corresponding header
    infomation
    """

    with fits.open(fits_file_name) as hdulist:     # importing fits file

        data = hdulist[extention].data
        if header_info == ():
            header_info_array = hdulist[extention].header

        else:
            header_info_array = []
            for i in range( len(header_info) ):
                header_info_array.append( hdulist[extention].header[header_info[i]])
        hdulist.close()

    return data, header_info_array

def create_gsgrid(figure_number, plot_gridding):
    fig = plt.figure(figure_number)
    #first number is the amount of rows, second is the amount of columns
    gs = gridspec.GridSpec(*plot_gridding)
    return fig, gs

def add_wcs_subplot(fits_header, gs, subplot_pos, fig=None,  subplot_size = (1,1), **kwargs):

    wcs_info = WCS(fits_header)
    ax = fig.add_subplot(gs[subplot_pos[0]:subplot_pos[0] + subplot_size[0], subplot_pos[1]:subplot_pos[1] + subplot_size[1]], projection=wcs_info, **kwargs)

    return ax


def set_wcs_type(header=None, wsc_type=None):

    if wsc_type is None:
        sky_type = header['CTYPE1'].lower()
        if 'glon' in sky_type or 'glat' in sky_type:
            # celesitial_type=['glon','glat']
            # wsc_type = 'galactic'
            set_wcs_grid = set_galactic_grid

        elif 'ra' in sky_type or 'dec' in sky_type:
            # celesitial_type=['RA','DEC']
            #wsc_type = 'fk5'
            set_wcs_grid = set_equatorial_grid
        else:
            raise ValueError('Can only accept \"galactic\" or \"fk5\"')

    else:
        if wsc_type.lower() == 'galactic' or wsc_type.lower() == 'gal':
            # celesitial_type=['glon','glat']
            set_wcs_grid = set_galactic_grid

        elif wsc_type.lower() == 'fk5' or wsc_type.lower() == 'j2000':
            # celesitial_type=['RA','DEC']
            set_wcs_grid = set_equatorial_grid
        else:
            raise NameError('Can only accept \"galactic\" or \"fk5\"')

    return set_wcs_grid

def subfigure_wsc(figure_number, fits_headers, fontsize='medium',\
                  tick_color='black', plot_gridding = (1,1), \
                  axis_on='all-b', wsc_type = None, sharexy=False):
    """
    This is a wraper function that creates a figure for wcs coordinates
    and allows for sub figures. The Subfigures all have equal gridding
    User needs to adjust the wsc type to match the header. Sky coordinates
    must be the same type i.e., all galactic, or all Ra and dec etc.
    To tweak any ticks, labels ect, see the documentation here:
    http://wcsaxes.readthedocs.io/en/latest/ticks_labels_grid.html#
    """
    #if just 1 plot, migh enter header not in a list format.
    if not isinstance(fits_headers, list):
        fits_headers = [fits_headers]
    set_wcs_grid = set_wcs_type(fits_headers[0], wsc_type)

    fig, gs = create_gsgrid(figure_number, plot_gridding)

    amount_of_plots = plot_gridding[0] * plot_gridding[1]
    axs = []
    x_axis = []
    y_axis = []
    axis_on_bool, tickloc = determine_axis_on_off(axis_on, plot_gridding, amount_of_plots)

    for i in range(amount_of_plots):
        try:
            if sharexy and i>0:
                ax =add_wcs_subplot(fits_headers[i], gs, (i//plot_gridding[1], i%plot_gridding[1]), fig=fig, sharex=ax, sharey=ax)
            else:
                ax =add_wcs_subplot(fits_headers[i], gs, (i//plot_gridding[1], i%plot_gridding[1]), fig=fig)

            y_ax, x_ax = set_wcs_grid(ax, fontsize=fontsize, tick_color=tick_color, tick_location=tickloc, vertical_subplot=axis_on_bool[i,0], horizontal_subplot=axis_on_bool[i,1])
            y_axis.append(y_ax)
            x_axis.append(x_ax)


        except IndexError: #if I want to plot non wcs plots. Assumes are bottom
            ax = fig.add_subplot(gs[i//plot_gridding[1], i%plot_gridding[1]])
            minor_tickers(ax)

        axs.append(ax)

    return x_axis, y_axis, axs, fig


def determine_axis_on_off(axis_on, plot_gridding, amount_of_plots):
    """
    This function says when a sublot grid should have axis values on or off
    """
    axis_on_bool_array = np.ones([amount_of_plots, 2], dtype=bool) # 0=x-axis, 1=y-axis

    if axis_on == 'all-b': #axis values around boarders, left, bottom
        axis_on_bool_array[::plot_gridding[1], 1] = False
        axis_on_bool_array[-plot_gridding[1]:, 0] = False

    if axis_on == 'all-t': #axis values around boarders, left, top
        axis_on_bool_array[::plot_gridding[1], 1] = False
        axis_on_bool_array[:plot_gridding[1], 0] = False


    if axis_on == 'one-b': # one set of axis values, on bottom
         axis_on_bool_array[-plot_gridding[1]] = False

    if axis_on == 'one-t': # one set of axis values, on top
         axis_on_bool_array[0] = False

    return axis_on_bool_array, axis_on[-1]

def set_galactic_grid(ax_object, fontsize='medium', celesitial_type=['glon','glat'], tick_color = 'black', vertical_subplot=False, horizontal_subplot=False, rotation=90, tick_location='b'):
    """
        This functionsets the wsc axis of an instance of a matplotlib axis
    """

    lon = ax_object.coords[celesitial_type[0]]
    lat = ax_object.coords[celesitial_type[1]]

    lon.set_axislabel_position(tick_location)
    lon.set_ticklabel_position(tick_location)

    if vertical_subplot:
        lon.set_ticklabel_visible(False)
    else:
        lon.set_axislabel( 'Galactic Longitude', fontsize=fontsize)

    if horizontal_subplot:
        lat.set_ticklabel_visible(False)
    else:
        lat.set_axislabel( 'Galactic Latitude', fontsize=fontsize, rotation=rotation)


    set_tickers(lon, fontsize, tick_color)
    set_tickers(lat, fontsize, tick_color, rotation=rotation)

    lon.set_major_formatter('d.dd')
    lat.set_major_formatter('d.dd')

    return lat, lon

def set_tickers(axis_object, fontsize='medium',color='black', rotation=0):
    """
    For wcs coordinates, this function sets them in a scientific way
    """
    from astropy import units as u # Move to top, lazy

    axis_object.set_major_formatter('d.dd')#u'%.2f\u00B0')#''d.dd')#u'\u00B0'))#'d:mm') # so only shows

    axis_object.set_ticks(color=color, size = 5, width=1.05)
    axis_object.set_ticklabel(fontsize=fontsize, exclude_overlapping=True, rotation=rotation)
    try:
        # axis_object.set_minor_frequency(5)
        axis_object.display_minor_ticks(True)
    except IndexError:
        print('Not enough ticks for minor ticks!')


def set_equatorial_grid(ax_object, fontsize='medium', celesitial_type=['RA','DEC'], tick_color='black', vertical_subplot=False, horizontal_subplot=False, rotation=90, tick_location='b'):

    RA = ax_object.coords[celesitial_type[0]]
    DEC = ax_object.coords[celesitial_type[1]]

    RA.set_axislabel_position(tick_location)
    RA.set_ticklabel_position(tick_location)

    if vertical_subplot:
        RA.set_ticklabel_visible(False)
    else:
        RA.set_axislabel('RA', fontsize=fontsize, minpad=0.5)

    if horizontal_subplot:
        DEC.set_ticklabel_visible(False)
    else:
        DEC.set_axislabel('DEC', fontsize=fontsize, minpad=0.5, rotation=rotation)

    set_tickers(RA, fontsize, tick_color)
    set_tickers(DEC, fontsize, tick_color, rotation=rotation)

    RA.set_major_formatter('dd:mm:ss.s')
    DEC.set_major_formatter('dd:mm:ss.s')


    return RA, DEC

def make_some_grid(amount):
    n = amount
    factors = list(set(reduce(list.__add__, ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0))))
    if n <=3:
        return factors
    if len(factors) < 3:
        n += 1
        factors = list(set(reduce(list.__add__, ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0))))

    num_factors = len(factors)
    if num_factors%2 == 1:
        grid = (num_factors-1)/2
        return [factors[grid],factors[grid]]
    else:
        grid = (num_factors)/2
        return [factors[grid-1],factors[grid]]

def remove_x_ticker_labels_for_subplots(*axes_objects):
    """
     Removes overlapping tickers from x
    """

    plt.setp([ax.get_xticklabels() for ax in axes_objects[:-1]], visible=False)

def remove_overlapping_tickers_for_vertical_subplots(bin_adjust=0, *axes_objects):
    '''
    This removes overlapping tickerlabels from subplots which have no hspace
    '''

    remove_x_ticker_labels_for_subplots(*axes_objects)

    No_of_bins = len(axes_objects[0].get_xticklabels())
    for ax in axes_objects[:]:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=No_of_bins-bin_adjust, prune='upper'))

def remove_y_ticker_labels_for_subplots(*axes_objects):
    """
    Removes overlapping tickers from y
    """

    plt.setp([ax.get_yticklabels() for ax in axes_objects[1:]], visible=False)

def remove_overlapping_tickers_for_horizontal_subplots(bin_adjust=0, *axes_objects):
    """
    This removes overlapping tickerlabels from subplots which have no hspace.
    """

    remove_y_ticker_labels_for_subplots(*axes_objects)

    No_of_bins = len(axes_objects[0].get_yticklabels())
    for ax in axes_objects[1:]:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=No_of_bins, prune='lower'))

def scientific_tickers(ax_object, **kwargs):
    """
    Requires the minor tickers to already be located
    """

    ax_object.tick_params(which='major', length=6,  **kwargs)
    ax_object.tick_params(which='minor', length=4)#,  **kwargs)

def minor_tickers(ax_object=None, **kwargs):

    if ax_object is None:
        ax_object = plt.gca()

    minorLocator = AutoMinorLocator()

    log_y = ax_object.get_yscale() == 'log'
    log_x = ax_object.get_xscale() == 'log'

    if log_y and log_x :
        pass

    elif log_y and not log_x:
        ax_object.xaxis.set_minor_locator(minorLocator)

    elif not log_y and log_x:
        ax_object.yaxis.set_minor_locator(minorLocator)

    else:
        ax_object.xaxis.set_minor_locator(minorLocator)
        minorLocator = AutoMinorLocator()
        ax_object.yaxis.set_minor_locator(minorLocator)

    scientific_tickers(ax_object, **kwargs)

def colorbar_subplot(ax, im, label='', fontsize='medium',  logged=False, width=0.02, nticks=4, use_gridspec=False, **kwargs):

    fig = plt.gcf()
    cax = fig.add_axes([ax.get_position().x1-0.001, ax.get_position().y0, width, ax.get_position().height])

    if logged:
        cbar = plt.colorbar(im, cax=cax, pad=0, ticks = LogLocator(subs=range(10)), use_gridspec=use_gridspec)
    else:
        cbar = plt.colorbar(im, cax=cax, pad=0, use_gridspec=use_gridspec)
        tick_locator = ticker.MaxNLocator(nbins=nticks)
        cbar.locator = tick_locator

    cbar.ax.tick_params(labelsize=fontsize, length=4, width=1.1, color='k', rotation=90)
    if label != '':
        cbar.set_label(label,fontsize=fontsize, labelpad=7)
    cbar.ax.tick_params(which='minor', color='k', rotation=90, length=2, width=1.1)

    cbar.update_ticks()
    if not logged:
        cbar.ax.minorticks_on()

    return cbar

def colorbar_top(ax, im, label, fontsize='medium',  logged=False, **kwargs):


    fig = plt.gcf()
    cax = fig.add_axes([ax.get_position().x0, ax.get_position().y1-0.02, ax.get_position().width, 0.02])

    if logged:
        cbar = plt.colorbar(im, cax=cax, ticks = LogLocator(subs=range(10)), orientation='horizontal')
    else:
        cbar = plt.colorbar(im, cax=cax, orientation='horizontal', use_gridspec=False)
        # tick_locator = ticker.MaxNLocator(nbins=4)
        # cbar.locator = tick_locator

    cbar.set_alpha(1)
    cbar.draw_all()

    cax.xaxis.tick_top()
    cbar.ax.xaxis.set_label_position('top')

    cbar.ax.tick_params(labelsize=fontsize, length=8, width=1.1, color='k')
    cbar.set_label(label,fontsize=fontsize)
    cbar.ax.tick_params(which='minor', color='k', length=4, width=1.1)


    cbar.update_ticks()
    if not logged:
        cbar.ax.minorticks_on()

    return cbar