# -*- coding: utf-8 -*-
"""
TO FINISH
"""

class Patch_astroContours:
    def __init__(self, data, levels, use_world_coords=None, header=None, convert_world=False, **kwargs):

        self.data = data
        self.levels = levels
        self.use_world_coords = use_world_coords
        self.header = header
        self.convert_world = convert_world
        if self.header is not None:
            self.world_coord_type = get_world_type(self.header)

    def get_contours(self):

        header = self.header
        c = maskContours.Levels_to_contours(data=self.data, levels=self.levels)

        if self.use_world_coords and header is not None:
            c.set_header(header)
            if self.convert_world:
                # new_header = convert_header_to_new_world(header, coord_types[self.world_coord_type])
                contours, hole_locations = c.convert_contours_wcs_to_new_world(config.coord_conversion[self.world_coord_type])
            else:
                contours, hole_locations = c.get_contours_wcs()
        else:
            contours, hole_locations = c.get_contours_yx()

        return contours, hole_locations

def get_colours_from_cmap(values, vmin, vmax, cmap):

    norm = colors.Normalize(vmin, vmax)
    colour_from_cmap = []
    for val in values:
        colour_from_cmap.append(cmap(norm(val)))

    return colour_from_cmap