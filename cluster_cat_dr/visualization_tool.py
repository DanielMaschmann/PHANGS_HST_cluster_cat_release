"""
Tool to visualize PHANGS imaging data
"""

import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground

from cluster_cat_dr import phot_data_access, helper_func

import multicolorfits as mcf
from matplotlib.colors import Normalize, LogNorm

import matplotlib.pyplot as plt
from matplotlib import patheffects


class PhotVisualize(phot_data_access.PhotAccess):
    """
    Class to plot cutouts in multiple bands
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def plot_multi_band_artifacts(self, ra, dec, str_line_1=None, str_line_2=None,
                                  cutout_size=2.0, circle_rad=0.16):
        """
        Function to create multi band panel for a single obejct

        Parameters
        ----------
        ra : float
            coordinate in deg without any units
        dec : float
            coordinate in deg without any units
        str_line_1 : str
            some descriptive string in line 1
        str_line_2 : str
            some descriptive string in line 1
        cutout_size : float
            size of cutout
        circle_rad : float
            size in arc seconds of the circle to mark the identified cluster

        Returns
        -------
        figure : ``matplotlib.pylab.Figure``
            figure
        """
        cutout_dict = self.get_band_cutout_dict(ra_cutout=ra, dec_cutout=dec, cutout_size=cutout_size)

        # get nircam and miri flags
        if self.target_name in self.phangs_nircam_obs_target_list:
            nircam_obs_flag = True
        else:
            nircam_obs_flag = False
        if self.target_name in self.phangs_miri_obs_target_list:
            miri_obs_flag = True
        else:
            miri_obs_flag = False

        # plot figure
        figure = plt.figure(figsize=(29, 17))
        # parameters for the axis alignment
        fontsize = 25
        axis_width = 0.27
        axis_height = 0.27
        axis_space_x = -0.11
        start_x = -0.02
        hst_ax_y = 0.64
        nircam_ax_y = 0.345
        miri_ax_y = 0.05

        # set up axis for each filter
        axis_dict = {}
        for band_index, hst_band in enumerate(self.get_hst_band_list()):
            axis_dict.update({
                'ax_%s' % hst_band: figure.add_axes([start_x + (axis_width + axis_space_x) * band_index, hst_ax_y,
                                                     axis_width, axis_height],
                                                    projection=cutout_dict['%s_img_cutout' % hst_band].wcs)
            })

            if nircam_obs_flag | miri_obs_flag:
                ra_axis_label = ' '
            else:
                ra_axis_label = 'R.A. (2000.0)'

            if band_index == 0:
                self.arr_axis_params(ax=axis_dict['ax_%s' % hst_band],
                                     ra_tick_label=np.invert((nircam_obs_flag | miri_obs_flag)),
                                     ra_axis_label=ra_axis_label, fontsize=fontsize, labelsize=fontsize)
            else:
                self.arr_axis_params(ax=axis_dict['ax_%s' % hst_band],
                                     ra_tick_label=np.invert((nircam_obs_flag | miri_obs_flag)),
                                     dec_tick_label=False, ra_axis_label=ra_axis_label, dec_axis_label=' ',
                                     fontsize=fontsize, labelsize=fontsize)
            axis_dict['ax_%s' % hst_band].set_title(hst_band, fontsize=fontsize)

        for band_index, nircam_band in enumerate(self.get_nircam_band_list()):
            if cutout_dict['%s_img_cutout' % nircam_band].data is not None:
                axis_dict.update({
                    'ax_%s' % nircam_band: figure.add_axes([start_x + (axis_width + axis_space_x) * band_index,
                                                            nircam_ax_y, axis_width, axis_height],
                                                           projection=cutout_dict['%s_img_cutout' % nircam_band].wcs)
                })
                if miri_obs_flag:
                    ra_axis_label = ' '
                else:
                    ra_axis_label = 'R.A. (2000.0)'

                if band_index == 0:
                    self.arr_axis_params(ax=axis_dict['ax_%s' % nircam_band],
                                         ra_tick_label=np.invert(miri_obs_flag),
                                         ra_axis_label=ra_axis_label, fontsize=fontsize, labelsize=fontsize)
                else:
                    self.arr_axis_params(ax=axis_dict['ax_%s' % nircam_band],
                                         ra_tick_label=np.invert(miri_obs_flag),
                                         dec_tick_label=False, ra_axis_label=ra_axis_label, dec_axis_label=' ',
                                         fontsize=fontsize, labelsize=fontsize)
                axis_dict['ax_%s' % nircam_band].set_title(nircam_band, fontsize=fontsize)

        for band_index, miri_band in enumerate(self.get_miri_band_list()):
            if cutout_dict['%s_img_cutout' % miri_band].data is not None:
                axis_dict.update({
                    'ax_%s' % miri_band: figure.add_axes([start_x + (axis_width + axis_space_x) * band_index, miri_ax_y,
                                                         axis_width, axis_height],
                                                         projection=cutout_dict['%s_img_cutout' % miri_band].wcs)
                })

                if band_index == 0:
                    self.arr_axis_params(ax=axis_dict['ax_%s' % miri_band], fontsize=fontsize, labelsize=fontsize,
                                         ra_tick_num=2)
                else:
                    self.arr_axis_params(ax=axis_dict['ax_%s' % miri_band], dec_tick_label=False, dec_axis_label=' ',
                                         fontsize=fontsize, labelsize=fontsize, ra_tick_num=2)
                axis_dict['ax_%s' % miri_band].set_title(miri_band, fontsize=fontsize)

        # now plot all bands
        for band in self.get_phangs_band_list(include_hst_ha=False):
            if cutout_dict['%s_img_cutout' % band].data is not None:

                v_min, v_max = self.get_image_scale_with_circle(img_data=cutout_dict['%s_img_cutout' % band].data,
                                                                img_wcs=cutout_dict['%s_img_cutout' % band].wcs,
                                                                ra=ra, dec=dec, circle_rad=circle_rad)
                axis_dict['ax_%s' % band].imshow(cutout_dict['%s_img_cutout' % band].data,
                                                 cmap='Greys', vmin=v_min, vmax=v_max)
                pos = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
                self.plot_circle_on_wcs_img(ax=axis_dict['ax_%s' % band], pos=pos, rad=circle_rad, color='r')

        # now plot all rgb images
        # specify bands
        # hst
        hst_red_band = 'F555W'
        if 'F438W' in self.phangs_hst_obs_band_dict[self.target_name]['wfc3_uvis_observed_bands']:
            hst_green_band = 'F438W'
        else:
            hst_green_band = 'F435W'
        hst_blue_band = 'F336W'
        # nircam
        nircam_red_band = 'F360M'
        nircam_green_band = 'F300M'
        nircam_blue_band = 'F200W'
        # miri
        miri_red_band = 'F1130W'
        miri_green_band = 'F1000W'
        miri_blue_band = 'F770W'

        # get hst rgb images
        axis_dict.update({
                'ax_hst_rgb': figure.add_axes([start_x+(axis_width+axis_space_x)*5, hst_ax_y, axis_width,
                                               axis_height], projection=cutout_dict['%s_img_cutout' % 'F336W'].wcs)
            })

        min_red_hst, max_red_hst = self.get_image_scale_with_circle(
            img_data=cutout_dict['%s_img_cutout' % hst_red_band].data,
            img_wcs=cutout_dict['%s_img_cutout' % hst_red_band].wcs, ra=ra, dec=dec, circle_rad=circle_rad)
        min_green_hst, max_green_hst = self.get_image_scale_with_circle(
            img_data=cutout_dict['%s_img_cutout' % hst_green_band].data,
            img_wcs=cutout_dict['%s_img_cutout' % hst_green_band].wcs, ra=ra, dec=dec, circle_rad=circle_rad)
        min_blue_hst, max_blue_hst = self.get_image_scale_with_circle(
            img_data=cutout_dict['%s_img_cutout' % hst_blue_band].data,
            img_wcs=cutout_dict['%s_img_cutout' % hst_blue_band].wcs, ra=ra, dec=dec, circle_rad=circle_rad)

        hst_rgb = self.get_rgb_img(data_r=cutout_dict['%s_img_cutout' % hst_red_band].data,
                                   data_g=cutout_dict['%s_img_cutout' % hst_green_band].data,
                                   data_b=cutout_dict['%s_img_cutout' % hst_blue_band].data,
                                   min_max_r=(min_red_hst, max_red_hst),
                                   min_max_g=(min_green_hst, max_green_hst),
                                   min_max_b=(min_blue_hst, max_blue_hst),
                                   scaletype_r='abs', scaletype_g='abs', scaletype_b='abs')
        axis_dict['ax_hst_rgb'].imshow(hst_rgb)
        if nircam_obs_flag | miri_obs_flag:
            ra_axis_label = ' '
        else:
            ra_axis_label = 'R.A. (2000.0)'
        self.arr_axis_params(ax=axis_dict['ax_hst_rgb'], ra_tick_label=np.invert((nircam_obs_flag | miri_obs_flag)),
                             dec_tick_label=False, ra_axis_label=ra_axis_label, dec_axis_label=' ',
                             fontsize=fontsize, labelsize=fontsize)
        axis_dict['ax_hst_rgb'].set_title(hst_red_band + '+' + hst_green_band + '+' + hst_blue_band, fontsize=fontsize)

        if nircam_obs_flag:
            if ((cutout_dict['%s_img_cutout' % nircam_red_band].data is not None) &
                    (cutout_dict['%s_img_cutout' % nircam_green_band].data is not None) &
                    (cutout_dict['%s_img_cutout' % nircam_blue_band].data is not None)):
                axis_dict.update({
                        'ax_nircam_rgb': figure.add_axes([start_x+(axis_width+axis_space_x)*5, nircam_ax_y, axis_width,
                                                          axis_height],
                                                         projection=cutout_dict['%s_img_cutout' % nircam_blue_band].wcs)
                    })

                blue_data_nircam = cutout_dict['%s_img_cutout' % nircam_blue_band].data
                wcs_nircam = cutout_dict['%s_img_cutout' % nircam_blue_band].wcs
                green_data_nircam = helper_func.reproject_image(
                    data=cutout_dict['%s_img_cutout' % nircam_green_band].data,
                    wcs=cutout_dict['%s_img_cutout' % nircam_green_band].wcs,
                    new_wcs=wcs_nircam, new_shape=blue_data_nircam.shape)
                red_data_nircam = helper_func.reproject_image(data=cutout_dict['%s_img_cutout' % nircam_red_band].data,
                                                              wcs=cutout_dict['%s_img_cutout' % nircam_red_band].wcs,
                                                              new_wcs=wcs_nircam, new_shape=blue_data_nircam.shape)

                min_red_nircam, max_red_nircam = self.get_image_scale_with_circle(
                    img_data=red_data_nircam,
                    img_wcs=wcs_nircam, ra=ra, dec=dec, circle_rad=circle_rad)
                min_green_nircam, max_green_nircam = self.get_image_scale_with_circle(
                    img_data=green_data_nircam,
                    img_wcs=wcs_nircam, ra=ra, dec=dec, circle_rad=circle_rad)
                min_blue_nircam, max_blue_nircam = self.get_image_scale_with_circle(
                    img_data=blue_data_nircam,
                    img_wcs=wcs_nircam, ra=ra, dec=dec, circle_rad=circle_rad)

                nircam_rgb = self.get_rgb_img(data_r=red_data_nircam,
                                              data_g=green_data_nircam,
                                              data_b=blue_data_nircam,
                                              min_max_r=(min_red_nircam, max_red_nircam),
                                              min_max_g=(min_green_nircam, max_green_nircam),
                                              min_max_b=(min_blue_nircam, max_blue_nircam),
                                              scaletype_r='abs', scaletype_g='abs', scaletype_b='abs')
                axis_dict['ax_nircam_rgb'].imshow(nircam_rgb)
                if miri_obs_flag:
                    ra_axis_label = ' '
                else:
                    ra_axis_label = 'R.A. (2000.0)'
                self.arr_axis_params(ax=axis_dict['ax_nircam_rgb'], ra_tick_label=np.invert(miri_obs_flag),
                                     dec_tick_label=False, ra_axis_label=ra_axis_label, dec_axis_label=' ',
                                     fontsize=fontsize, labelsize=fontsize, ra_tick_num=2)
                axis_dict['ax_nircam_rgb'].set_title('%s+%s+%s' %
                                                     (nircam_red_band, nircam_green_band, nircam_blue_band),
                                                     fontsize=fontsize)

        if miri_obs_flag:
            if ((cutout_dict['%s_img_cutout' % miri_red_band].data is not None) &
                    (cutout_dict['%s_img_cutout' % miri_green_band].data is not None) &
                    (cutout_dict['%s_img_cutout' % miri_blue_band].data is not None)):
                axis_dict.update({
                        'ax_miri_rgb': figure.add_axes([start_x+(axis_width+axis_space_x)*5, miri_ax_y, axis_width,
                                                        axis_height],
                                                       projection=cutout_dict['%s_img_cutout' % miri_blue_band].wcs)
                    })

                blue_data_miri = cutout_dict['%s_img_cutout' % miri_blue_band].data
                wcs_miri = cutout_dict['%s_img_cutout' % miri_blue_band].wcs
                green_data_miri = helper_func.reproject_image(data=cutout_dict['%s_img_cutout' % miri_green_band].data,
                                                              wcs=cutout_dict['%s_img_cutout' % miri_green_band].wcs,
                                                              new_wcs=wcs_miri,
                                                              new_shape=blue_data_miri.shape)
                red_data_miri = helper_func.reproject_image(data=cutout_dict['%s_img_cutout' % miri_red_band].data,
                                                            wcs=cutout_dict['%s_img_cutout' % miri_red_band].wcs,
                                                            new_wcs=wcs_miri,
                                                            new_shape=blue_data_miri.shape)

                min_red_miri, max_red_miri = self.get_image_scale_with_circle(
                    img_data=red_data_miri,
                    img_wcs=wcs_miri, ra=ra, dec=dec, circle_rad=circle_rad)
                min_green_miri, max_green_miri = self.get_image_scale_with_circle(
                    img_data=green_data_miri,
                    img_wcs=wcs_miri, ra=ra, dec=dec, circle_rad=circle_rad)
                min_blue_miri, max_blue_miri = self.get_image_scale_with_circle(
                    img_data=blue_data_miri,
                    img_wcs=wcs_miri, ra=ra, dec=dec, circle_rad=circle_rad)

                miri_rgb = self.get_rgb_img(data_r=red_data_miri, data_g=green_data_miri, data_b=blue_data_miri,
                                            min_max_r=(min_red_miri, max_red_miri),
                                            min_max_g=(min_green_miri, max_green_miri),
                                            min_max_b=(min_blue_miri, max_blue_miri),
                                            scaletype_r='abs', scaletype_g='abs', scaletype_b='abs')
                axis_dict['ax_miri_rgb'].imshow(miri_rgb)
                if miri_obs_flag:
                    ra_axis_label = ' '
                else:
                    ra_axis_label = 'R.A. (2000.0)'
                self.arr_axis_params(ax=axis_dict['ax_miri_rgb'], ra_tick_label=np.invert(miri_obs_flag),
                                     dec_tick_label=False, ra_axis_label=ra_axis_label, dec_axis_label=' ',
                                     fontsize=fontsize, labelsize=fontsize, ra_tick_num=2)
                axis_dict['ax_miri_rgb'].set_title('%s+%s+%s' % (miri_red_band, miri_green_band, miri_blue_band),
                                                   fontsize=fontsize)

        # add descriptive string
        plt.figtext(0.05, 0.97, str_line_1, fontsize=fontsize)
        plt.figtext(0.05, 0.94, str_line_2, fontsize=fontsize)

        return figure


    def plot_overview_zoom_in_panel(self, ra_region, dec_region, include_h_alpha=True, ha_cont_sub_stamp_band='ha',
                                    include_nircam=True,
                                    include_miri=True,
                                    env_cutout_size=(10, 10), circle_rad_region=1.25, circle_rad_obj=0.25,
                                    stamp_cutout_size=(2.5, 2.5), cbar_log_scale=True,
                                    fontsize_large=30, fontsize_small=20,
                                    ):

        # load all bands needed
        band_list = []
        # load HST bands
        band_list += self.get_hst_band_list()
        hst_stamp_band_list = band_list.copy()
        # get BVI filters
        # specify color the bands
        hst_bvi_band_red = 'F814W'
        hst_bvi_band_green = 'F555W'
        if 'F438W' in band_list:
            hst_bvi_band_blue = 'F438W'
        else:
            hst_bvi_band_blue = 'F435W'

        # load ha native and continuum subtracted
        if include_h_alpha:
            band_list += [ha_cont_sub_stamp_band]
            # get HST-H-alpha
            if 'F657N' in self.phangs_hst_ha_obs_band_dict[self.hst_ha_target_name]['ha_observed']:
                hst_habu_band_red = 'F657N'
                band_list += ['F657N']
                hst_ha_stamp_band = 'F657N'
            elif 'F658N' in self.phangs_hst_ha_obs_band_dict[self.hst_ha_target_name]['ha_observed']:
                hst_habu_band_red = 'F658N'
                band_list += ['F658N']
                hst_ha_stamp_band = 'F658N'
            else:
                hst_ha_stamp_band = None
                hst_habu_band_red = None
            hst_habu_band_green = hst_bvi_band_blue
            hst_habu_band_blue = 'F336W'
        else:
            hst_ha_stamp_band = None
            hst_habu_band_red = None
            hst_habu_band_green = None
            hst_habu_band_blue = None
        hst_stamp_band_list = self.sort_band_list(hst_stamp_band_list)
        if include_nircam:
            band_list += self.get_nircam_band_list()
            nircam_stamp_band_list = self.get_nircam_band_list()
            nircam_band_red = 'F300M'
            nircam_band_green = 'F335M'
            nircam_band_blue = 'F200W'
        else:
            nircam_stamp_band_list = []
            nircam_band_red = None
            nircam_band_green = None
            nircam_band_blue = None

        if include_miri:
            band_list += self.get_miri_band_list()
            miri_stamp_band_list = self.get_miri_band_list()
            miri_band_red = 'F1130W'
            miri_band_green = 'F1000W'
            miri_band_blue = 'F770W'
        else:
            miri_stamp_band_list = []
            miri_band_red = None
            miri_band_green = None
            miri_band_blue = None

        print(band_list)

        # load all bands into constructor
        self.load_phangs_bands(band_list=band_list, flux_unit='MJy/sr', load_err=False)

        # get hst_bvi_zoom_in
        img_hst_bvi_overview, wcs_hst_bvi_overview = self.get_target_overview_rgb_img(red_band=hst_bvi_band_red,
                                                                                          green_band=hst_bvi_band_green,
                                                                                          blue_band=hst_bvi_band_blue,
                                                                                          overview_img_size=(500, 500))
        img_hst_bvi_zoom_in, wcs_hst_bvi_zoom_in = self.get_rgb_zoom_in(ra=ra_region, dec=dec_region,
                                                                        cutout_size=env_cutout_size,
                                                                        circle_rad=circle_rad_region,
                                                                        band_red=hst_bvi_band_red,
                                                                        band_green=hst_bvi_band_green,
                                                                        band_blue=hst_bvi_band_blue)

        if include_h_alpha:
            img_hst_habu_zoom_in, wcs_hst_habu_zoom_in = self.get_rgb_zoom_in(ra=ra_region, dec=dec_region,
                                                                              cutout_size=env_cutout_size,
                                                                              circle_rad=circle_rad_region,
                                                                              band_red=hst_habu_band_red,
                                                                              band_green=hst_habu_band_green,
                                                                              band_blue=hst_habu_band_blue)

        else:
            img_hst_habu_zoom_in, wcs_hst_habu_zoom_in = None, None
        if include_nircam:
            img_nircam_zoom_in, wcs_nircam_zoom_in = self.get_rgb_zoom_in(ra=ra_region, dec=dec_region,
                                                                          cutout_size=env_cutout_size,
                                                                          circle_rad=circle_rad_region,
                                                                          band_red=nircam_band_red,
                                                                          band_green=nircam_band_green,
                                                                          band_blue=nircam_band_blue)
        else:
            img_nircam_zoom_in, wcs_nircam_zoom_in = None, None
        if include_miri:
            img_miri_zoom_in, wcs_miri_zoom_in = self.get_rgb_zoom_in(ra=ra_region, dec=dec_region,
                                                                      cutout_size=env_cutout_size,
                                                                      circle_rad=circle_rad_region,
                                                                      band_red=miri_band_red,
                                                                      band_green=miri_band_green,
                                                                      band_blue=miri_band_blue)
        else:
            img_miri_zoom_in, wcs_miri_zoom_in = None, None

        # load cutout stamps
        cutout_dict_stamp = self.get_band_cutout_dict(ra_cutout=ra_region, dec_cutout=dec_region,
                                                      cutout_size=stamp_cutout_size,
                                                      band_list=band_list)

        # plotting
        figure = plt.figure(figsize=(30, 25))


        # limits for the overview image
        overview_width = 0.51
        overview_height = 0.51
        overview_left_align = -0.0
        overview_bottom_align = 0.44

        # limits for the zoom in panels
        zoom_in_width = 0.29
        zoom_in_height = 0.29
        zoom_in_left_align = 0.47
        zoom_in_bottom_align = 0.40
        zoom_in_space_horizontal = 0.01
        zoom_in_space_vertical = -0.04

        # limits for the stamps
        stamp_width = 0.12
        stamp_height = 0.12
        stamp_left_align = 0.01
        stamp_bottom_align = 0.2
        stamp_space_horizontal = 0.05
        stamp_space_vertical = -0.017


        ax_hst_bvi_overview = figure.add_axes([overview_left_align, overview_bottom_align,
                                               overview_width, overview_height],
                                              projection=wcs_hst_bvi_overview)

        ax_hst_bvi_zoom_in = figure.add_axes([zoom_in_left_align,
                                              zoom_in_bottom_align + zoom_in_height + zoom_in_space_horizontal,
                                              zoom_in_width, zoom_in_height], projection=wcs_hst_bvi_zoom_in)

        if include_h_alpha:
            ax_hst_habu_zoom_in = figure.add_axes([zoom_in_left_align + zoom_in_width + zoom_in_space_vertical,
                                                   zoom_in_bottom_align + zoom_in_height + zoom_in_space_horizontal,
                                                   zoom_in_width, zoom_in_height], projection=wcs_hst_habu_zoom_in)
        else:
            ax_hst_habu_zoom_in = None
        if include_nircam:
            ax_nircam_zoom_in = figure.add_axes([zoom_in_left_align,
                                                 zoom_in_bottom_align, zoom_in_width, zoom_in_height],
                                                projection=wcs_nircam_zoom_in)
        else:
            ax_nircam_zoom_in = None
        if include_miri:
            ax_miri_zoom_in = figure.add_axes([zoom_in_left_align + zoom_in_width + zoom_in_space_vertical, zoom_in_bottom_align,
                                               zoom_in_width, zoom_in_height], projection=wcs_miri_zoom_in)
        else:
            ax_miri_zoom_in = None

        # plot the RGB images
        ax_hst_bvi_overview.imshow(img_hst_bvi_overview)
        pe = [patheffects.withStroke(linewidth=3, foreground="w")]
        ax_hst_bvi_overview.text(0.02, 0.98, 'HST', horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='white',
                                 transform=ax_hst_bvi_overview.transAxes, path_effects=pe)
        ax_hst_bvi_overview.text(0.02, 0.95, hst_bvi_band_red.upper(), horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='red',
                                 transform=ax_hst_bvi_overview.transAxes, path_effects=pe)
        ax_hst_bvi_overview.text(0.02, 0.92, hst_bvi_band_green.upper(), horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='green',
                                 transform=ax_hst_bvi_overview.transAxes, path_effects=pe)
        ax_hst_bvi_overview.text(0.02, 0.89, hst_bvi_band_blue.upper(), horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='blue',
                                 transform=ax_hst_bvi_overview.transAxes, path_effects=pe)
        ax_hst_bvi_overview.set_title(self.target_name.upper(), fontsize=fontsize_large + 10)
        self.arr_axis_params(ax=ax_hst_bvi_overview, ra_tick_label=True, dec_tick_label=True,
                        ra_axis_label='R.A. (2000.0)', dec_axis_label='DEC. (2000.0)',
                        ra_minpad=0.8, dec_minpad=0.8, tick_color='white', label_color='k',
                        fontsize=fontsize_large, labelsize=fontsize_large, ra_tick_num=3, dec_tick_num=3)
        helper_func.draw_box(ax=ax_hst_bvi_overview, wcs=wcs_hst_bvi_overview,
                             coord=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                             box_size=env_cutout_size, color='red', line_width=2, line_style='--')
        helper_func.plot_coord_circle(ax=ax_hst_bvi_zoom_in, pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                      rad=circle_rad_region, color='white', line_style='-', line_width=2, alpha=1., fill=False)

        ax_hst_bvi_zoom_in.imshow(img_hst_bvi_zoom_in)
        ax_hst_bvi_zoom_in.text(0.02, 0.98, 'HST', horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='white',
                                 transform=ax_hst_bvi_zoom_in.transAxes, path_effects=pe)
        ax_hst_bvi_zoom_in.text(0.02, 0.92, hst_bvi_band_red.upper(), horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='red',
                                 transform=ax_hst_bvi_zoom_in.transAxes, path_effects=pe)
        ax_hst_bvi_zoom_in.text(0.02, 0.86, hst_bvi_band_green.upper(), horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='green',
                                 transform=ax_hst_bvi_zoom_in.transAxes, path_effects=pe)
        ax_hst_bvi_zoom_in.text(0.02, 0.80, hst_bvi_band_blue.upper(), horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='blue',
                                 transform=ax_hst_bvi_zoom_in.transAxes, path_effects=pe)
        helper_func.draw_box(ax=ax_hst_bvi_zoom_in, wcs=wcs_hst_bvi_zoom_in,
                             coord=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                             box_size=stamp_cutout_size, color='red', line_width=2, line_style='--')
        self.arr_axis_params(ax=ax_hst_bvi_zoom_in, ra_tick_label=False, dec_tick_label=True,
                        ra_axis_label=' ', dec_axis_label=' ',
                        ra_minpad=0.8, dec_minpad=0.8, tick_color='white', label_color='k',
                        fontsize=fontsize_large, labelsize=fontsize_large, ra_tick_num=3, dec_tick_num=3)
        if include_h_alpha:
            ax_hst_habu_zoom_in.imshow(img_hst_habu_zoom_in)
            self.arr_axis_params(ax=ax_hst_habu_zoom_in, ra_tick_label=False, dec_tick_label=False,
                            ra_axis_label=' ', dec_axis_label=' ',
                            ra_minpad=0.8, dec_minpad=0.8, tick_color='white', label_color='k',
                            fontsize=fontsize_large, labelsize=fontsize_large, ra_tick_num=3, dec_tick_num=3)
            ax_hst_habu_zoom_in.imshow(img_hst_habu_zoom_in)
            ax_hst_habu_zoom_in.text(0.02, 0.98, 'HST', horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='white',
                                     transform=ax_hst_habu_zoom_in.transAxes, path_effects=pe)
            ax_hst_habu_zoom_in.text(0.02, 0.92, hst_habu_band_red.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='red',
                                     transform=ax_hst_habu_zoom_in.transAxes, path_effects=pe)
            ax_hst_habu_zoom_in.text(0.02, 0.86, hst_habu_band_green.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='green',
                                     transform=ax_hst_habu_zoom_in.transAxes, path_effects=pe)
            ax_hst_habu_zoom_in.text(0.02, 0.80, hst_habu_band_blue.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='blue',
                                     transform=ax_hst_habu_zoom_in.transAxes, path_effects=pe)
            helper_func.draw_box(ax=ax_hst_habu_zoom_in, wcs=wcs_hst_habu_zoom_in,
                                 coord=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                 box_size=stamp_cutout_size, color='red', line_width=2, line_style='--')
            helper_func.plot_coord_circle(ax=ax_hst_habu_zoom_in, pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                          rad=circle_rad_region, color='white', line_style='-', line_width=2, alpha=1., fill=False)
        if include_nircam:
            ax_nircam_zoom_in.imshow(img_nircam_zoom_in)
            self.arr_axis_params(ax=ax_nircam_zoom_in, ra_tick_label=True, dec_tick_label=True,
                            ra_axis_label=' ', dec_axis_label=' ',
                            ra_minpad=0.8, dec_minpad=0.8, tick_color='white', label_color='k',
                            fontsize=fontsize_large, labelsize=fontsize_large, ra_tick_num=3, dec_tick_num=3)
            ax_nircam_zoom_in.text(0.02, 0.98, 'NIRCAM', horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='white',
                                     transform=ax_nircam_zoom_in.transAxes, path_effects=pe)
            ax_nircam_zoom_in.text(0.02, 0.92, nircam_band_red.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='red',
                                     transform=ax_nircam_zoom_in.transAxes, path_effects=pe)
            ax_nircam_zoom_in.text(0.02, 0.86, nircam_band_green.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='green',
                                     transform=ax_nircam_zoom_in.transAxes, path_effects=pe)
            ax_nircam_zoom_in.text(0.02, 0.80, nircam_band_blue.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='blue',
                                     transform=ax_nircam_zoom_in.transAxes, path_effects=pe)
            helper_func.draw_box(ax=ax_nircam_zoom_in, wcs=wcs_nircam_zoom_in,
                                 coord=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                 box_size=stamp_cutout_size, color='red', line_width=2, line_style='--')
            helper_func.plot_coord_circle(ax=ax_nircam_zoom_in, pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                          rad=circle_rad_region, color='white', line_style='-', line_width=2, alpha=1., fill=False)
        if include_miri:
            ax_miri_zoom_in.imshow(img_miri_zoom_in)
            self.arr_axis_params(ax=ax_miri_zoom_in, ra_tick_label=True, dec_tick_label=False,
                            ra_axis_label=' ', dec_axis_label=' ',
                            ra_minpad=0.8, dec_minpad=0.8, tick_color='white', label_color='k',
                            fontsize=fontsize_large, labelsize=fontsize_large, ra_tick_num=3, dec_tick_num=3)
            ax_miri_zoom_in.text(0.02, 0.98, 'MIRI', horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='white',
                                     transform=ax_miri_zoom_in.transAxes, path_effects=pe)
            ax_miri_zoom_in.text(0.02, 0.92, miri_band_red.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='red',
                                     transform=ax_miri_zoom_in.transAxes, path_effects=pe)
            ax_miri_zoom_in.text(0.02, 0.86, miri_band_green.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='green',
                                     transform=ax_miri_zoom_in.transAxes, path_effects=pe)
            ax_miri_zoom_in.text(0.02, 0.80, miri_band_blue.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='blue',
                                     transform=ax_miri_zoom_in.transAxes, path_effects=pe)
            helper_func.draw_box(ax=ax_miri_zoom_in, wcs=wcs_miri_zoom_in,
                                 coord=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                 box_size=stamp_cutout_size, color='red', line_width=2, line_style='--')
            helper_func.plot_coord_circle(ax=ax_miri_zoom_in, pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                          rad=circle_rad_region, color='white', line_style='-', line_width=2, alpha=1., fill=False)

        # add stamp axis
        stamp_row_count_down = 0
        # add stamp axis for hst
        ax_hst_stamp_list = []
        hst_stamp_index = 0

        for hst_stamp_index, hst_stamp_band in enumerate(hst_stamp_band_list):

            ax_hst_stamp_list.append(
                figure.add_axes([stamp_left_align + hst_stamp_index*(stamp_width + stamp_space_vertical),
                                 stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                 stamp_width, stamp_height],
                                projection=cutout_dict_stamp['%s_img_cutout' % hst_stamp_band].wcs))
            norm_hst_stamp = helper_func.compute_cbar_norm(
                cutout_list=cutout_dict_stamp['%s_img_cutout' % hst_stamp_band].data, log_scale=cbar_log_scale)
            ax_hst_stamp_list[hst_stamp_index].imshow(cutout_dict_stamp['%s_img_cutout' % hst_stamp_band].data,
                                            norm=norm_hst_stamp, cmap='Greys')
            ax_hst_stamp_list[hst_stamp_index].set_title(hst_stamp_band.upper(), fontsize=fontsize_large)
            if hst_stamp_index == 0:
                ra_tick_label, dec_tick_label = True, True
            else:
                ra_tick_label, dec_tick_label = False, False
            self.arr_axis_params(ax=ax_hst_stamp_list[hst_stamp_index], ra_tick_label=ra_tick_label,
                                 dec_tick_label=dec_tick_label, ra_axis_label=' ', dec_axis_label=' ',
                                 fontsize=fontsize_small, labelsize=fontsize_small)

        if include_h_alpha:
            ax_hst_ha_stamp = figure.add_axes([stamp_left_align + (hst_stamp_index + 2)*(stamp_width + stamp_space_vertical),
                                 stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                 stamp_width, stamp_height],
                                projection=cutout_dict_stamp['%s_img_cutout' % hst_ha_stamp_band].wcs)
            norm_hst_ha_stamp = helper_func.compute_cbar_norm(
                cutout_list=cutout_dict_stamp['%s_img_cutout' % hst_ha_stamp_band].data, log_scale=cbar_log_scale)
            ax_hst_ha_stamp.imshow(cutout_dict_stamp['%s_img_cutout' % hst_ha_stamp_band].data,
                                   norm=norm_hst_ha_stamp, cmap='Greys')
            ax_hst_ha_stamp.set_title(hst_ha_stamp_band.upper(), fontsize=fontsize_large)
            self.arr_axis_params(ax=ax_hst_ha_stamp, ra_tick_label=True, dec_tick_label=True,
                        ra_axis_label=' ', dec_axis_label=' ',
                        ra_minpad=0.8, dec_minpad=0.8, tick_color='k', label_color='k',
                        fontsize=fontsize_small, labelsize=fontsize_small, ra_tick_num=3, dec_tick_num=3)

            ax_hst_ha_cont_sub_stamp = figure.add_axes([stamp_left_align + (hst_stamp_index + 4)*(stamp_width + stamp_space_vertical),
                                 stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                 stamp_width, stamp_height],
                                projection=cutout_dict_stamp['%s_img_cutout' % ha_cont_sub_stamp_band].wcs)
            norm_hst_ha_cont_sub_stamp = helper_func.compute_cbar_norm(
                cutout_list=cutout_dict_stamp['%s_img_cutout' % ha_cont_sub_stamp_band].data, log_scale=cbar_log_scale)
            ax_hst_ha_cont_sub_stamp.imshow(cutout_dict_stamp['%s_img_cutout' % ha_cont_sub_stamp_band].data,
                                   norm=norm_hst_ha_cont_sub_stamp, cmap='Greys')
            ax_hst_ha_cont_sub_stamp.set_title(ha_cont_sub_stamp_band.upper(), fontsize=fontsize_large)
            self.arr_axis_params(ax=ax_hst_ha_cont_sub_stamp, ra_tick_label=True, dec_tick_label=True,
                        ra_axis_label=' ', dec_axis_label=' ',
                        ra_minpad=0.8, dec_minpad=0.8, tick_color='k', label_color='k',
                        fontsize=fontsize_small, labelsize=fontsize_small, ra_tick_num=3, dec_tick_num=3)

        nircam_stamp_index = 0
        if include_nircam:
            stamp_row_count_down -= 1
            ax_nircam_stamp_list = []
            for nircam_stamp_index, nircam_stamp_band in enumerate(nircam_stamp_band_list):
                if np.sum(cutout_dict_stamp['%s_img_cutout' % nircam_stamp_band].data) == 0:
                    continue
                ax_nircam_stamp_list.append(
                    figure.add_axes([stamp_left_align + nircam_stamp_index*(stamp_width + stamp_space_vertical),
                                     stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                     stamp_width, stamp_height],
                                    projection=cutout_dict_stamp['%s_img_cutout' % nircam_stamp_band].wcs))
                norm_nircam_stamp = helper_func.compute_cbar_norm(
                cutout_list=cutout_dict_stamp['%s_img_cutout' % nircam_stamp_band].data, log_scale=cbar_log_scale)
                ax_nircam_stamp_list[nircam_stamp_index].imshow(cutout_dict_stamp['%s_img_cutout' % nircam_stamp_band].data,
                                                norm=norm_nircam_stamp, cmap='Greys')
                ax_nircam_stamp_list[nircam_stamp_index].set_title(nircam_stamp_band.upper(), fontsize=fontsize_large)
                if nircam_stamp_index == 0:
                    ra_tick_label, dec_tick_label = True, True
                else:
                    ra_tick_label, dec_tick_label = False, False
                self.arr_axis_params(ax=ax_nircam_stamp_list[nircam_stamp_index], ra_tick_label=ra_tick_label,
                                     dec_tick_label=dec_tick_label, ra_axis_label=' ', dec_axis_label=' ',
                                     fontsize=fontsize_small, labelsize=fontsize_small)

        if include_miri:
            # stamp_row_count_down -= 1
            ax_miri_stamp_list = []
            for miri_stamp_index, miri_stamp_band in enumerate(miri_stamp_band_list):
                if np.sum(cutout_dict_stamp['%s_img_cutout' % miri_stamp_band].data) == 0:
                    continue
                ax_miri_stamp_list.append(
                    figure.add_axes([stamp_left_align + (nircam_stamp_index + miri_stamp_index + 2)*(stamp_width + stamp_space_vertical),
                                     stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                     stamp_width, stamp_height],
                                    projection=cutout_dict_stamp['%s_img_cutout' % miri_stamp_band].wcs))
                norm_miri_stamp = helper_func.compute_cbar_norm(
                cutout_list=cutout_dict_stamp['%s_img_cutout' % miri_stamp_band].data, log_scale=cbar_log_scale)
                ax_miri_stamp_list[miri_stamp_index].imshow(cutout_dict_stamp['%s_img_cutout' % miri_stamp_band].data,
                                                norm=norm_miri_stamp, cmap='Greys')
                ax_miri_stamp_list[miri_stamp_index].set_title(miri_stamp_band.upper(), fontsize=fontsize_large)
                if miri_stamp_index == 0:
                    ra_tick_label, dec_tick_label = True, True
                else:
                    ra_tick_label, dec_tick_label = False, False
                self.arr_axis_params(ax=ax_miri_stamp_list[miri_stamp_index], ra_tick_label=ra_tick_label,
                                     dec_tick_label=dec_tick_label, ra_axis_label=' ', dec_axis_label=' ',
                                     fontsize=fontsize_small, labelsize=fontsize_small)

        return figure

    def plot_over_view_zoom_in_panel(self,
                                     figure,
                                     ra_region, dec_region,
                                     hst_overview_dict, hst_zoom_in_dict, stamp_dict, hst_habu_zoom_in_dict=None,
                                     nircam_zoom_in_dict=None, miri_zoom_in_dict=None,
                                     radius_dict=None,
                                     env_cutout_size=(10, 10), circle_rad_region=None,
                                     stamp_cutout_size=(2.5, 2.5),
                                     fontsize_large=30, fontsize_small=20,):
        # plot the RGB images
        hst_overview_dict['ax_hst_bvi_overview'].imshow(hst_overview_dict['img_hst_bvi_overview'])
        pe = [patheffects.withStroke(linewidth=3, foreground="w")]
        hst_overview_dict['ax_hst_bvi_overview'].text(0.02, 0.98, 'HST', horizontalalignment='left',
                                                      verticalalignment='top', fontsize=fontsize_large, color='white',
                                                      transform=hst_overview_dict['ax_hst_bvi_overview'].transAxes,
                                                      path_effects=pe)
        hst_overview_dict['ax_hst_bvi_overview'].text(0.02, 0.95, hst_overview_dict['hst_bvi_band_red'].upper(),
                                                      horizontalalignment='left', verticalalignment='top',
                                                      fontsize=fontsize_large, color='red',
                                                      transform=hst_overview_dict['ax_hst_bvi_overview'].transAxes,
                                                      path_effects=pe)
        hst_overview_dict['ax_hst_bvi_overview'].text(0.02, 0.92, hst_overview_dict['hst_bvi_band_green'].upper(),
                                                      horizontalalignment='left', verticalalignment='top',
                                                      fontsize=fontsize_large, color='green',
                                                      transform=hst_overview_dict['ax_hst_bvi_overview'].transAxes,
                                                      path_effects=pe)
        hst_overview_dict['ax_hst_bvi_overview'].text(0.02, 0.89, hst_overview_dict['hst_bvi_band_blue'].upper(),
                                                      horizontalalignment='left', verticalalignment='top',
                                                      fontsize=fontsize_large, color='blue',
                                                      transform=hst_overview_dict['ax_hst_bvi_overview'].transAxes,
                                                      path_effects=pe)
        hst_overview_dict['ax_hst_bvi_overview'].set_title(self.target_name.upper(), fontsize=fontsize_large + 10)
        self.arr_axis_params(ax=hst_overview_dict['ax_hst_bvi_overview'],  tick_color='white',
                             fontsize=fontsize_large, labelsize=fontsize_large)
        helper_func.draw_box(ax=hst_overview_dict['ax_hst_bvi_overview'], wcs=hst_overview_dict['wcs_hst_bvi_overview'],
                             coord=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                             box_size=env_cutout_size, color='red', line_style='--')

        hst_zoom_in_dict['ax_hst_bvi_zoom_in'].imshow(hst_zoom_in_dict['img_hst_bvi_zoom_in'])
        hst_zoom_in_dict['ax_hst_bvi_zoom_in'].text(0.02, 0.98, 'HST', horizontalalignment='left',
                                                    verticalalignment='top', fontsize=fontsize_large, color='white',
                                                    transform=hst_zoom_in_dict['ax_hst_bvi_zoom_in'].transAxes,
                                                    path_effects=pe)
        hst_zoom_in_dict['ax_hst_bvi_zoom_in'].text(0.02, 0.92, hst_zoom_in_dict['hst_bvi_band_red'].upper(),
                                                    horizontalalignment='left', verticalalignment='top',
                                                    fontsize=fontsize_large, color='red',
                                                    transform=hst_zoom_in_dict['ax_hst_bvi_zoom_in'].transAxes,
                                                    path_effects=pe)
        hst_zoom_in_dict['ax_hst_bvi_zoom_in'].text(0.02, 0.86, hst_zoom_in_dict['hst_bvi_band_green'].upper(),
                                                    horizontalalignment='left', verticalalignment='top',
                                                    fontsize=fontsize_large, color='green',
                                                    transform=hst_zoom_in_dict['ax_hst_bvi_zoom_in'].transAxes,
                                                    path_effects=pe)
        hst_zoom_in_dict['ax_hst_bvi_zoom_in'].text(0.02, 0.80, hst_zoom_in_dict['hst_bvi_band_blue'].upper(),
                                                    horizontalalignment='left', verticalalignment='top',
                                                    fontsize=fontsize_large, color='blue',
                                                    transform=hst_zoom_in_dict['ax_hst_bvi_zoom_in'].transAxes,
                                                    path_effects=pe)
        helper_func.draw_box(ax=hst_zoom_in_dict['ax_hst_bvi_zoom_in'], wcs=hst_zoom_in_dict['wcs_hst_bvi_zoom_in'],
                             coord=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                             box_size=stamp_cutout_size, color='cyan', line_style='--')
        if circle_rad_region is not None:
            helper_func.plot_coord_circle(ax=hst_zoom_in_dict['ax_hst_bvi_zoom_in'], pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                          rad=circle_rad_region, color='white', line_style='-', line_width=2, alpha=1., fill=False)
        self.arr_axis_params(ax=hst_zoom_in_dict['ax_hst_bvi_zoom_in'], ra_tick_label=False,  ra_axis_label=' ',
                             dec_axis_label=' ', tick_color='white', fontsize=fontsize_large, labelsize=fontsize_large)
        if hst_habu_zoom_in_dict is not None:
            hst_habu_zoom_in_dict['ax_hst_habu_zoom_in'].imshow(hst_habu_zoom_in_dict['img_hst_habu_zoom_in'])
            self.arr_axis_params(ax=hst_habu_zoom_in_dict['ax_hst_habu_zoom_in'], ra_tick_label=False,
                                 dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ', tick_color='white',
                                 fontsize=fontsize_large, labelsize=fontsize_large)
            hst_habu_zoom_in_dict['ax_hst_habu_zoom_in'].imshow(hst_habu_zoom_in_dict['img_hst_habu_zoom_in'])
            hst_habu_zoom_in_dict['ax_hst_habu_zoom_in'].text(0.02, 0.98, 'HST', horizontalalignment='left',
                                                              verticalalignment='top', fontsize=fontsize_large,
                                                              color='white',
                                                              transform=
                                                              hst_habu_zoom_in_dict['ax_hst_habu_zoom_in'].transAxes,
                                                              path_effects=pe)
            hst_habu_zoom_in_dict['ax_hst_habu_zoom_in'].text(0.02, 0.92,
                                                              hst_habu_zoom_in_dict['hst_habu_band_red'].upper(),
                                                              horizontalalignment='left', verticalalignment='top',
                                                              fontsize=fontsize_large, color='red',
                                                              transform=
                                                              hst_habu_zoom_in_dict['ax_hst_habu_zoom_in'].transAxes,
                                                              path_effects=pe)
            hst_habu_zoom_in_dict['ax_hst_habu_zoom_in'].text(0.02, 0.86,
                                                              hst_habu_zoom_in_dict['hst_habu_band_green'].upper(),
                                                              horizontalalignment='left', verticalalignment='top',
                                                              fontsize=fontsize_large, color='green',
                                                              transform=
                                                              hst_habu_zoom_in_dict['ax_hst_habu_zoom_in'].transAxes,
                                                              path_effects=pe)
            hst_habu_zoom_in_dict['ax_hst_habu_zoom_in'].text(0.02, 0.80,
                                                              hst_habu_zoom_in_dict['hst_habu_band_blue'].upper(),
                                                              horizontalalignment='left', verticalalignment='top',
                                                              fontsize=fontsize_large, color='blue',
                                                              transform=
                                                              hst_habu_zoom_in_dict['ax_hst_habu_zoom_in'].transAxes,
                                                              path_effects=pe)
            helper_func.draw_box(ax=hst_habu_zoom_in_dict['ax_hst_habu_zoom_in'],
                                 wcs=hst_habu_zoom_in_dict['wcs_hst_habu_zoom_in'],
                                 coord=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                 box_size=stamp_cutout_size, color='cyan', line_style='--')
            if circle_rad_region is not None:
                helper_func.plot_coord_circle(ax=hst_habu_zoom_in_dict['ax_hst_habu_zoom_in'],
                                              pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                              rad=circle_rad_region, color='white', line_width=2)
        if nircam_zoom_in_dict is not None:
            nircam_zoom_in_dict['ax_nircam_zoom_in'].imshow(nircam_zoom_in_dict['img_nircam_zoom_in'])
            self.arr_axis_params(ax=nircam_zoom_in_dict['ax_nircam_zoom_in'],  ra_axis_label=' ', dec_axis_label=' ',
                                 tick_color='white', fontsize=fontsize_large, labelsize=fontsize_large)
            nircam_zoom_in_dict['ax_nircam_zoom_in'].text(0.02, 0.98, 'NIRCAM', horizontalalignment='left',
                                                          verticalalignment='top', fontsize=fontsize_large,
                                                          color='white', transform=
                                                          nircam_zoom_in_dict['ax_nircam_zoom_in'].transAxes,
                                                          path_effects=pe)
            nircam_zoom_in_dict['ax_nircam_zoom_in'].text(0.02, 0.92, nircam_zoom_in_dict['nircam_band_red'].upper(),
                                                          horizontalalignment='left', verticalalignment='top',
                                                          fontsize=fontsize_large, color='red', transform=
                                                          nircam_zoom_in_dict['ax_nircam_zoom_in'].transAxes,
                                                          path_effects=pe)
            nircam_zoom_in_dict['ax_nircam_zoom_in'].text(0.02, 0.86, nircam_zoom_in_dict['nircam_band_green'].upper(),
                                                          horizontalalignment='left', verticalalignment='top',
                                                          fontsize=fontsize_large, color='green', transform=
                                                          nircam_zoom_in_dict['ax_nircam_zoom_in'].transAxes,
                                                          path_effects=pe)
            nircam_zoom_in_dict['ax_nircam_zoom_in'].text(0.02, 0.80, nircam_zoom_in_dict['nircam_band_blue'].upper(),
                                                          horizontalalignment='left', verticalalignment='top',
                                                          fontsize=fontsize_large, color='blue', transform=
                                                          nircam_zoom_in_dict['ax_nircam_zoom_in'].transAxes,
                                                          path_effects=pe)
            helper_func.draw_box(ax=nircam_zoom_in_dict['ax_nircam_zoom_in'], wcs=nircam_zoom_in_dict['wcs_nircam_zoom_in'],
                                 coord=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg), box_size=stamp_cutout_size,
                                 color='cyan', line_style='--')
            if circle_rad_region is not None:
                helper_func.plot_coord_circle(ax=nircam_zoom_in_dict['ax_nircam_zoom_in'],
                                              pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                              rad=circle_rad_region, color='white', line_width=2)
        if miri_zoom_in_dict is not None:
            miri_zoom_in_dict['ax_miri_zoom_in'].imshow(miri_zoom_in_dict['img_miri_zoom_in'])
            self.arr_axis_params(ax=miri_zoom_in_dict['ax_miri_zoom_in'], dec_tick_label=False, ra_axis_label=' ',
                                 dec_axis_label=' ', tick_color='white',  fontsize=fontsize_large,
                                 labelsize=fontsize_large)
            miri_zoom_in_dict['ax_miri_zoom_in'].text(0.02, 0.98, 'MIRI', horizontalalignment='left',
                                                      verticalalignment='top', fontsize=fontsize_large, color='white',
                                                      transform=miri_zoom_in_dict['ax_miri_zoom_in'].transAxes,
                                                      path_effects=pe)
            miri_zoom_in_dict['ax_miri_zoom_in'].text(0.02, 0.92, miri_zoom_in_dict['miri_band_red'].upper(),
                                                      horizontalalignment='left', verticalalignment='top',
                                                      fontsize=fontsize_large, color='red',
                                                      transform=miri_zoom_in_dict['ax_miri_zoom_in'].transAxes,
                                                      path_effects=pe)
            miri_zoom_in_dict['ax_miri_zoom_in'].text(0.02, 0.86, miri_zoom_in_dict['miri_band_green'].upper(),
                                                      horizontalalignment='left', verticalalignment='top',
                                                      fontsize=fontsize_large, color='green',
                                                      transform=miri_zoom_in_dict['ax_miri_zoom_in'].transAxes,
                                                      path_effects=pe)
            miri_zoom_in_dict['ax_miri_zoom_in'].text(0.02, 0.80, miri_zoom_in_dict['miri_band_blue'].upper(),
                                                      horizontalalignment='left', verticalalignment='top',
                                                      fontsize=fontsize_large, color='blue',
                                                      transform=miri_zoom_in_dict['ax_miri_zoom_in'].transAxes,
                                                      path_effects=pe)
            helper_func.draw_box(ax=miri_zoom_in_dict['ax_miri_zoom_in'], wcs=miri_zoom_in_dict['wcs_miri_zoom_in'],
                                 coord=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                 box_size=stamp_cutout_size, color='cyan', line_style='--')
            if circle_rad_region is not None:
                helper_func.plot_coord_circle(ax=miri_zoom_in_dict['ax_miri_zoom_in'],
                                              pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                              rad=circle_rad_region, color='white', line_width=2)

        # add stamp axis
        stamp_row_count_down = 0
        # add stamp axis for hst
        ax_hst_stamp_list = []
        hst_stamp_index = 0

        for hst_stamp_index, hst_stamp_band in enumerate(stamp_dict['hst_stamp_band_list']):

            ax_hst_stamp_list.append(
                figure.add_axes([stamp_dict['stamp_left_align'] +
                                 hst_stamp_index*(stamp_dict['stamp_width'] + stamp_dict['stamp_space_vertical']),
                                 stamp_dict['stamp_bottom_align'] +
                                 stamp_row_count_down*(stamp_dict['stamp_height'] +
                                                       stamp_dict['stamp_space_horizontal']),
                                 stamp_dict['stamp_width'], stamp_dict['stamp_height']],
                                projection=stamp_dict['cutout_dict_stamp']['%s_img_cutout' % hst_stamp_band].wcs))
            norm_hst_stamp = helper_func.compute_cbar_norm(
                cutout_list=stamp_dict['cutout_dict_stamp']['%s_img_cutout' % hst_stamp_band].data,
                log_scale=stamp_dict['cbar_log_scale'])
            ax_hst_stamp_list[hst_stamp_index].imshow(
                stamp_dict['cutout_dict_stamp']['%s_img_cutout' % hst_stamp_band].data, norm=norm_hst_stamp,
                cmap='Greys')
            if radius_dict is not None:
                helper_func.plot_coord_circle(ax=ax_hst_stamp_list[hst_stamp_index],
                                              pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                              rad=radius_dict['aperture_%s_ee50' % hst_stamp_band],
                                              color='cyan', line_width=2)
                helper_func.plot_coord_circle(ax=ax_hst_stamp_list[hst_stamp_index],
                                              pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                              rad=radius_dict['aperture_%s_ee80' % hst_stamp_band],
                                              color='red', line_width=2)

            helper_func.plot_coord_croshair(ax=ax_hst_stamp_list[hst_stamp_index],
                                              pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                            wcs=stamp_dict['cutout_dict_stamp']['%s_img_cutout' % hst_stamp_band].wcs,
                                              rad=0.3, hair_length=0.3,
                                              color='red', line_width=2)

            ax_hst_stamp_list[hst_stamp_index].set_title(hst_stamp_band.upper(), fontsize=fontsize_large,
                                                         color='tab:blue')
            if hst_stamp_index == 0:
                ra_tick_label, dec_tick_label = True, True
            else:
                ra_tick_label, dec_tick_label = False, False
            self.arr_axis_params(ax=ax_hst_stamp_list[hst_stamp_index], ra_tick_label=ra_tick_label,
                                 dec_tick_label=dec_tick_label, ra_axis_label=' ', dec_axis_label=' ',
                                 fontsize=fontsize_small, labelsize=fontsize_small)

        if hst_habu_zoom_in_dict is not None:
            ax_hst_ha_stamp = figure.add_axes([stamp_dict['stamp_left_align'] +
                                               (hst_stamp_index + 2)*(stamp_dict['stamp_width'] +
                                                                      stamp_dict['stamp_space_vertical']),
                                               stamp_dict['stamp_bottom_align'] +
                                               stamp_row_count_down*(stamp_dict['stamp_height'] +
                                                                     stamp_dict['stamp_space_horizontal']),
                                               stamp_dict['stamp_width'], stamp_dict['stamp_height']],
                                              projection=
                                              stamp_dict['cutout_dict_stamp']['%s_img_cutout' %
                                                                              stamp_dict['hst_ha_stamp_band']].wcs)
            norm_hst_ha_stamp = helper_func.compute_cbar_norm(
                cutout_list=stamp_dict['cutout_dict_stamp']['%s_img_cutout' % stamp_dict['hst_ha_stamp_band']].data,
                log_scale=stamp_dict['cbar_log_scale'])
            ax_hst_ha_stamp.imshow(
                stamp_dict['cutout_dict_stamp']['%s_img_cutout' % stamp_dict['hst_ha_stamp_band']].data,
                norm=norm_hst_ha_stamp, cmap='Greys')

            if radius_dict is not None:
                helper_func.plot_coord_circle(ax=ax_hst_ha_stamp,
                                              pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                              rad=radius_dict['aperture_%s_ee50' % stamp_dict['hst_ha_stamp_band']],
                                              color='cyan', line_width=2)
                helper_func.plot_coord_circle(ax=ax_hst_ha_stamp,
                                              pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                              rad=radius_dict['aperture_%s_ee80' % stamp_dict['hst_ha_stamp_band']],
                                              color='red', line_width=2)

            helper_func.plot_coord_croshair(ax=ax_hst_ha_stamp,
                                              pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                            wcs=stamp_dict['cutout_dict_stamp']['%s_img_cutout' % stamp_dict['hst_ha_stamp_band']].wcs,
                                              rad=0.3, hair_length=0.3,
                                              color='red', line_width=2)
            ax_hst_ha_stamp.set_title(stamp_dict['hst_ha_stamp_band'].upper(), fontsize=fontsize_large, color='tab:red')
            self.arr_axis_params(ax=ax_hst_ha_stamp, ra_axis_label=' ', dec_axis_label=' ', fontsize=fontsize_small,
                                 labelsize=fontsize_small)

            ax_hst_ha_cont_sub_stamp = figure.add_axes([stamp_dict['stamp_left_align'] +
                                                        (hst_stamp_index + 4)*(stamp_dict['stamp_width'] +
                                                                               stamp_dict['stamp_space_vertical']),
                                                        stamp_dict['stamp_bottom_align'] +
                                                        stamp_row_count_down*(stamp_dict['stamp_height'] +
                                                                              stamp_dict['stamp_space_horizontal']),
                                                        stamp_dict['stamp_width'], stamp_dict['stamp_height']],
                                                       projection=
                                                       stamp_dict['cutout_dict_stamp']
                                                       ['%s_img_cutout' % stamp_dict['ha_cont_sub_stamp_band']].wcs)
            norm_hst_ha_cont_sub_stamp = helper_func.compute_cbar_norm(
                cutout_list=stamp_dict['cutout_dict_stamp']['%s_img_cutout' %
                                                            stamp_dict['ha_cont_sub_stamp_band']].data,
                log_scale=stamp_dict['cbar_log_scale'])
            ax_hst_ha_cont_sub_stamp.imshow(stamp_dict['cutout_dict_stamp']['%s_img_cutout' %
                                                                            stamp_dict['ha_cont_sub_stamp_band']].data,
                                            norm=norm_hst_ha_cont_sub_stamp, cmap='Greys')
            helper_func.plot_coord_croshair(ax=ax_hst_ha_cont_sub_stamp,
                                              pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                            wcs=stamp_dict['cutout_dict_stamp']['%s_img_cutout' % stamp_dict['ha_cont_sub_stamp_band']].wcs,
                                              rad=0.3, hair_length=0.3,
                                              color='red', line_width=2)
            ax_hst_ha_cont_sub_stamp.set_title(r'Cont. sub H$\alpha$', fontsize=fontsize_large,
                                               color='tab:red')
            self.arr_axis_params(ax=ax_hst_ha_cont_sub_stamp, ra_axis_label=' ', dec_axis_label=' ',
                                 fontsize=fontsize_small, labelsize=fontsize_small)

        nircam_stamp_index = 0
        stamp_row_count_down -= 1
        if nircam_zoom_in_dict is not None:
            ax_nircam_stamp_list = []
            for nircam_stamp_index, nircam_stamp_band in enumerate(stamp_dict['nircam_stamp_band_list']):
                if np.sum(stamp_dict['cutout_dict_stamp']['%s_img_cutout' % nircam_stamp_band].data) == 0:
                    continue
                ax_nircam_stamp_list.append(
                    figure.add_axes([stamp_dict['stamp_left_align'] +
                                     nircam_stamp_index*(stamp_dict['stamp_width'] +
                                                         stamp_dict['stamp_space_vertical']),
                                     stamp_dict['stamp_bottom_align'] +
                                     stamp_row_count_down*(stamp_dict['stamp_height'] +
                                                           stamp_dict['stamp_space_horizontal']),
                                     stamp_dict['stamp_width'], stamp_dict['stamp_height']],
                                    projection=stamp_dict['cutout_dict_stamp']['%s_img_cutout' %
                                                                               nircam_stamp_band].wcs))
                norm_nircam_stamp = helper_func.compute_cbar_norm(
                    cutout_list=stamp_dict['cutout_dict_stamp']['%s_img_cutout' % nircam_stamp_band].data,
                    log_scale=stamp_dict['cbar_log_scale'])
                ax_nircam_stamp_list[nircam_stamp_index].imshow(stamp_dict['cutout_dict_stamp']['%s_img_cutout' %
                                                                                                nircam_stamp_band].data,
                                                                norm=norm_nircam_stamp, cmap='Greys')
                if radius_dict is not None:
                    helper_func.plot_coord_circle(ax=ax_nircam_stamp_list[nircam_stamp_index],
                                                  pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                                  rad=radius_dict['aperture_%s_ee50' % nircam_stamp_band],
                                                  color='cyan', line_width=2)
                    helper_func.plot_coord_circle(ax=ax_nircam_stamp_list[nircam_stamp_index],
                                                  pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                                  rad=radius_dict['aperture_%s_ee80' % nircam_stamp_band],
                                                  color='red', line_width=2)
                helper_func.plot_coord_croshair(ax=ax_nircam_stamp_list[nircam_stamp_index],
                                              pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                            wcs=stamp_dict['cutout_dict_stamp']['%s_img_cutout' % nircam_stamp_band].wcs,
                                              rad=0.3, hair_length=0.3,
                                              color='red', line_width=2)
                ax_nircam_stamp_list[nircam_stamp_index].set_title(nircam_stamp_band.upper(), fontsize=fontsize_large,
                                                                   color='tab:green')
                if nircam_stamp_index == 0:
                    ra_tick_label, dec_tick_label = True, True
                else:
                    ra_tick_label, dec_tick_label = False, False
                self.arr_axis_params(ax=ax_nircam_stamp_list[nircam_stamp_index], ra_tick_label=ra_tick_label,
                                     dec_tick_label=dec_tick_label, ra_axis_label=' ', dec_axis_label=' ',
                                     fontsize=fontsize_small, labelsize=fontsize_small)

        if miri_zoom_in_dict is not None:
            # stamp_row_count_down -= 1
            ax_miri_stamp_list = []
            for miri_stamp_index, miri_stamp_band in enumerate(stamp_dict['miri_stamp_band_list']):
                if np.sum(stamp_dict['cutout_dict_stamp']['%s_img_cutout' % miri_stamp_band].data) == 0:
                    continue
                ax_miri_stamp_list.append(
                    figure.add_axes([stamp_dict['stamp_left_align'] +
                                     (nircam_stamp_index + miri_stamp_index + 2)*(stamp_dict['stamp_width'] +
                                                                                  stamp_dict['stamp_space_vertical']),
                                     stamp_dict['stamp_bottom_align'] +
                                     stamp_row_count_down*(stamp_dict['stamp_height'] +
                                                           stamp_dict['stamp_space_horizontal']),
                                     stamp_dict['stamp_width'], stamp_dict['stamp_height']],
                                    projection=stamp_dict['cutout_dict_stamp']['%s_img_cutout' % miri_stamp_band].wcs))
                norm_miri_stamp = helper_func.compute_cbar_norm(
                    cutout_list=stamp_dict['cutout_dict_stamp']['%s_img_cutout' % miri_stamp_band].data,
                    log_scale=stamp_dict['cbar_log_scale'])
                ax_miri_stamp_list[miri_stamp_index].imshow(
                    stamp_dict['cutout_dict_stamp']['%s_img_cutout' % miri_stamp_band].data,
                    norm=norm_miri_stamp, cmap='Greys')
                if radius_dict is not None:
                    helper_func.plot_coord_circle(ax=ax_miri_stamp_list[miri_stamp_index],
                                                  pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                                  rad=radius_dict['aperture_%s_ee50' % miri_stamp_band],
                                                  color='cyan', line_width=2)
                    helper_func.plot_coord_circle(ax=ax_miri_stamp_list[miri_stamp_index],
                                                  pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                                  rad=radius_dict['aperture_%s_ee80' % miri_stamp_band],
                                                  color='red', line_width=2)
                helper_func.plot_coord_croshair(ax=ax_miri_stamp_list[miri_stamp_index],
                                              pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                            wcs=stamp_dict['cutout_dict_stamp']['%s_img_cutout' % miri_stamp_band].wcs,
                                              rad=0.3, hair_length=0.3,
                                              color='red', line_width=2)
                ax_miri_stamp_list[miri_stamp_index].set_title(miri_stamp_band.upper(), fontsize=fontsize_large,
                                                               color='tab:purple')
                if miri_stamp_index == 0:
                    ra_tick_label, dec_tick_label = True, True
                else:
                    ra_tick_label, dec_tick_label = False, False
                self.arr_axis_params(ax=ax_miri_stamp_list[miri_stamp_index], ra_tick_label=ra_tick_label,
                                     dec_tick_label=dec_tick_label, ra_axis_label=' ', dec_axis_label=' ',
                                     fontsize=fontsize_small, labelsize=fontsize_small)


    def plot_overview_zoom_in_co_env_panel(self, ra_region, dec_region, mol_gas_dens_dict,
                                           string_dict,  color_color_dict,
                                           include_h_alpha=True, ha_cont_sub_stamp_band='ha',
                                           include_nircam=True,
                                           include_miri=True,
                                           env_cutout_size=(10, 10), circle_rad_region=1.25, circle_rad_obj=0.25,
                                           stamp_cutout_size=(2.5, 2.5), cbar_log_scale=True,
                                           fontsize_large=30, fontsize_small=20,
                                           ):

        # load all bands needed
        band_list = []
        # load HST bands
        band_list += self.get_hst_band_list()
        hst_stamp_band_list = band_list.copy()
        # get BVI filters
        # specify color the bands
        hst_bvi_band_red = 'F814W'
        hst_bvi_band_green = 'F555W'
        if 'F438W' in band_list:
            hst_bvi_band_blue = 'F438W'
        else:
            hst_bvi_band_blue = 'F435W'

        # load ha native and continuum subtracted
        if include_h_alpha:
            band_list += [ha_cont_sub_stamp_band]
            # get HST-H-alpha
            if 'F657N' in self.phangs_hst_ha_obs_band_dict[self.hst_ha_target_name]['ha_observed']:
                hst_habu_band_red = 'F657N'
                band_list += ['F657N']
                hst_ha_stamp_band = 'F657N'
            elif 'F658N' in self.phangs_hst_ha_obs_band_dict[self.hst_ha_target_name]['ha_observed']:
                hst_habu_band_red = 'F658N'
                band_list += ['F658N']
                hst_ha_stamp_band = 'F658N'
            else:
                hst_ha_stamp_band = None
                hst_habu_band_red = None
            hst_habu_band_green = hst_bvi_band_blue
            hst_habu_band_blue = 'F336W'
        else:
            hst_ha_stamp_band = None
            hst_habu_band_red = None
            hst_habu_band_green = None
            hst_habu_band_blue = None
        hst_stamp_band_list = self.sort_band_list(hst_stamp_band_list)
        if include_nircam:
            band_list += self.get_nircam_band_list()
            nircam_stamp_band_list = self.get_nircam_band_list()
            nircam_band_red = 'F300M'
            nircam_band_green = 'F335M'
            nircam_band_blue = 'F200W'
        else:
            nircam_stamp_band_list = []
            nircam_band_red = None
            nircam_band_green = None
            nircam_band_blue = None

        if include_miri:
            band_list += self.get_miri_band_list()
            miri_stamp_band_list = self.get_miri_band_list()
            miri_band_red = 'F1130W'
            miri_band_green = 'F1000W'
            miri_band_blue = 'F770W'
        else:
            miri_stamp_band_list = []
            miri_band_red = None
            miri_band_green = None
            miri_band_blue = None

        # load all bands into constructor
        self.load_phangs_bands(band_list=band_list, flux_unit='MJy/sr', load_err=False)

        # get hst_bvi_zoom_in
        img_hst_bvi_overview, wcs_hst_bvi_overview = self.get_target_overview_rgb_img(red_band=hst_bvi_band_red,
                                                                                      green_band=hst_bvi_band_green,
                                                                                      blue_band=hst_bvi_band_blue,
                                                                                      overview_img_size=(500, 500))

        img_hst_bvi_zoom_in, wcs_hst_bvi_zoom_in = self.get_rgb_zoom_in(ra=ra_region, dec=dec_region,
                                                                        cutout_size=env_cutout_size,
                                                                        circle_rad=circle_rad_region,
                                                                        band_red=hst_bvi_band_red,
                                                                        band_green=hst_bvi_band_green,
                                                                        band_blue=hst_bvi_band_blue)

        if include_h_alpha:
            img_hst_habu_zoom_in, wcs_hst_habu_zoom_in = self.get_rgb_zoom_in(ra=ra_region, dec=dec_region,
                                                                              cutout_size=env_cutout_size,
                                                                              circle_rad=circle_rad_region,
                                                                              band_red=hst_habu_band_red,
                                                                              band_green=hst_habu_band_green,
                                                                              band_blue=hst_habu_band_blue)
        else:
            img_hst_habu_zoom_in, wcs_hst_habu_zoom_in = None, None
        if include_nircam:
            img_nircam_zoom_in, wcs_nircam_zoom_in = self.get_rgb_zoom_in(ra=ra_region, dec=dec_region,
                                                                          cutout_size=env_cutout_size,
                                                                          circle_rad=circle_rad_region,
                                                                          band_red=nircam_band_red,
                                                                          band_green=nircam_band_green,
                                                                          band_blue=nircam_band_blue)
        else:
            img_nircam_zoom_in, wcs_nircam_zoom_in = None, None
        if include_miri:
            img_miri_zoom_in, wcs_miri_zoom_in = self.get_rgb_zoom_in(ra=ra_region, dec=dec_region,
                                                                      cutout_size=env_cutout_size,
                                                                      circle_rad=circle_rad_region,
                                                                      band_red=miri_band_red,
                                                                      band_green=miri_band_green,
                                                                      band_blue=miri_band_blue)
        else:
            img_miri_zoom_in, wcs_miri_zoom_in = None, None

        # load cutout stamps
        cutout_dict_stamp = self.get_band_cutout_dict(ra_cutout=ra_region, dec_cutout=dec_region,
                                                      cutout_size=stamp_cutout_size,
                                                      band_list=band_list)

        # plotting
        figure = plt.figure(figsize=(30, 50))

        # limits for the overview image
        overview_width = 0.45
        overview_height = 0.45
        overview_left_align = 0.035
        overview_bottom_align = 0.628

        # limits for the zoom in panels
        zoom_in_width = 0.24
        zoom_in_height = 0.24
        zoom_in_left_align = 0.51
        zoom_in_bottom_align = 0.66
        zoom_in_space_horizontal = -0.095
        zoom_in_space_vertical = 0.005

        # limits for the stamps
        stamp_width = 0.1
        stamp_height = 0.1
        stamp_left_align = 0.035
        stamp_bottom_align = 0.60
        stamp_space_horizontal = -0.02
        stamp_space_vertical = 0.005

        stamp_dict = {
            'hst_stamp_band_list': hst_stamp_band_list,
            'nircam_stamp_band_list': nircam_stamp_band_list,
            'miri_stamp_band_list': miri_stamp_band_list,
            'cutout_dict_stamp': cutout_dict_stamp, 'stamp_width': stamp_width, 'stamp_height': stamp_height,
            'stamp_left_align': stamp_left_align, 'stamp_bottom_align': stamp_bottom_align,
            'stamp_space_horizontal': stamp_space_horizontal, 'stamp_space_vertical': stamp_space_vertical,
            'hst_ha_stamp_band': hst_ha_stamp_band, 'ha_cont_sub_stamp_band': ha_cont_sub_stamp_band,
            'cbar_log_scale': cbar_log_scale
        }

        ax_hst_bvi_overview = figure.add_axes([overview_left_align, overview_bottom_align,
                                               overview_width, overview_height],
                                              projection=wcs_hst_bvi_overview)

        ax_hst_bvi_zoom_in = figure.add_axes([zoom_in_left_align,
                                              zoom_in_bottom_align + zoom_in_height + zoom_in_space_horizontal,
                                              zoom_in_width, zoom_in_height], projection=wcs_hst_bvi_zoom_in)
        hst_overview_dict = {
            'ax_hst_bvi_overview': ax_hst_bvi_overview,
            'img_hst_bvi_overview': img_hst_bvi_overview,
            'wcs_hst_bvi_overview': wcs_hst_bvi_overview,
            'hst_bvi_band_red': hst_bvi_band_red,
            'hst_bvi_band_green': hst_bvi_band_green,
            'hst_bvi_band_blue': hst_bvi_band_blue,
        }
        hst_zoom_in_dict = {
            'ax_hst_bvi_zoom_in': ax_hst_bvi_zoom_in,
            'img_hst_bvi_zoom_in': img_hst_bvi_zoom_in,
            'wcs_hst_bvi_zoom_in': wcs_hst_bvi_zoom_in,
            'hst_bvi_band_red': hst_bvi_band_red,
            'hst_bvi_band_green': hst_bvi_band_green,
            'hst_bvi_band_blue': hst_bvi_band_blue,
        }

        if include_h_alpha:
            ax_hst_habu_zoom_in = figure.add_axes([zoom_in_left_align + zoom_in_width + zoom_in_space_vertical,
                                                   zoom_in_bottom_align + zoom_in_height + zoom_in_space_horizontal,
                                                   zoom_in_width, zoom_in_height], projection=wcs_hst_habu_zoom_in)
            hst_habu_zoom_in_dict = {
                'ax_hst_habu_zoom_in': ax_hst_habu_zoom_in,
                'img_hst_habu_zoom_in': img_hst_habu_zoom_in,
                'wcs_hst_habu_zoom_in': wcs_hst_habu_zoom_in,
                'hst_habu_band_red': hst_habu_band_red,
                'hst_habu_band_green': hst_habu_band_green,
                'hst_habu_band_blue': hst_habu_band_blue,
            }
        else:
            hst_habu_zoom_in_dict = None
        if include_nircam & (img_nircam_zoom_in is not None):
            ax_nircam_zoom_in = figure.add_axes([zoom_in_left_align,
                                                 zoom_in_bottom_align, zoom_in_width, zoom_in_height],
                                                projection=wcs_nircam_zoom_in)
            nircam_zoom_in_dict = {
                'ax_nircam_zoom_in': ax_nircam_zoom_in,
                'img_nircam_zoom_in': img_nircam_zoom_in,
                'wcs_nircam_zoom_in': wcs_nircam_zoom_in,
                'nircam_band_red': nircam_band_red,
                'nircam_band_green': nircam_band_green,
                'nircam_band_blue': nircam_band_blue,
            }
        else:
            nircam_zoom_in_dict = None
        if include_miri & (img_miri_zoom_in is not None):
            ax_miri_zoom_in = figure.add_axes([zoom_in_left_align + zoom_in_width + zoom_in_space_vertical, zoom_in_bottom_align,
                                               zoom_in_width, zoom_in_height], projection=wcs_miri_zoom_in)
            miri_zoom_in_dict = {
                'ax_miri_zoom_in': ax_miri_zoom_in,
                'img_miri_zoom_in': img_miri_zoom_in,
                'wcs_miri_zoom_in': wcs_miri_zoom_in,
                'miri_band_red': miri_band_red,
                'miri_band_green': miri_band_green,
                'miri_band_blue': miri_band_blue,
            }
        else:
            miri_zoom_in_dict = None

        # get 50 and 80 % encirceled energy radius
        radius_dict = {}
        for band in band_list:
            if band == 'ha':
                continue
            if band in self.phangs_hst_obs_band_dict[self.target_name]['wfc3_uvis_observed_bands']:
                radius_dict.update({'aperture_%s_ee50' % band: self.hst_encircle_apertures_wfc3_uvis2_arcsec[band]['ee50']})
                radius_dict.update({'aperture_%s_ee80' % band: self.hst_encircle_apertures_wfc3_uvis2_arcsec[band]['ee80']})
            elif band in self.phangs_hst_obs_band_dict[self.target_name]['acs_wfc1_observed_bands']:
                radius_dict.update({'aperture_%s_ee50' % band: self.hst_encircle_apertures_acs_wfc1_arcsec[band]['ee50']})
                radius_dict.update({'aperture_%s_ee80' % band: self.hst_encircle_apertures_acs_wfc1_arcsec[band]['ee80']})
            elif band in ['F657N', 'F658N']:
                radius_dict.update({'aperture_%s_ee50' % band: self.hst_encircle_apertures_wfc3_uvis2_arcsec[band]['ee50']})
                radius_dict.update({'aperture_%s_ee80' % band: self.hst_encircle_apertures_wfc3_uvis2_arcsec[band]['ee80']})
            elif band in self.nircam_bands:
                radius_dict.update({'aperture_%s_ee50' % band: self.nircam_encircle_apertures_arcsec[band]['ee50']})
                radius_dict.update({'aperture_%s_ee80' % band: self.nircam_encircle_apertures_arcsec[band]['ee80']})
            elif band in self.miri_bands:
                radius_dict.update({'aperture_%s_ee50' % band: self.miri_encircle_apertures_arcsec[band]['ee50']})
                radius_dict.update({'aperture_%s_ee80' % band: self.miri_encircle_apertures_arcsec[band]['ee80']})

        self.plot_over_view_zoom_in_panel(figure=figure, ra_region=ra_region, dec_region=dec_region,
                                          hst_overview_dict=hst_overview_dict, hst_zoom_in_dict=hst_zoom_in_dict,
                                          stamp_dict=stamp_dict, hst_habu_zoom_in_dict=hst_habu_zoom_in_dict,
                                          nircam_zoom_in_dict=nircam_zoom_in_dict, miri_zoom_in_dict=miri_zoom_in_dict,
                                          radius_dict=None,
                                          env_cutout_size=env_cutout_size, circle_rad_region=None,
                                          stamp_cutout_size=stamp_cutout_size, fontsize_large=fontsize_large,
                                          fontsize_small=fontsize_small)

        # plot SED
        # load large cutout dict for flux density estimation
        self.change_phangs_band_units(band_list=band_list, new_unit='mJy')
        cutout_dict_zoom_in = self.get_band_cutout_dict(ra_cutout=ra_region, dec_cutout=dec_region,
                                                        cutout_size=env_cutout_size,
                                                        band_list=band_list)


        # get fluxes from circular aperture
        flux_dict_ee50 = {}
        flux_dict_ee80 = {}
        pos_obj = SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg)
        for band in band_list:
            if band == 'ha':
                continue
            if cutout_dict_zoom_in['%s_img_cutout' % band].data is None:
                continue
            flux_ee50, flux_ee50_err = helper_func.extract_flux_from_circ_aperture(
                data=cutout_dict_zoom_in['%s_img_cutout' % band].data,
                wcs=cutout_dict_zoom_in['%s_img_cutout' % band].wcs,
                pos=pos_obj,
                aperture_rad=radius_dict['aperture_%s_ee50' % band])
            flux_dict_ee50.update({'flux_%s' % band: flux_ee50, 'flux_err_%s' % band: flux_ee50_err})
            flux_ee80, flux_ee80_err = helper_func.extract_flux_from_circ_aperture(
                data=cutout_dict_zoom_in['%s_img_cutout' % band].data,
                wcs=cutout_dict_zoom_in['%s_img_cutout' % band].wcs,
                pos=pos_obj,
                aperture_rad=radius_dict['aperture_%s_ee80' % band])
            flux_dict_ee80.update({'flux_%s' % band: flux_ee80, 'flux_err_%s' % band: flux_ee80_err})

        # plot flux-density distribution
        ax_sed = figure.add_axes([0.065, 0.32, 0.9, 0.15])

        wave_ee50_list = []
        wave_ee80_list = []
        flux_ee50_list = []
        flux_ee80_list = []
        for band in band_list:
            if 'ha' in band:
                continue
            if 'flux_%s' % band not in flux_dict_ee50:
                continue
            if band in self.get_hst_band_list():
                color = 'tab:blue'
            elif band in self.get_hst_ha_band_list():
                color = 'tab:red'
            elif band in self.get_nircam_band_list():
                color = 'tab:green'
            elif band in self.get_miri_band_list():
                color = 'tab:purple'
            elif band in self.get_astrosat_band_list():
                color = 'tab:pink'
            else:
                color = ''
            mean_wave = self.get_band_wave(band=band)
            min_wave = self.get_band_wave(band=band, wave_estimator='min_wave')
            max_wave = self.get_band_wave(band=band, wave_estimator='max_wave')
            if ((flux_dict_ee50['flux_%s' % band] < 0) |
                    (flux_dict_ee50['flux_%s' % band] < 3*flux_dict_ee50['flux_err_%s' % band])):
                ax_sed.errorbar(mean_wave, 3*flux_dict_ee50['flux_err_%s' % band],
                                yerr=flux_dict_ee50['flux_err_%s' % band],
                                xerr=[[mean_wave-min_wave], [max_wave-mean_wave]],
                                ecolor=color, elinewidth=5, capsize=10, uplims=True, xlolims=False)
            else:
                wave_ee50_list.append(mean_wave)
                flux_ee50_list.append(flux_dict_ee50['flux_%s' % band])
                ax_sed.errorbar(mean_wave, flux_dict_ee50['flux_%s' % band],
                                         xerr=[[mean_wave-min_wave], [max_wave-mean_wave]],
                                         yerr=flux_dict_ee50['flux_err_%s' % band],
                                         fmt='o', color=color, ms=20)

            if ((flux_dict_ee80['flux_%s' % band] < 0) |
                    (flux_dict_ee80['flux_%s' % band] < 3*flux_dict_ee80['flux_err_%s' % band])):
                ax_sed.errorbar(mean_wave, 3*flux_dict_ee80['flux_err_%s' % band],
                                yerr=flux_dict_ee80['flux_err_%s' % band],
                                xerr=[[mean_wave-min_wave], [max_wave-mean_wave]],
                                ecolor=color, elinewidth=5, capsize=10, uplims=True, xlolims=False)
            else:
                wave_ee80_list.append(mean_wave)
                flux_ee80_list.append(flux_dict_ee80['flux_%s' % band])
                ax_sed.errorbar(mean_wave, flux_dict_ee80['flux_%s' % band],
                                         xerr=[[mean_wave-min_wave], [max_wave-mean_wave]],
                                         yerr=flux_dict_ee80['flux_err_%s' % band],
                                         fmt='*', color=color, ms=20)

        sort_ee50 = np.argsort(wave_ee50_list)
        sort_ee80 = np.argsort(wave_ee80_list)
        ax_sed.plot(np.array(wave_ee50_list)[sort_ee50], np.array(flux_ee50_list)[sort_ee50], color='k', linewidth=2, label='EE 50%')
        ax_sed.plot(np.array(wave_ee80_list)[sort_ee80], np.array(flux_ee80_list)[sort_ee80], color='gray', linewidth=2, label='EE 80%')
        ax_sed.legend(frameon=False, fontsize=fontsize_large)

        ax_sed.set_yscale('log')
        ax_sed.set_xscale('log')
        ax_sed.set_xlabel(r'Wavelength [$\mu$m]', fontsize=fontsize_large)
        ax_sed.set_ylabel(r'Flux [mJy]', fontsize=fontsize_large)
        ax_sed.tick_params(axis='both', which='both', width=2, direction='in', labelsize=fontsize_large)

        # put up strings as title
        title = None
        for title_string in string_dict.keys():
            if title is None:
                title = ''
            else:
                title += '\n'
            title += string_dict[title_string]

        ax_sed.set_title(title, fontsize=fontsize_large, loc='left')

        cutout_alma = helper_func.get_img_cutout(img=mol_gas_dens_dict['mol_gas_dens_data'],
                                                 wcs=mol_gas_dens_dict['mol_gas_dens_wcs'],
                                                 coord=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                                 cutout_size=env_cutout_size)
        if not np.isnan(cutout_alma.data).all():
            ax_alma = figure.add_axes([0.05, -0.05, 0.4, 0.4], projection=cutout_alma.wcs)
            ax_cbar_alma = figure.add_axes([0.1, 0.28, 0.3, 0.01])
            cmap_alma = 'inferno'

            min_alma_value = np.nanmin(cutout_alma.data)
            max_alma_value = np.nanmax(cutout_alma.data)
            if min_alma_value <= 0:
                min_alma_value = max_alma_value / 100
            norm = LogNorm(min_alma_value, max_alma_value)
            helper_func.create_cbar(ax_cbar=ax_cbar_alma, cmap=cmap_alma, norm=norm,
                                    cbar_label=r'log($\Sigma_{\rm H2}$/[M$_{\odot}$ kpc$^{-2}$])', fontsize=fontsize_large,
                                    ticks=None, labelpad=2, tick_width=2, orientation='horizontal', extend='neither')

            ax_alma.imshow(cutout_alma.data, norm=norm, cmap=cmap_alma)
            coord_cloud_pix = cutout_alma.wcs.world_to_pixel(mol_gas_dens_dict['coord_cloud_world'])
            for idx, pos in enumerate(mol_gas_dens_dict['coord_cloud_world']):
                rad = mol_gas_dens_dict['rad_cloud_arcsec'][idx]
                if np.isnan(rad):
                    continue
                if ((coord_cloud_pix[0][idx] < 0) | (coord_cloud_pix[1][idx] < 0) |
                        (coord_cloud_pix[0][idx] > cutout_alma.data.shape[0]) |
                        (coord_cloud_pix[1][idx] > cutout_alma.data.shape[1]) ):
                    continue
                helper_func.plot_coord_circle(ax=ax_alma, pos=pos, rad=rad, color='blue', line_style='-', line_width=2)

            # central_coord_pixel = cutout_alma.wcs.world_to_pixel(SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg))
            # ax_alma.scatter(central_coord_pixel[0], central_coord_pixel[1], s=180, marker='*', color='k')
            # ax_alma.set_xlim(0, cutout_alma.data.shape[0])
            # ax_alma.set_ylim(0, cutout_alma.data.shape[1])
            helper_func.plot_coord_circle(ax=ax_alma, pos=SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg),
                                          rad=0.1, color='k', line_style='-', line_width=2, fill=True)
            self.arr_axis_params(ax=ax_alma, fontsize=fontsize_large, labelsize=fontsize_large)
            # ax_alma.set_title(info_string, fontsize=fontsize_large)


        ax_ccd = figure.add_axes([0.55, 0.03, 0.4, 0.25])

        vmax = np.nanmax(color_color_dict['gauss_dict_ubvi_hum_12']['gauss_map'])
        ax_ccd.imshow(color_color_dict['gauss_dict_ubvi_hum_12']['gauss_map'], origin='lower',
                             extent=(color_color_dict['x_lim_vi'][0], color_color_dict['x_lim_vi'][1],
                                     color_color_dict['y_lim_ub'][1], color_color_dict['y_lim_ub'][0]),
                             interpolation='nearest', aspect='auto', cmap='Greys', vmin=0+vmax/10, vmax=vmax/1.1)



        ax_ccd.plot(color_color_dict['vi_color_ycl_hull'], color_color_dict['ub_color_ycl_hull'], color='blue', linewidth=3)
        ax_ccd.plot(color_color_dict['vi_color_map_hull'], color_color_dict['ub_color_map_hull'], color='green', linewidth=3)
        ax_ccd.plot(color_color_dict['vi_color_ogcc_hull'], color_color_dict['ub_color_ogcc_hull'], color='red', linewidth=3)

        ax_ccd.plot(color_color_dict['model_vi_sol'], color_color_dict['model_ub_sol'], color='tab:cyan', linewidth=5, label=r'BC03, Z$_{\odot}$')
        ax_ccd.scatter(color_color_dict['model_vi_sol'][color_color_dict['age_mod_sol'] == 1],
                    color_color_dict['model_ub_sol'][color_color_dict['age_mod_sol'] == 1], color='k', s=150)
        ax_ccd.scatter(color_color_dict['model_vi_sol'][color_color_dict['age_mod_sol'] == 5],
                    color_color_dict['model_ub_sol'][color_color_dict['age_mod_sol'] == 5], color='k', s=150)
        ax_ccd.scatter(color_color_dict['model_vi_sol'][color_color_dict['age_mod_sol'] == 10],
                    color_color_dict['model_ub_sol'][color_color_dict['age_mod_sol'] == 10], color='k', s=150)
        ax_ccd.scatter(color_color_dict['model_vi_sol'][color_color_dict['age_mod_sol'] == 100],
                    color_color_dict['model_ub_sol'][color_color_dict['age_mod_sol'] == 100], color='k', s=150)
        ax_ccd.scatter(color_color_dict['model_vi_sol'][color_color_dict['age_mod_sol'] == 1000],
                    color_color_dict['model_ub_sol'][color_color_dict['age_mod_sol'] == 1000], color='k', s=150)
        ax_ccd.scatter(color_color_dict['model_vi_sol'][color_color_dict['age_mod_sol'] == 13750],
                    color_color_dict['model_ub_sol'][color_color_dict['age_mod_sol'] == 13750], color='k', s=150)

        ax_ccd.text(color_color_dict['model_vi_sol'][color_color_dict['age_mod_sol'] == 1],
                    color_color_dict['model_ub_sol'][color_color_dict['age_mod_sol'] == 1] - 0.05,
                    '1 Myr', horizontalalignment='center', verticalalignment='bottom',
                    color='k', fontsize=fontsize_large)
        ax_ccd.text(color_color_dict['model_vi_sol'][color_color_dict['age_mod_sol'] == 5] - 0.05,
                    color_color_dict['model_ub_sol'][color_color_dict['age_mod_sol'] == 5],
                    '5 Myr', horizontalalignment='right', verticalalignment='center',
                    color='k', fontsize=fontsize_large)
        ax_ccd.text(color_color_dict['model_vi_sol'][color_color_dict['age_mod_sol'] == 10] + 0.05,
                    color_color_dict['model_ub_sol'][color_color_dict['age_mod_sol'] == 10],
                    '10 Myr', horizontalalignment='left', verticalalignment='center',
                    color='k', fontsize=fontsize_large)
        ax_ccd.text(color_color_dict['model_vi_sol'][color_color_dict['age_mod_sol'] == 100] - 0.05,
                    color_color_dict['model_ub_sol'][color_color_dict['age_mod_sol'] == 100],
                    '100 Myr', horizontalalignment='right', verticalalignment='center',
                    color='k', fontsize=fontsize_large)
        ax_ccd.text(color_color_dict['model_vi_sol'][color_color_dict['age_mod_sol'] == 1000],
                    color_color_dict['model_ub_sol'][color_color_dict['age_mod_sol'] == 1000] - 0.05,
                    '1 Gyr', horizontalalignment='center', verticalalignment='bottom',
                    color='k', fontsize=fontsize_large)
        ax_ccd.text(color_color_dict['model_vi_sol'][color_color_dict['age_mod_sol'] == 13750],
                    color_color_dict['model_ub_sol'][color_color_dict['age_mod_sol'] == 13750] + 0.05,
                    '13.7 Gyr', horizontalalignment='center', verticalalignment='top',
                    color='k', fontsize=fontsize_large)

        ax_ccd.errorbar(color_color_dict['color_vi'], color_color_dict['color_ub'],
                        xerr=color_color_dict['color_vi_err'],
                        yerr=color_color_dict['color_ub_err'],
                        fmt='o', ms=10, capsize=5, markeredgewidth=2, elinewidth=3, color='red')

        helper_func.plot_reddening_vect(ax=ax_ccd,
                        x_color_int=0.8, y_color_int=-1.6, av_val=1,
                        linewidth=3, line_color='k',
                        text=True, fontsize=fontsize_large, text_color='k', x_text_offset=0.1, y_text_offset=-0.05)

        ax_ccd.legend(frameon=False, loc=3, fontsize=fontsize_large)
        ax_ccd.set_title(string_dict['color_color_str'], fontsize=fontsize_large)
        ax_ccd.set_xlim(color_color_dict['x_lim_vi'])
        ax_ccd.set_ylim(color_color_dict['y_lim_ub'])

        ax_ccd.set_ylabel('U (F336W) - B (F438W/F435W'+'$^*$'+')', fontsize=fontsize_large)
        ax_ccd.set_xlabel('V (F555W) - I (F814W)', fontsize=fontsize_large)

        ax_ccd.tick_params(axis='both', which='both', width=1.5, length=4, right=True, top=True, direction='in', labelsize=fontsize_large)


        # color_color_dict


        return figure

    def plot_ultimate_zoom_in_panel(self, ra_region, dec_region, include_h_alpha=True, ha_cont_sub_stamp_band='ha',
                                    include_nircam=True,
                                    include_miri=True,
                                    muse_fit_dict=None,
                                    env_cutout_size=(10, 10), circle_rad_region=1.25, circle_rad_obj=0.4,
                                    stamp_cutout_size=(2.5, 2.5),
                                    overview_img_pixel_dim=(500, 500),
                                    cbar_log_scale=True,
                                    fontsize_large=30, fontsize_small=20,
                                    ):

        # load all bands needed
        band_list = []
        # load HST bands
        band_list += self.get_hst_band_list()
        hst_stamp_band_list = band_list.copy()
        # get BVI filters
        # specify color the bands
        hst_bvi_band_red = 'F814W'
        hst_bvi_band_green = 'F555W'
        if 'F438W' in band_list:
            hst_bvi_band_blue = 'F438W'
        else:
            hst_bvi_band_blue = 'F435W'

        # load ha native and continuum subtracted
        if include_h_alpha:
            band_list += [ha_cont_sub_stamp_band]
            # get HST-H-alpha
            if 'F657N' in self.phangs_hst_ha_obs_band_dict[self.hst_ha_target_name]['ha_observed']:
                hst_habu_band_red = 'F657N'
                band_list += ['F657N']
                hst_ha_stamp_band = 'F657N'
            elif 'F658N' in self.phangs_hst_ha_obs_band_dict[self.hst_ha_target_name]['ha_observed']:
                hst_habu_band_red = 'F658N'
                band_list += ['F658N']
                hst_ha_stamp_band = 'F658N'
            else:
                hst_ha_stamp_band = None
                hst_habu_band_red = None
            hst_habu_band_green = hst_bvi_band_blue
            hst_habu_band_blue = 'F336W'
        else:
            hst_ha_stamp_band = None
            hst_habu_band_red = None
            hst_habu_band_green = None
            hst_habu_band_blue = None
        hst_stamp_band_list = self.sort_band_list(hst_stamp_band_list)
        if include_nircam:
            band_list += self.get_nircam_band_list()
            nircam_stamp_band_list = self.get_nircam_band_list()
            nircam_band_red = 'F300M'
            nircam_band_green = 'F335M'
            nircam_band_blue = 'F200W'
        else:
            nircam_stamp_band_list = []
            nircam_band_red = None
            nircam_band_green = None
            nircam_band_blue = None

        if include_miri:
            band_list += self.get_miri_band_list()
            miri_stamp_band_list = self.get_miri_band_list()
            miri_band_red = 'F1130W'
            miri_band_green = 'F1000W'
            miri_band_blue = 'F770W'
        else:
            miri_stamp_band_list = []
            miri_band_red = None
            miri_band_green = None
            miri_band_blue = None

        # include astrosat FUV observation
        # identify the bluest astrosat filter if available
        if self.target_name in self.astrosat_targets.keys():
            available_astrosat_bands = self.astrosat_targets[self.target_name]['observed_bands']
            if isinstance(available_astrosat_bands, list):
                bluest_astrosat_band = available_astrosat_bands[0]
            else:
                bluest_astrosat_band = available_astrosat_bands
            if bluest_astrosat_band in self.astrosat_fuv_bands:
                astrosat_band = bluest_astrosat_band
                band_list += [astrosat_band]
            else:
                astrosat_band = None
        else:
            astrosat_band = None

        # load all bands into constructor
        self.load_phangs_bands(band_list=band_list, flux_unit='MJy/sr', load_err=False)

        #define the position of object
        pos_obj = SkyCoord(ra=ra_region*u.deg, dec=dec_region*u.deg)
        # get hst_bvi_zoom_in
        img_hst_bvi_overview, wcs_hst_bvi_overview = self.get_target_overview_rgb_img(red_band=hst_bvi_band_red,
                                                                                      green_band=hst_bvi_band_green,
                                                                                      blue_band=hst_bvi_band_blue,
                                                                                      overview_img_size=
                                                                                      overview_img_pixel_dim)
        img_hst_bvi_zoom_in, wcs_hst_bvi_zoom_in = self.get_rgb_zoom_in(ra=ra_region, dec=dec_region,
                                                                        cutout_size=env_cutout_size,
                                                                        circle_rad=circle_rad_region,
                                                                        band_red=hst_bvi_band_red,
                                                                        band_green=hst_bvi_band_green,
                                                                        band_blue=hst_bvi_band_blue)

        if include_h_alpha:
            img_hst_habu_zoom_in, wcs_hst_habu_zoom_in = self.get_rgb_zoom_in(ra=ra_region, dec=dec_region,
                                                                              cutout_size=env_cutout_size,
                                                                              circle_rad=circle_rad_region,
                                                                              band_red=hst_habu_band_red,
                                                                              band_green=hst_habu_band_green,
                                                                              band_blue=hst_habu_band_blue)
        else:
            img_hst_habu_zoom_in, wcs_hst_habu_zoom_in = None, None
        if include_nircam:
            img_nircam_zoom_in, wcs_nircam_zoom_in = self.get_rgb_zoom_in(ra=ra_region, dec=dec_region,
                                                                          cutout_size=env_cutout_size,
                                                                          circle_rad=circle_rad_region,
                                                                          band_red=nircam_band_red,
                                                                          band_green=nircam_band_green,
                                                                          band_blue=nircam_band_blue)
        else:
            img_nircam_zoom_in, wcs_nircam_zoom_in = None, None
        if include_miri:
            img_miri_zoom_in, wcs_miri_zoom_in = self.get_rgb_zoom_in(ra=ra_region, dec=dec_region,
                                                                      cutout_size=env_cutout_size,
                                                                      circle_rad=circle_rad_region,
                                                                      band_red=miri_band_red,
                                                                      band_green=miri_band_green,
                                                                      band_blue=miri_band_blue)
        else:
            img_miri_zoom_in, wcs_miri_zoom_in = None, None

        # load cutout stamps
        cutout_dict_stamp = self.get_band_cutout_dict(ra_cutout=ra_region, dec_cutout=dec_region,
                                                      cutout_size=stamp_cutout_size,
                                                      band_list=band_list)

        # load large cutout dict for flux density estimation
        self.change_phangs_band_units(band_list=band_list, new_unit='erg A-1 cm-2 s-1')
        cutout_dict_zoom_in = self.get_band_cutout_dict(ra_cutout=ra_region, dec_cutout=dec_region,
                                                        cutout_size=env_cutout_size,
                                                        band_list=band_list)

        # get fluxes from circular aperture
        flux_dict = {}
        flux_dict_small_target = {}
        for band in band_list:
            if band == 'ha':
                continue
            if cutout_dict_zoom_in['%s_img_cutout' % band].data is None:
                continue
            flux, flux_err = helper_func.extract_flux_from_circ_aperture(
                data=cutout_dict_zoom_in['%s_img_cutout' % band].data,
                wcs=cutout_dict_zoom_in['%s_img_cutout' % band].wcs,
                pos=pos_obj, aperture_rad=circle_rad_region, data_err=None)
            flux_dict.update({'flux_%s' % band: flux, 'flux_err_%s' % band: flux_err})

            flux_small, flux_small_err = helper_func.extract_flux_from_circ_aperture(
                data=cutout_dict_zoom_in['%s_img_cutout' % band].data,
                wcs=cutout_dict_zoom_in['%s_img_cutout' % band].wcs,
                pos=pos_obj, aperture_rad=circle_rad_obj, data_err=None)
            flux_dict_small_target.update({'flux_%s' % band: flux_small, 'flux_err_%s' % band: flux_small_err})

        # plotting
        figure = plt.figure(figsize=(30, 50))


        # limits for the overview image
        overview_width = 0.45
        overview_height = 0.45
        overview_left_align = 0.035
        overview_bottom_align = 0.628

        # limits for the zoom in panels
        zoom_in_width = 0.24
        zoom_in_height = 0.24
        zoom_in_left_align = 0.51
        zoom_in_bottom_align = 0.66
        zoom_in_space_horizontal = -0.095
        zoom_in_space_vertical = 0.005

        # limits for the stamps
        stamp_width = 0.1
        stamp_height = 0.1
        stamp_left_align = 0.035
        stamp_bottom_align = 0.60
        stamp_space_horizontal = -0.02
        stamp_space_vertical = 0.005


        ax_hst_bvi_overview = figure.add_axes([overview_left_align, overview_bottom_align,
                                               overview_width, overview_height],
                                              projection=wcs_hst_bvi_overview)

        ax_hst_bvi_zoom_in = figure.add_axes([zoom_in_left_align,
                                              zoom_in_bottom_align + zoom_in_height + zoom_in_space_horizontal,
                                              zoom_in_width, zoom_in_height], projection=wcs_hst_bvi_zoom_in)

        if include_h_alpha:
            ax_hst_habu_zoom_in = figure.add_axes([zoom_in_left_align + zoom_in_width + zoom_in_space_vertical,
                                                   zoom_in_bottom_align + zoom_in_height + zoom_in_space_horizontal,
                                                   zoom_in_width, zoom_in_height], projection=wcs_hst_habu_zoom_in)
        else:
            ax_hst_habu_zoom_in = None
        if include_nircam:
            ax_nircam_zoom_in = figure.add_axes([zoom_in_left_align,
                                                 zoom_in_bottom_align, zoom_in_width, zoom_in_height],
                                                projection=wcs_nircam_zoom_in)
        else:
            ax_nircam_zoom_in = None
        if include_miri:
            ax_miri_zoom_in = figure.add_axes([zoom_in_left_align + zoom_in_width + zoom_in_space_vertical, zoom_in_bottom_align,
                                               zoom_in_width, zoom_in_height], projection=wcs_miri_zoom_in)
        else:
            ax_miri_zoom_in = None

        # plot the RGB images
        ax_hst_bvi_overview.imshow(img_hst_bvi_overview)
        pe = [patheffects.withStroke(linewidth=3, foreground="w")]
        ax_hst_bvi_overview.text(0.02, 0.98, 'HST', horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='white',
                                 transform=ax_hst_bvi_overview.transAxes, path_effects=pe)
        ax_hst_bvi_overview.text(0.02, 0.95, hst_bvi_band_red.upper(), horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='red',
                                 transform=ax_hst_bvi_overview.transAxes, path_effects=pe)
        ax_hst_bvi_overview.text(0.02, 0.92, hst_bvi_band_green.upper(), horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='green',
                                 transform=ax_hst_bvi_overview.transAxes, path_effects=pe)
        ax_hst_bvi_overview.text(0.02, 0.89, hst_bvi_band_blue.upper(), horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='blue',
                                 transform=ax_hst_bvi_overview.transAxes, path_effects=pe)
        ax_hst_bvi_overview.set_title(self.target_name.upper(), fontsize=fontsize_large + 10)
        self.arr_axis_params(ax=ax_hst_bvi_overview, ra_tick_label=True, dec_tick_label=True,
                        ra_axis_label='R.A. (2000.0)', dec_axis_label='DEC. (2000.0)',
                        ra_minpad=0.8, dec_minpad=0.8, tick_color='white', label_color='k',
                        fontsize=fontsize_large, labelsize=fontsize_large, ra_tick_num=3, dec_tick_num=3)
        helper_func.draw_box(ax=ax_hst_bvi_overview, wcs=wcs_hst_bvi_overview,
                             coord=pos_obj,
                             box_size=env_cutout_size, color='red', line_width=2, line_style='--')
        helper_func.plot_coord_circle(ax=ax_hst_bvi_zoom_in, pos=pos_obj,
                                      rad=circle_rad_region, color='white', line_style='-', line_width=2, alpha=1., fill=False)

        ax_hst_bvi_zoom_in.imshow(img_hst_bvi_zoom_in)
        ax_hst_bvi_zoom_in.text(0.02, 0.98, 'HST', horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='white',
                                 transform=ax_hst_bvi_zoom_in.transAxes, path_effects=pe)
        ax_hst_bvi_zoom_in.text(0.02, 0.92, hst_bvi_band_red.upper(), horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='red',
                                 transform=ax_hst_bvi_zoom_in.transAxes, path_effects=pe)
        ax_hst_bvi_zoom_in.text(0.02, 0.86, hst_bvi_band_green.upper(), horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='green',
                                 transform=ax_hst_bvi_zoom_in.transAxes, path_effects=pe)
        ax_hst_bvi_zoom_in.text(0.02, 0.80, hst_bvi_band_blue.upper(), horizontalalignment='left', verticalalignment='top',
                                 fontsize=fontsize_large, color='blue',
                                 transform=ax_hst_bvi_zoom_in.transAxes, path_effects=pe)
        helper_func.draw_box(ax=ax_hst_bvi_zoom_in, wcs=wcs_hst_bvi_zoom_in,
                             coord=pos_obj,
                             box_size=stamp_cutout_size, color='red', line_width=2, line_style='--')
        self.arr_axis_params(ax=ax_hst_bvi_zoom_in, ra_tick_label=False, dec_tick_label=True,
                        ra_axis_label=' ', dec_axis_label=' ',
                        ra_minpad=0.8, dec_minpad=0.8, tick_color='white', label_color='k',
                        fontsize=fontsize_large, labelsize=fontsize_large, ra_tick_num=3, dec_tick_num=3)
        if include_h_alpha:
            ax_hst_habu_zoom_in.imshow(img_hst_habu_zoom_in)
            self.arr_axis_params(ax=ax_hst_habu_zoom_in, ra_tick_label=False, dec_tick_label=False,
                            ra_axis_label=' ', dec_axis_label=' ',
                            ra_minpad=0.8, dec_minpad=0.8, tick_color='white', label_color='k',
                            fontsize=fontsize_large, labelsize=fontsize_large, ra_tick_num=3, dec_tick_num=3)
            ax_hst_habu_zoom_in.imshow(img_hst_habu_zoom_in)
            ax_hst_habu_zoom_in.text(0.02, 0.98, 'HST', horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='white',
                                     transform=ax_hst_habu_zoom_in.transAxes, path_effects=pe)
            ax_hst_habu_zoom_in.text(0.02, 0.92, hst_habu_band_red.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='red',
                                     transform=ax_hst_habu_zoom_in.transAxes, path_effects=pe)
            ax_hst_habu_zoom_in.text(0.02, 0.86, hst_habu_band_green.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='green',
                                     transform=ax_hst_habu_zoom_in.transAxes, path_effects=pe)
            ax_hst_habu_zoom_in.text(0.02, 0.80, hst_habu_band_blue.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='blue',
                                     transform=ax_hst_habu_zoom_in.transAxes, path_effects=pe)
            helper_func.draw_box(ax=ax_hst_habu_zoom_in, wcs=wcs_hst_habu_zoom_in,
                                 coord=pos_obj,
                                 box_size=stamp_cutout_size, color='red', line_width=2, line_style='--')
            helper_func.plot_coord_circle(ax=ax_hst_habu_zoom_in, pos=pos_obj,
                                          rad=circle_rad_region, color='white', line_style='-', line_width=2, alpha=1., fill=False)
        if include_nircam and (img_nircam_zoom_in is not None):
            ax_nircam_zoom_in.imshow(img_nircam_zoom_in)
            self.arr_axis_params(ax=ax_nircam_zoom_in, ra_tick_label=True, dec_tick_label=True,
                            ra_axis_label=' ', dec_axis_label=' ',
                            ra_minpad=0.8, dec_minpad=0.8, tick_color='white', label_color='k',
                            fontsize=fontsize_large, labelsize=fontsize_large, ra_tick_num=3, dec_tick_num=3)
            ax_nircam_zoom_in.text(0.02, 0.98, 'NIRCAM', horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='white',
                                     transform=ax_nircam_zoom_in.transAxes, path_effects=pe)
            ax_nircam_zoom_in.text(0.02, 0.92, nircam_band_red.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='red',
                                     transform=ax_nircam_zoom_in.transAxes, path_effects=pe)
            ax_nircam_zoom_in.text(0.02, 0.86, nircam_band_green.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='green',
                                     transform=ax_nircam_zoom_in.transAxes, path_effects=pe)
            ax_nircam_zoom_in.text(0.02, 0.80, nircam_band_blue.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='blue',
                                     transform=ax_nircam_zoom_in.transAxes, path_effects=pe)
            helper_func.draw_box(ax=ax_nircam_zoom_in, wcs=wcs_nircam_zoom_in,
                                 coord=pos_obj,
                                 box_size=stamp_cutout_size, color='red', line_width=2, line_style='--')
            helper_func.plot_coord_circle(ax=ax_nircam_zoom_in, pos=pos_obj,
                                          rad=circle_rad_region, color='white', line_style='-', line_width=2, alpha=1., fill=False)
        if include_miri and (img_miri_zoom_in is not None):
            ax_miri_zoom_in.imshow(img_miri_zoom_in)
            self.arr_axis_params(ax=ax_miri_zoom_in, ra_tick_label=True, dec_tick_label=False,
                            ra_axis_label=' ', dec_axis_label=' ',
                            ra_minpad=0.8, dec_minpad=0.8, tick_color='white', label_color='k',
                            fontsize=fontsize_large, labelsize=fontsize_large, ra_tick_num=3, dec_tick_num=3)
            ax_miri_zoom_in.text(0.02, 0.98, 'MIRI', horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='white',
                                     transform=ax_miri_zoom_in.transAxes, path_effects=pe)
            ax_miri_zoom_in.text(0.02, 0.92, miri_band_red.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='red',
                                     transform=ax_miri_zoom_in.transAxes, path_effects=pe)
            ax_miri_zoom_in.text(0.02, 0.86, miri_band_green.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='green',
                                     transform=ax_miri_zoom_in.transAxes, path_effects=pe)
            ax_miri_zoom_in.text(0.02, 0.80, miri_band_blue.upper(), horizontalalignment='left', verticalalignment='top',
                                     fontsize=fontsize_large, color='blue',
                                     transform=ax_miri_zoom_in.transAxes, path_effects=pe)
            helper_func.draw_box(ax=ax_miri_zoom_in, wcs=wcs_miri_zoom_in,
                                 coord=pos_obj,
                                 box_size=stamp_cutout_size, color='red', line_width=2, line_style='--')
            helper_func.plot_coord_circle(ax=ax_miri_zoom_in, pos=pos_obj,
                                          rad=circle_rad_region, color='white', line_style='-', line_width=2, alpha=1., fill=False)

        # add stamp axis
        stamp_row_count_down = 0
        # add stamp axis for hst
        ax_hst_stamp_list = []
        hst_stamp_index = 0

        for hst_stamp_index, hst_stamp_band in enumerate(hst_stamp_band_list):

            ax_hst_stamp_list.append(
                figure.add_axes([stamp_left_align + hst_stamp_index*(stamp_width + stamp_space_vertical),
                                 stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                 stamp_width, stamp_height],
                                projection=cutout_dict_stamp['%s_img_cutout' % hst_stamp_band].wcs))
            norm_hst_stamp = helper_func.compute_cbar_norm(
                cutout_list=cutout_dict_stamp['%s_img_cutout' % hst_stamp_band].data, log_scale=cbar_log_scale)

            ax_hst_stamp_list[hst_stamp_index].imshow(cutout_dict_stamp['%s_img_cutout' % hst_stamp_band].data,
                                            norm=norm_hst_stamp, cmap='Greys')
            ax_hst_stamp_list[hst_stamp_index].set_title(hst_stamp_band.upper(), fontsize=fontsize_large,
                                                         color='tab:blue')
            helper_func.plot_coord_circle(ax=ax_hst_stamp_list[hst_stamp_index],
                                          pos=pos_obj,
                                          rad=circle_rad_region, color='red',
                                          line_style='--', line_width=2, alpha=1., fill=False)
            helper_func.plot_coord_circle(ax=ax_hst_stamp_list[hst_stamp_index],
                                          pos=pos_obj,
                                          rad=circle_rad_obj, color='blue',
                                          line_style='-', line_width=2, alpha=1., fill=False)
            if hst_stamp_index == 0:
                ra_tick_label, dec_tick_label = True, True
            else:
                ra_tick_label, dec_tick_label = False, False
            self.arr_axis_params(ax=ax_hst_stamp_list[hst_stamp_index], ra_tick_label=ra_tick_label,
                                 dec_tick_label=dec_tick_label, ra_axis_label=' ', dec_axis_label=' ',
                                 fontsize=fontsize_small, labelsize=fontsize_small)

        if include_h_alpha:
            ax_hst_ha_stamp = figure.add_axes([stamp_left_align + (hst_stamp_index + 2)*(stamp_width + stamp_space_vertical),
                                 stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                 stamp_width, stamp_height],
                                projection=cutout_dict_stamp['%s_img_cutout' % hst_ha_stamp_band].wcs)
            norm_hst_ha_stamp = helper_func.compute_cbar_norm(
                cutout_list=cutout_dict_stamp['%s_img_cutout' % hst_ha_stamp_band].data, log_scale=cbar_log_scale)
            ax_hst_ha_stamp.imshow(cutout_dict_stamp['%s_img_cutout' % hst_ha_stamp_band].data,
                                   norm=norm_hst_ha_stamp, cmap='Greys')
            ax_hst_ha_stamp.set_title(r'H$\alpha$ ' + hst_ha_stamp_band.upper(), fontsize=fontsize_large, color='tab:red')
            helper_func.plot_coord_circle(ax=ax_hst_ha_stamp,
                                          pos=pos_obj,
                                          rad=circle_rad_region, color='red',
                                          line_style='--', line_width=2, alpha=1., fill=False)
            helper_func.plot_coord_circle(ax=ax_hst_ha_stamp,
                                          pos=pos_obj,
                                          rad=circle_rad_obj, color='blue',
                                          line_style='-', line_width=2, alpha=1., fill=False)
            self.arr_axis_params(ax=ax_hst_ha_stamp, ra_tick_label=True, dec_tick_label=True,
                        ra_axis_label=' ', dec_axis_label=' ',
                        ra_minpad=0.8, dec_minpad=0.8, tick_color='k', label_color='k',
                        fontsize=fontsize_small, labelsize=fontsize_small, ra_tick_num=3, dec_tick_num=3)

            ax_hst_ha_cont_sub_stamp = figure.add_axes([stamp_left_align + (hst_stamp_index + 4)*(stamp_width + stamp_space_vertical),
                                 stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                 stamp_width, stamp_height],
                                projection=cutout_dict_stamp['%s_img_cutout' % ha_cont_sub_stamp_band].wcs)
            norm_hst_ha_cont_sub_stamp = helper_func.compute_cbar_norm(
                cutout_list=cutout_dict_stamp['%s_img_cutout' % ha_cont_sub_stamp_band].data, log_scale=cbar_log_scale)
            ax_hst_ha_cont_sub_stamp.imshow(cutout_dict_stamp['%s_img_cutout' % ha_cont_sub_stamp_band].data,
                                            norm=norm_hst_ha_cont_sub_stamp, cmap='Greys')
            ax_hst_ha_cont_sub_stamp.set_title(r'H$\alpha$ cont. sub.', fontsize=fontsize_large, color='tab:red')
            helper_func.plot_coord_circle(ax=ax_hst_ha_cont_sub_stamp,
                                          pos=pos_obj,
                                          rad=circle_rad_region, color='red',
                                          line_style='--', line_width=2, alpha=1., fill=False)
            helper_func.plot_coord_circle(ax=ax_hst_ha_cont_sub_stamp,
                                          pos=pos_obj,
                                          rad=circle_rad_obj, color='blue',
                                          line_style='-', line_width=2, alpha=1., fill=False)
            self.arr_axis_params(ax=ax_hst_ha_cont_sub_stamp, ra_tick_label=True, dec_tick_label=True,
                        ra_axis_label=' ', dec_axis_label=' ',
                        ra_minpad=0.8, dec_minpad=0.8, tick_color='k', label_color='k',
                        fontsize=fontsize_small, labelsize=fontsize_small, ra_tick_num=3, dec_tick_num=3)


        nircam_stamp_index = 0
        if include_nircam:
            stamp_row_count_down -= 1
            ax_nircam_stamp_list = []
            for nircam_stamp_index, nircam_stamp_band in enumerate(nircam_stamp_band_list):
                if ((np.sum(cutout_dict_stamp['%s_img_cutout' % nircam_stamp_band].data) == 0) or
                        (cutout_dict_stamp['%s_img_cutout' % nircam_stamp_band].data is None)):
                    continue
                ax_nircam_stamp_list.append(
                    figure.add_axes([stamp_left_align + nircam_stamp_index*(stamp_width + stamp_space_vertical),
                                     stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                     stamp_width, stamp_height],
                                    projection=cutout_dict_stamp['%s_img_cutout' % nircam_stamp_band].wcs))
                norm_nircam_stamp = helper_func.compute_cbar_norm(
                cutout_list=cutout_dict_stamp['%s_img_cutout' % nircam_stamp_band].data, log_scale=cbar_log_scale)
                ax_nircam_stamp_list[nircam_stamp_index].imshow(cutout_dict_stamp['%s_img_cutout' % nircam_stamp_band].data,
                                                norm=norm_nircam_stamp, cmap='Greys')
                ax_nircam_stamp_list[nircam_stamp_index].set_title(nircam_stamp_band.upper(), fontsize=fontsize_large,
                                                                   color='tab:green')
                helper_func.plot_coord_circle(ax=ax_nircam_stamp_list[nircam_stamp_index],
                                          pos=pos_obj,
                                          rad=circle_rad_region, color='red',
                                          line_style='--', line_width=2, alpha=1., fill=False)
                helper_func.plot_coord_circle(ax=ax_nircam_stamp_list[nircam_stamp_index],
                                          pos=pos_obj,
                                          rad=circle_rad_obj, color='blue',
                                          line_style='-', line_width=2, alpha=1., fill=False)
                if nircam_stamp_index == 0:
                    ra_tick_label, dec_tick_label = True, True
                else:
                    ra_tick_label, dec_tick_label = False, False
                self.arr_axis_params(ax=ax_nircam_stamp_list[nircam_stamp_index], ra_tick_label=ra_tick_label,
                                     dec_tick_label=dec_tick_label, ra_axis_label=' ', dec_axis_label=' ',
                                     fontsize=fontsize_small, labelsize=fontsize_small)

        if include_miri:
            # stamp_row_count_down -= 1
            ax_miri_stamp_list = []
            for miri_stamp_index, miri_stamp_band in enumerate(miri_stamp_band_list):
                if ((np.sum(cutout_dict_stamp['%s_img_cutout' % miri_stamp_band].data) == 0) or
                        (cutout_dict_stamp['%s_img_cutout' % miri_stamp_band].data is None)):
                    continue
                ax_miri_stamp_list.append(
                    figure.add_axes([stamp_left_align + (nircam_stamp_index + miri_stamp_index + 2)*(stamp_width + stamp_space_vertical),
                                     stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                     stamp_width, stamp_height],
                                    projection=cutout_dict_stamp['%s_img_cutout' % miri_stamp_band].wcs))
                norm_miri_stamp = helper_func.compute_cbar_norm(
                cutout_list=cutout_dict_stamp['%s_img_cutout' % miri_stamp_band].data, log_scale=cbar_log_scale)
                ax_miri_stamp_list[miri_stamp_index].imshow(cutout_dict_stamp['%s_img_cutout' % miri_stamp_band].data,
                                                norm=norm_miri_stamp, cmap='Greys')
                ax_miri_stamp_list[miri_stamp_index].set_title(miri_stamp_band.upper(), fontsize=fontsize_large,
                                                               color='tab:purple')
                helper_func.plot_coord_circle(ax=ax_miri_stamp_list[miri_stamp_index],
                                              pos=pos_obj,
                                              rad=circle_rad_region, color='red',
                                              line_style='--', line_width=2, alpha=1., fill=False)
                helper_func.plot_coord_circle(ax=ax_miri_stamp_list[miri_stamp_index],
                                              pos=pos_obj,
                                              rad=circle_rad_obj, color='blue',
                                              line_style='-', line_width=2, alpha=1., fill=False)
                if miri_stamp_index == 0:
                    ra_tick_label, dec_tick_label = True, True
                else:
                    ra_tick_label, dec_tick_label = False, False
                self.arr_axis_params(ax=ax_miri_stamp_list[miri_stamp_index], ra_tick_label=ra_tick_label,
                                     dec_tick_label=dec_tick_label, ra_axis_label=' ', dec_axis_label=' ',
                                     fontsize=fontsize_small, labelsize=fontsize_small)

        # plot flux-density distribution
        ax_flux_density = figure.add_axes([0.065, 0.27, 0.9, 0.15])

        for band in band_list:
            if 'ha' in band:
                continue
            if 'flux_%s' % band not in flux_dict:
                continue

            if band in self.get_hst_band_list():
                color = 'tab:blue'
            elif band in self.get_hst_ha_band_list():
                color = 'tab:red'
            elif band in self.get_nircam_band_list():
                color = 'tab:green'
            elif band in self.get_miri_band_list():
                color = 'tab:purple'
            elif band in self.get_astrosat_band_list():
                color = 'tab:pink'
            else:
                color = ''

            mean_wave = self.get_band_wave(band=band)
            min_wave = self.get_band_wave(band=band, wave_estimator='min_wave')
            max_wave = self.get_band_wave(band=band, wave_estimator='max_wave')

            ax_flux_density.errorbar(mean_wave, flux_dict['flux_%s' % band],
                                     xerr=[[mean_wave-min_wave], [max_wave-mean_wave]],
                                     yerr=flux_dict['flux_err_%s' % band],
                                     fmt='o', color=color, ms=20)
            ax_flux_density.errorbar(mean_wave, flux_dict_small_target['flux_%s' % band],
                                     xerr=[[mean_wave-min_wave], [max_wave-mean_wave]],
                                     yerr=flux_dict_small_target['flux_err_%s' % band],
                                     fmt='o', color='k', ms=20)

        ax_flux_density.set_yscale('log')
        ax_flux_density.set_xscale('log')
        ax_flux_density.set_xlabel(r'Wavelength [$\mu$m]', fontsize=fontsize_large)
        ax_flux_density.set_ylabel(r'Flux density [erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$]', fontsize=fontsize_large)
        ax_flux_density.tick_params(axis='both', which='both', width=2, direction='in', labelsize=fontsize_large)

        ax_flux_density.set_title('NUV flux_ratio: %.2f' %
                                  (flux_dict_small_target['flux_F275W'] / flux_dict['flux_F275W']),
                                  fontsize=fontsize_large)

        if astrosat_band is not None:

            # show astrosat image
            ax_astrosat = figure.add_axes([0.09, 0.40, 0.15, 0.15],
                                          projection=cutout_dict_zoom_in['%s_img_cutout' % astrosat_band].wcs)
            ax_astrosat_cbar = figure.add_axes([0.245, 0.43, 0.015, 0.09])
            norm_astrosat = helper_func.compute_cbar_norm(
                    cutout_list=cutout_dict_stamp['%s_img_cutout' % astrosat_band].data, log_scale=cbar_log_scale)
            ax_astrosat.imshow(cutout_dict_zoom_in['%s_img_cutout' % astrosat_band].data, norm=norm_astrosat, cmap='bone_r')
            ax_astrosat.set_title(astrosat_band.upper(), fontsize=fontsize_large, color='tab:pink')
            helper_func.plot_coord_circle(ax=ax_astrosat,
                                              pos=pos_obj,
                                              rad=circle_rad_region, color='red',
                                              line_style='--', line_width=2, alpha=1., fill=False)
            self.arr_axis_params(ax=ax_astrosat, ra_axis_label=' ', dec_axis_label=' ',
                                 fontsize=fontsize_small, labelsize=fontsize_small)
            helper_func.create_cbar(ax_cbar=ax_astrosat_cbar, cmap='bone_r', norm=norm_astrosat,
                                    cbar_label=r'S [erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$]',
                                    fontsize=fontsize_small, ticks=None, labelpad=2, tick_width=2,
                                    orientation='vertical', extend='neither')

        # hst_flux_title_str = ('id_phangs_cluster = %i, '
        #                   'HUM_CLASS = %s, '
        #                   'R.A.=%.4f, DEC.=%.4f,    '
        #                       'age= %i, E(B-V)=%.2f, log(M$_*$/M$_{\odot}$)=%.2f '
        #                   % (cluster_dict['phangs_cluster_id'],
        #                      cluster_dict['cluster_class'],
        #                      cluster_dict['ra'],
        #                      cluster_dict['dec'],
        #                      cluster_dict['age'],
        #                      cluster_dict['ebv'],
        #                      np.log10(cluster_dict['mass'])))
        #
        # ax_hst_spec.set_title(hst_flux_title_str, fontsize=fontsize_large)
        # ax_hst_spec.tick_params(axis='both', which='both', width=2, direction='in', labelsize=fontsize_large)

        # add axis
        # get MUSE H-alpha cutout
        if muse_fit_dict is not None:
            muse_ha_cutout = helper_func.get_img_cutout(img=muse_fit_dict['muse_ha_map_dict']['muse_map_data'],
                                                        wcs=muse_fit_dict['muse_ha_map_dict']['muse_map_wcs'],
                                                        coord=pos_obj, cutout_size=env_cutout_size)

            ax_muse_spec = figure.add_axes([0.065, 0.025, 0.90, 0.2])
            ax_muse_zoom_in = figure.add_axes([0.10, 0.13, 0.18, 0.18], projection=muse_ha_cutout.wcs)
            ax_red_bump_zoom_in = figure.add_axes([0.50, 0.07, 0.30, 0.10])

            # plot spectrum
            ax_muse_spec.plot(muse_fit_dict['wavelength'], muse_fit_dict['total_flux'] * 1e16, color='k',
                              label='Obs. Spectrum')
            ax_muse_spec.plot(muse_fit_dict['wavelength'], muse_fit_dict['best_fit'] * 1e16, color='tab:red',
                              label='Best total fit')
            ax_muse_spec.plot(muse_fit_dict['wavelength'], muse_fit_dict['continuum_best_fit'] * 1e16, color='tab:orange',
                              label='Best Continuum fit')
            muse_title_str = (r'age= %.2f Myr, M/H=%.3f, Mass-to-light=%.2f, A$_v (star)=$%.2f mag' %
                              (10 ** (muse_fit_dict['ages']) * 1e-6,
                               muse_fit_dict['met'],
                               muse_fit_dict['mass2light'],
                               muse_fit_dict['star_red'])
                              + '\n' +
                              'A$_v (gas)=$%.2f mag, '
                              r'EW(H$\alpha$)=%.1f$\AA$, EW(H$\beta$)=%.1f$\AA$, 12 + log (O/H) = %.2f'
                              % (muse_fit_dict['gas_red'],
                                 muse_fit_dict['ha_ew'],
                                 muse_fit_dict['hb_ew'],
                                 muse_fit_dict['gas_phase_met']))
            ax_muse_spec.set_title(muse_title_str, fontsize=fontsize_large, loc='right')
            ax_muse_spec.tick_params(axis='both', which='both', width=2, direction='in', labelsize=fontsize_large)
            ax_muse_spec.set_xlabel(r'Wavelength [$\AA$]', fontsize=fontsize_large)
            ax_muse_spec.set_ylabel(r'Flux density [10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$]', fontsize=fontsize_large)
            ax_muse_spec.legend(frameon=False, fontsize=fontsize_large)
            ax_muse_spec.set_xlim(np.nanmin(muse_fit_dict['wavelength']), np.nanmax(muse_fit_dict['wavelength']))

            # plot ha plot cutout
            norm_muse_ha = helper_func.compute_cbar_norm(cutout_list=muse_ha_cutout.data, log_scale=cbar_log_scale)
            ax_muse_zoom_in.imshow(muse_ha_cutout.data, norm=norm_muse_ha, cmap='cividis')
            ax_muse_zoom_in.set_title(r'MUSE H$\alpha$', fontsize=fontsize_large, color='k')
            helper_func.plot_coord_circle(ax=ax_muse_zoom_in, pos=pos_obj, rad=circle_rad_region, color='red',
                                          line_style='--', line_width=2, alpha=1., fill=False)
            self.arr_axis_params(ax=ax_muse_zoom_in, ra_axis_label=' ', dec_axis_label=' ',
                                 fontsize=fontsize_small, labelsize=fontsize_small)

            # plot red bump
            min_wave_red_bump = 5550
            max_wave_red_bump = 6050
            ax_red_bump_zoom_in.plot(muse_fit_dict['wavelength'], muse_fit_dict['total_flux'] * 1e16, color='k')
            ax_red_bump_zoom_in.plot(muse_fit_dict['wavelength'], muse_fit_dict['best_fit'] * 1e16, color='tab:orange')
            ax_red_bump_zoom_in.plot(muse_fit_dict['wavelength'], muse_fit_dict['continuum_best_fit'] * 1e16, color='tab:orange')
            ax_red_bump_zoom_in.tick_params(axis='both', which='both', width=2, direction='in', labelsize=fontsize_large)
            ax_red_bump_zoom_in.set_xlim(min_wave_red_bump, max_wave_red_bump)
            min_y_value_red_bump = np.nanmin((muse_fit_dict['total_flux'] * 1e16)[(muse_fit_dict['wavelength'] > min_wave_red_bump) &
                                                                                  (muse_fit_dict['wavelength'] < max_wave_red_bump) &
                                                                                  (muse_fit_dict['total_flux'] > 0)])
            max_y_value_red_bump = np.nanmax((muse_fit_dict['total_flux'] * 1e16)[(muse_fit_dict['wavelength'] > min_wave_red_bump) &
                                             (muse_fit_dict['wavelength'] < max_wave_red_bump)])
            add_space = (max_y_value_red_bump - min_y_value_red_bump) / 10
            ax_red_bump_zoom_in.set_ylim(min_y_value_red_bump - add_space, max_y_value_red_bump + add_space)
            ax_red_bump_zoom_in.set_title('Red bump', fontsize=fontsize_large)

        return figure

    def get_rgb_zoom_in(self, ra, dec, cutout_size, circle_rad, band_red, band_green, band_blue, ref_band='blue'):
        """
        Function to create an RGB image of a zoom in region of PHANGS observations

        Parameters
        ----------
        ra, dec : float
            coordinates in degree
        cutout_size: tuple
            cutout size in arcsec
        circle_rad : float
            radius of circle in which the object of interest is located
        band_red, band_green, band_blue : str
            filter names
        ref_band : str
            can be blue, green or red. In case the images are not the same size they get reprojected to one frame

        Returns
        -------
        rgb_image, wcs : ``numpy.ndarray``, ``astropy.wcs.WCS``
        """

        self.load_phangs_bands(band_list=[band_red, band_green, band_blue],
                               flux_unit='MJy/sr', load_err=False)

        cutout = self.get_band_cutout_dict(ra_cutout=ra, dec_cutout=dec, cutout_size=cutout_size,
                                           band_list=[band_red, band_green, band_blue])



        if ((cutout['%s_img_cutout' % band_red].data is None) or
                (cutout['%s_img_cutout' % band_green].data is None) or
                (cutout['%s_img_cutout' % band_blue].data is None)):
            return None, None

        ref_wcs = cutout['%s_img_cutout' % eval('band_%s' % ref_band)].wcs

        if not (cutout['%s_img_cutout' % band_red].data.shape ==
                cutout['%s_img_cutout' % band_green].data.shape ==
                cutout['%s_img_cutout' % band_blue].data.shape):
            new_shape = cutout['%s_img_cutout' % eval('band_%s' % ref_band)].data.shape
            if ref_band == 'red':
                cutout_data_red = cutout['%s_img_cutout' % band_red].data
            else:
                cutout_data_red = helper_func.reproject_image(data=cutout['%s_img_cutout' % band_red].data,
                                                              wcs=cutout['%s_img_cutout' % band_red].wcs,
                                                              new_wcs=ref_wcs, new_shape=new_shape)
            if ref_band == 'green':
                cutout_data_green = cutout['%s_img_cutout' % band_green].data
            else:
                cutout_data_green = helper_func.reproject_image(data=cutout['%s_img_cutout' % band_green].data,
                                                                wcs=cutout['%s_img_cutout' % band_green].wcs,
                                                                new_wcs=ref_wcs, new_shape=new_shape)
            if ref_band == 'blue':
                cutout_data_blue = cutout['%s_img_cutout' % band_blue].data
            else:
                cutout_data_blue = helper_func.reproject_image(data=cutout['%s_img_cutout' % band_blue].data,
                                                                wcs=cutout['%s_img_cutout' % band_blue].wcs,
                                                                new_wcs=ref_wcs, new_shape=new_shape)
        else:
            cutout_data_red = cutout['%s_img_cutout' % band_red].data
            cutout_data_green = cutout['%s_img_cutout' % band_green].data
            cutout_data_blue = cutout['%s_img_cutout' % band_blue].data

        # get rgb image
        min_red_hst, max_red_hst = self.get_image_scale_with_circle(
            img_data=cutout_data_red,
            img_wcs=ref_wcs, ra=ra, dec=dec, circle_rad=circle_rad)
        min_green_hst, max_green_hst = self.get_image_scale_with_circle(
            img_data=cutout_data_green,
            img_wcs=ref_wcs, ra=ra, dec=dec, circle_rad=circle_rad)
        min_blue_hst, max_blue_hst = self.get_image_scale_with_circle(
            img_data=cutout_data_blue,
            img_wcs=ref_wcs, ra=ra, dec=dec, circle_rad=circle_rad)

        cutout_rgb_img = self.get_rgb_img(data_r=cutout_data_red,
                                   data_g=cutout_data_green,
                                   data_b=cutout_data_blue,
                                   min_max_r=(min_red_hst, max_red_hst),
                                   min_max_g=(min_green_hst, max_green_hst),
                                   min_max_b=(min_blue_hst, max_blue_hst),
                                   scaletype_r='abs', scaletype_g='abs', scaletype_b='abs')

        return cutout_rgb_img, ref_wcs

    def get_rgb_img_(self, ra, dec, cutout_size, band_red, band_green, band_blue, ref_band='blue',
                     color_r='#FF4433', color_g='#40E0D0', color_b='#1F51FF',
                     overview_img_size=(500, 500),
                     min_max_r=(0.3, 99.5), min_max_g=(0.3, 99.5), min_max_b=(0.3, 99.5),
                                   gamma_r=17.5, gamma_g=17.5, gamma_b=17.5,
                                   gamma_corr_r=17.5, gamma_corr_g=17.5, gamma_corr_b=17.5,
                                   combined_gamma=17.5):
        """
        Function to create an RGB image of a zoom in region of PHANGS observations

        Parameters
        ----------
        ra, dec : float
            coordinates in degree
        cutout_size: tuple
            cutout size in arcsec
        circle_rad : float
            radius of circle in which the object of interest is located
        band_red, band_green, band_blue : str
            filter names
        ref_band : str
            can be blue, green or red. In case the images are not the same size they get reprojected to one frame

        Returns
        -------
        rgb_image, wcs : ``numpy.ndarray``, ``astropy.wcs.WCS``
        """

        self.load_phangs_bands(band_list=[band_red, band_green, band_blue],
                               flux_unit='MJy/sr', load_err=False)

        cutout = self.get_band_cutout_dict(ra_cutout=ra, dec_cutout=dec, cutout_size=cutout_size,
                                           band_list=[band_red, band_green, band_blue])

        ref_cutout_shape = cutout['%s_img_cutout' % eval('band_%s' % ref_band)].data.shape

        ra_max = cutout['%s_img_cutout' % eval('band_%s' % ref_band)].wcs.pixel_to_world(
            0, 0).ra.value
        ra_min = cutout['%s_img_cutout' % eval('band_%s' % ref_band)].wcs.pixel_to_world(
            ref_cutout_shape[0], ref_cutout_shape[1]).ra.value
        dec_min = cutout['%s_img_cutout' % eval('band_%s' % ref_band)].wcs.pixel_to_world(
            0, 0).dec.value
        dec_max = cutout['%s_img_cutout' % eval('band_%s' % ref_band)].wcs.pixel_to_world(
            ref_cutout_shape[0], ref_cutout_shape[1]).dec.value

        new_wcs = helper_func.construct_wcs(ra_min=ra_min, ra_max=ra_max, dec_min=dec_max, dec_max=dec_min,
                                            img_shape=overview_img_size, quadratic_image=False)


        img_data_red = helper_func.reproject_image(data=cutout['%s_img_cutout' % band_red].data,
                                                   wcs=cutout['%s_img_cutout' % band_red].wcs,
                                                   new_wcs=new_wcs, new_shape=overview_img_size)
        img_data_green = helper_func.reproject_image(data=cutout['%s_img_cutout' % band_green].data,
                                                     wcs=cutout['%s_img_cutout' % band_green].wcs,
                                                     new_wcs=new_wcs, new_shape=overview_img_size)
        img_data_blue = helper_func.reproject_image(data=cutout['%s_img_cutout' % band_blue].data,
                                                    wcs=cutout['%s_img_cutout' % band_blue].wcs,
                                                    new_wcs=new_wcs, new_shape=overview_img_size)

        img_data_red[img_data_red == 0] = np.nan
        img_data_green[img_data_green == 0] = np.nan
        img_data_blue[img_data_blue == 0] = np.nan



        rgb_img = self.get_rgb_img(data_r=img_data_red, data_g=img_data_green, data_b=img_data_blue,
                                   color_r=color_r, color_g=color_g, color_b=color_b,
                                   min_max_r=min_max_r, min_max_g=min_max_g, min_max_b=min_max_b,
                                   gamma_r=gamma_r, gamma_g=gamma_g, gamma_b=gamma_b,
                                   gamma_corr_r=gamma_corr_r, gamma_corr_g=gamma_corr_g, gamma_corr_b=gamma_corr_b,
                                   combined_gamma=combined_gamma)

        return rgb_img, new_wcs

    def get_target_overview_rgb_img(self, red_band, green_band, blue_band, ref_band='red',
                                    overview_img_size=(500, 500)):
        """
        Function to create an overview RGB image of PHANGS HST observations

        Parameters
        ----------
        red_band, green_band, blue_band : str
            Can be specified to any hst band
        ref_band: str
            Band which is used for the image limits. can be red green or blue
        overview_img_size : tuple
            denotes the shape of the new image

        Returns
        -------
        rgb_image, wcs : ``numpy.ndarray``, ``astropy.wcs.WCS``
        """

        # band list need to be loaded
        self.load_phangs_bands(band_list=[red_band, green_band, blue_band], flux_unit='MJy/sr', load_err=False)
        # get overview image
        non_zero_elements = np.where(self.hst_bands_data['%s_data_img' % eval('%s_band' % ref_band)] != 0)

        min_index_ra_axis_x_val = int(np.mean(non_zero_elements[1][non_zero_elements[1] ==
                                                                   np.min(non_zero_elements[1])]))
        min_index_ra_axis_y_val = int(np.mean(non_zero_elements[0][non_zero_elements[1] ==
                                                                   np.min(non_zero_elements[1])]))
        max_index_ra_axis_x_val = int(np.mean(non_zero_elements[1][non_zero_elements[1] ==
                                                                   np.max(non_zero_elements[1])]))
        max_index_ra_axis_y_val = int(np.mean(non_zero_elements[0][non_zero_elements[1] ==
                                                                   np.max(non_zero_elements[1])]))
        min_index_dec_axis_x_val = int(np.mean(non_zero_elements[1][non_zero_elements[0] ==
                                                                    np.min(non_zero_elements[0])]))
        min_index_dec_axis_y_val = int(np.mean(non_zero_elements[0][non_zero_elements[0] ==
                                                                    np.min(non_zero_elements[0])]))
        max_index_dec_axis_x_val = int(np.mean(non_zero_elements[1][non_zero_elements[0] ==
                                                                    np.max(non_zero_elements[0])]))
        max_index_dec_axis_y_val = int(np.mean(non_zero_elements[0][non_zero_elements[0] ==
                                                                    np.max(non_zero_elements[0])]))

        ra_max = self.hst_bands_data['%s_wcs_img' % eval('%s_band' % ref_band)].pixel_to_world(
            min_index_ra_axis_x_val, min_index_ra_axis_y_val).ra.value
        ra_min = self.hst_bands_data['%s_wcs_img' % eval('%s_band' % ref_band)].pixel_to_world(
            max_index_ra_axis_x_val, max_index_ra_axis_y_val).ra.value

        dec_min = self.hst_bands_data['%s_wcs_img' % eval('%s_band' % ref_band)].pixel_to_world(
            min_index_dec_axis_x_val, min_index_dec_axis_y_val).dec.value
        dec_max = self.hst_bands_data['%s_wcs_img' % eval('%s_band' % ref_band)].pixel_to_world(
            max_index_dec_axis_x_val, max_index_dec_axis_y_val).dec.value

        new_wcs = helper_func.construct_wcs(ra_min=ra_min, ra_max=ra_max, dec_min=dec_max, dec_max=dec_min,
                                            img_shape=overview_img_size, quadratic_image=True)

        img_data_red = helper_func.reproject_image(data=self.hst_bands_data['%s_data_img' % red_band],
                                                   wcs=self.hst_bands_data['%s_wcs_img' % red_band],
                                                   new_wcs=new_wcs, new_shape=overview_img_size)
        img_data_green = helper_func.reproject_image(data=self.hst_bands_data['%s_data_img' % green_band],
                                                     wcs=self.hst_bands_data['%s_wcs_img' % green_band],
                                                     new_wcs=new_wcs, new_shape=overview_img_size)
        img_data_blue = helper_func.reproject_image(data=self.hst_bands_data['%s_data_img' % blue_band],
                                                    wcs=self.hst_bands_data['%s_wcs_img' % blue_band],
                                                    new_wcs=new_wcs, new_shape=overview_img_size)

        img_data_red[img_data_red == 0] = np.nan
        img_data_green[img_data_green == 0] = np.nan
        img_data_blue[img_data_blue == 0] = np.nan

        hst_rgb = self.get_rgb_img(data_r=img_data_red, data_g=img_data_green, data_b=img_data_blue,
                                   min_max_r=(0.3, 99.5), min_max_g=(0.3, 99.5), min_max_b=(0.3, 99.5),
                                   gamma_r=17.5, gamma_g=17.5, gamma_b=17.5,
                                   gamma_corr_r=17.5, gamma_corr_g=17.5, gamma_corr_b=17.5,
                                   combined_gamma=17.5)

        return hst_rgb, new_wcs

    @staticmethod
    def get_image_scale_with_circle(img_data, img_wcs, ra, dec, circle_rad=0.16, box_scaling=1/5, filter_scaling=1/10):
        """
        Function calculate scaling in image by taking median background of image and the maximum value inside a circle

        Parameters
        ----------
        img_data : ``numpy.ndarray``
            data image
        img_wcs : ``astropy.wcs.WCS``
            world coordinate system
        ra, dec : float
            coordinates
        circle_rad : float
            radius of circle
        box_scaling : float
            factor by how much the box estimation for the Background2D should be relative to the input image
        filter_scaling : float
            factor by how much the filter for the Background2D should be relative to the input image


        Returns
        -------
        (bkg_median, max_value) : tuple
            median background and maximum value inside the circle
        """
        # get background value as minimum
        sigma_clip = SigmaClip()
        bkg_estimator = MedianBackground()
        box_size = list([int(img_data.shape[0]*box_scaling), int(img_data.shape[0]*box_scaling)])
        filter_size = list([int(img_data.shape[0]*filter_scaling), int(img_data.shape[0]*filter_scaling)])
        # assure that filter has an odd size
        if filter_size[0] % 2 == 0:
            filter_size[0] += 1
        if filter_size[1] % 2 == 0:
            filter_size[1] += 1

        # get background estimation
        bkg = Background2D(img_data, box_size=box_size, filter_size=filter_size,
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        bkg_median = bkg.background_median

        # get coordinates and radius in pixel scale
        central_pos_world = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        central_pos_pixel = img_wcs.world_to_pixel(central_pos_world)
        circle_rad_pix = helper_func.transform_world2pix_scale(length_in_arcsec=circle_rad, wcs=img_wcs)
        # get meshgrid of image
        mesh_x, mesh_y = np.meshgrid(np.linspace(0, img_data.shape[1]-1, img_data.shape[1]),
                                     np.linspace(0, img_data.shape[0]-1, img_data.shape[0]))
        # mask pixels inside the radius
        mask_inside_circle = np.sqrt((mesh_x - central_pos_pixel[0]) ** 2 +
                                     (mesh_y - central_pos_pixel[1])**2) < circle_rad_pix
        max_value = np.nanmax(img_data[mask_inside_circle])

        # what if the background is higher than the maximal value in the circle ?
        if bkg_median > max_value:
            return max_value/10, max_value
        else:
            return bkg_median, max_value

    @staticmethod
    def get_rgb_img(data_r, data_g, data_b, color_r='#FF4433', color_g='#0FFF50', color_b='#1F51FF',
                    min_max_r=None, min_max_g=None, min_max_b=None,
                    rescalefn='asinh',
                    scaletype_r='perc', scaletype_g='perc', scaletype_b='perc',
                    gamma_r=2.2, gamma_g=2.2, gamma_b=2.2,
                    gamma_corr_r=2.2, gamma_corr_g=2.2, gamma_corr_b=2.2, combined_gamma=2.2):
        """
        Function to create an RGB image

        Parameters
        ----------
        data_r, data_g, data_b : ``numpy.ndarray``, array
            color images. Must be all same shape
        color_r, color_g, color_b: str
            hex code for color
        min_max_r, min_max_g, min_max_b : tuple or None
            denotes the percentages till where the data is used
        rescalefn : str
            scale function can be linear sqrt squared log power sinh asinh
        scaletype_r, scaletype_g, scaletype_b : str
        'abs' for absolute values, 'perc' for percentiles
        gamma_r, gamma_g, gamma_b : float
            gamma factor for each individual color band
        gamma_corr_r, gamma_corr_g, gamma_corr_b : float
            gamma correction factor for each grey scale image
        combined_gamma : float
            gamma factor of resulting rgb image

        Returns
        -------
        rgb image : ``numpy.ndarray``
            of shape (N,N, 3)
        """
        if min_max_r is None:
            min_max_r = [0., 100.]
        if min_max_g is None:
            min_max_g = [0., 100.]
        if min_max_b is None:
            min_max_b = [0., 100.]

        grey_r = mcf.greyRGBize_image(data_r, rescalefn=rescalefn, scaletype=scaletype_r, min_max=min_max_r,
                                      gamma=gamma_r)
        grey_g = mcf.greyRGBize_image(data_g, rescalefn=rescalefn, scaletype=scaletype_g, min_max=min_max_g,
                                      gamma=gamma_g)
        grey_b = mcf.greyRGBize_image(data_b, rescalefn=rescalefn, scaletype=scaletype_b, min_max=min_max_b,
                                      gamma=gamma_b)
        r = mcf.colorize_image(grey_r, color_r, colorintype='hex', gammacorr_color=gamma_corr_r)
        g = mcf.colorize_image(grey_g, color_g, colorintype='hex', gammacorr_color=gamma_corr_g)
        b = mcf.colorize_image(grey_b, color_b, colorintype='hex', gammacorr_color=gamma_corr_b)
        return mcf.combine_multicolor([r, g, b], gamma=combined_gamma)

    @staticmethod
    def plot_circle_on_wcs_img(ax, pos, rad, color, line_style='-', line_width=2., alpha=1., fill=False):
        """
        plots circle on image using coordinates and WCS to orientate

        Parameters
        ----------
        ax : ``astropy.visualization.wcsaxes.core.WCSAxes``
            axis for plotting
        pos : ``astropy.coordinates.SkyCoord``
            position in form of Sky coordinates
        rad : float
            radius in arc seconds of circle
        color : str
            matplotlib color
        line_style : str
            matplotlib line style
        line_width : float
        alpha : float
            matplotlib alpha factor
        fill : bool
            flag whether circle should be filled or not

        Returns
        -------
        None
        """

        if fill:
            face_color = color
        else:
            face_color = 'none'

        if isinstance(pos, list):
            if not isinstance(rad, list):
                rad = [rad] * len(pos)
            if not isinstance(color, list):
                color = [color] * len(pos)
            if not isinstance(line_style, list):
                line_style = [line_style] * len(pos)
            if not isinstance(line_width, list):
                line_width = [line_width] * len(pos)
            if not isinstance(alpha, list):
                alpha = [alpha] * len(pos)
            for pos_i, rad_i, color_i, line_style_i, line_width_i, alpha_i in zip(pos, rad, color, line_style,
                                                                                  line_width, alpha):
                circle = SphericalCircle(pos_i, rad_i * u.arcsec, edgecolor=color_i, facecolor=face_color,
                                         linewidth=line_width_i,
                                         linestyle=line_style_i, alpha=alpha_i, transform=ax.get_transform('fk5'))
                ax.add_patch(circle)
        else:
            circle = SphericalCircle(pos, rad * u.arcsec, edgecolor=color, facecolor=face_color, linewidth=line_width,
                                     linestyle=line_style, alpha=alpha, transform=ax.get_transform('fk5'))
            ax.add_patch(circle)

    @staticmethod
    def arr_axis_params(ax, ra_tick_label=True, dec_tick_label=True,
                        ra_axis_label='R.A. (2000.0)', dec_axis_label='DEC. (2000.0)',
                        ra_minpad=0.8, dec_minpad=0.8, tick_color='k', label_color='k',
                        fontsize=15., labelsize=14., ra_tick_num=3, dec_tick_num=3):
        """
        plots circle on image using coordinates and WCS to orientate

        Parameters
        ----------
        ax : ``astropy.visualization.wcsaxes.core.WCSAxes``
            axis for plotting
        ra_tick_label, dec_tick_label : bool
        ra_axis_label, dec_axis_label : str
        ra_minpad, dec_minpad : float
        tick_color, label_color : str
        fontsize, labelsize : float
        ra_tick_num, dec_tick_num : int

        Returns
        -------
        None
        """
        ax.tick_params(which='both', width=1.5, length=7, direction='in', color=tick_color, labelsize=labelsize)

        if not ra_tick_label:
            ax.coords['ra'].set_ticklabel_visible(False)
            ax.coords['ra'].set_axislabel(' ', color=label_color)
        else:
            ax.coords['ra'].set_ticklabel(rotation=0, color=label_color)
            ax.coords['ra'].set_axislabel(ra_axis_label, minpad=ra_minpad, color=label_color, fontsize=fontsize)

        if not dec_tick_label:
            ax.coords['dec'].set_ticklabel_visible(False)
            ax.coords['dec'].set_axislabel(' ', color=label_color)
        else:
            ax.coords['dec'].set_ticklabel(rotation=90, color=label_color)
            ax.coords['dec'].set_axislabel(dec_axis_label, minpad=dec_minpad, color=label_color, fontsize=fontsize)

        ax.coords['ra'].set_ticks(number=ra_tick_num)
        ax.coords['ra'].display_minor_ticks(True)
        ax.coords['dec'].set_ticks(number=dec_tick_num)
        ax.coords['dec'].display_minor_ticks(True)
