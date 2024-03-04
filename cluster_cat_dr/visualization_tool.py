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

    def plot_ultimate_zoom_in_panel(self, ra_region, dec_region, include_h_alpha=True, ha_cont_sub_stamp_band='ha',
                                    include_nircam=True,
                                    include_miri=True,
                                    muse_fit_dict=None, muse_spec_dict=None,
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
        self.load_phangs_bands(band_list=band_list, flux_unit='MJy/sr',)

        # get hst_bvi_zoom_in
        img_hst_bvi_overview, wcs_hst_bvi_overview = self.get_target_hst_overview_rgb_img(red_band=hst_bvi_band_red,
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
        # get scales for cutout stamps
        norm_hst_stamps = helper_func.compute_cbar_norm(cutout_list=[cutout_dict_stamp['%s_img_cutout' % band].data
                                                                     for band in hst_stamp_band_list],
                                                        log_scale=cbar_log_scale)

        if include_h_alpha:
            norm_ha_native_stamp = helper_func.compute_cbar_norm(cutout_list=cutout_dict_stamp['%s_img_cutout' % hst_ha_stamp_band].data,
                                                               log_scale=cbar_log_scale)
            norm_ha_cont_sub_stamp = helper_func.compute_cbar_norm(cutout_list=cutout_dict_stamp['%s_img_cutout' % ha_cont_sub_stamp_band].data,
                                                               log_scale=cbar_log_scale)
        else:
            norm_ha_native_stamp = None
            norm_ha_cont_sub_stamp = None

        if include_nircam:
            norm_nircam_stamps = helper_func.compute_cbar_norm(cutout_list=[cutout_dict_stamp['%s_img_cutout' % band].data
                                                                            for band in nircam_stamp_band_list],
                                                               log_scale=cbar_log_scale)
        else:
            norm_nircam_stamps = None
        if include_nircam:
            norm_miri_stamps = helper_func.compute_cbar_norm(cutout_list=[cutout_dict_stamp['%s_img_cutout' % band].data
                                                                          for band in miri_stamp_band_list],
                                                             log_scale=cbar_log_scale)
        else:
            norm_miri_stamp = None

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
        stamp_bottom_align = 0.58
        stamp_space_horizontal = 0.005
        stamp_space_vertical = 0.005

        cbar_width = stamp_width * 2
        cbar_width_hst_ha = stamp_width
        cbar_height = 0.008
        c_bar_left_offset = 0.04
        c_bar_left_offset_hst_ha = 0.0
        c_bar_top_offset = -0.015

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
        if cbar_log_scale:
            cbar_label = r'log(S /[MJy / sr])'
        else:
            cbar_label = r'S [MJy / sr]'

        for hst_stamp_index, hst_stamp_band in enumerate(hst_stamp_band_list):

            ax_hst_stamp_list.append(
                figure.add_axes([stamp_left_align + hst_stamp_index*(stamp_width + stamp_space_vertical),
                                 stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                 stamp_width, stamp_height],
                                projection=cutout_dict_stamp['%s_img_cutout' % hst_stamp_band].wcs))
            ax_hst_stamp_list[hst_stamp_index].imshow(cutout_dict_stamp['%s_img_cutout' % hst_stamp_band].data,
                                            norm=norm_hst_stamps, cmap='Greys')
            ax_hst_stamp_list[hst_stamp_index].set_title(hst_stamp_band.upper(), fontsize=fontsize_large)
            if hst_stamp_index == 0:
                ra_tick_label, dec_tick_label = True, True
            else:
                ra_tick_label, dec_tick_label = False, False
            self.arr_axis_params(ax=ax_hst_stamp_list[hst_stamp_index], ra_tick_label=ra_tick_label,
                                 dec_tick_label=dec_tick_label, ra_axis_label=' ', dec_axis_label=' ',
                                 fontsize=fontsize_small, labelsize=fontsize_small)

        ax_hst_cbar = figure.add_axes([stamp_left_align + c_bar_left_offset,
                                       stamp_bottom_align + (stamp_row_count_down+1)*(stamp_height + stamp_space_horizontal) + c_bar_top_offset,
                                       cbar_width, cbar_height])
        helper_func.create_cbar(ax_cbar=ax_hst_cbar, cmap='Greys', norm=norm_hst_stamps, cbar_label=cbar_label,
                                fontsize=fontsize_large, ticks=None, labelpad=2, tick_width=2, orientation='horizontal',
                                extend='neither')

        if include_h_alpha:
            ax_hst_ha_stamp = figure.add_axes([stamp_left_align + (hst_stamp_index + 2)*(stamp_width + stamp_space_vertical),
                                 stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                 stamp_width, stamp_height],
                                projection=cutout_dict_stamp['%s_img_cutout' % hst_ha_stamp_band].wcs)
            ax_hst_ha_stamp.imshow(cutout_dict_stamp['%s_img_cutout' % hst_ha_stamp_band].data,
                                   norm=norm_ha_native_stamp, cmap='Greys')
            ax_hst_ha_stamp.set_title(hst_ha_stamp_band.upper(), fontsize=fontsize_large)
            self.arr_axis_params(ax=ax_hst_ha_stamp, ra_tick_label=True, dec_tick_label=True,
                        ra_axis_label=' ', dec_axis_label=' ',
                        ra_minpad=0.8, dec_minpad=0.8, tick_color='k', label_color='k',
                        fontsize=fontsize_small, labelsize=fontsize_small, ra_tick_num=3, dec_tick_num=3)
            ax_hst_ha_cbar = figure.add_axes([stamp_left_align + (hst_stamp_index + 2)*(stamp_width + stamp_space_vertical) + c_bar_left_offset_hst_ha,
                                              stamp_bottom_align + (stamp_row_count_down+1)*(stamp_height + stamp_space_horizontal) + c_bar_top_offset,
                                              cbar_width_hst_ha, cbar_height])
            print(norm_ha_native_stamp.vmin)
            print(norm_ha_native_stamp.vmax)
            helper_func.create_cbar(ax_cbar=ax_hst_ha_cbar, cmap='Greys', norm=norm_ha_native_stamp,
                                    cbar_label=cbar_label, fontsize=fontsize_large, orientation='horizontal',
                                    extend='neither')


            ax_hst_ha_cont_sub_stamp = figure.add_axes([stamp_left_align + (hst_stamp_index + 4)*(stamp_width + stamp_space_vertical),
                                 stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                 stamp_width, stamp_height],
                                projection=cutout_dict_stamp['%s_img_cutout' % ha_cont_sub_stamp_band].wcs)
            ax_hst_ha_cont_sub_stamp.imshow(cutout_dict_stamp['%s_img_cutout' % ha_cont_sub_stamp_band].data,
                                   norm=norm_ha_cont_sub_stamp, cmap='Greys')
            ax_hst_ha_cont_sub_stamp.set_title(ha_cont_sub_stamp_band.upper(), fontsize=fontsize_large)
            self.arr_axis_params(ax=ax_hst_ha_cont_sub_stamp, ra_tick_label=True, dec_tick_label=True,
                        ra_axis_label=' ', dec_axis_label=' ',
                        ra_minpad=0.8, dec_minpad=0.8, tick_color='k', label_color='k',
                        fontsize=fontsize_small, labelsize=fontsize_small, ra_tick_num=3, dec_tick_num=3)
            ax_hst_ha_cont_sub_cbar = figure.add_axes([stamp_left_align + (hst_stamp_index + 4)*(stamp_width + stamp_space_vertical) + c_bar_left_offset_hst_ha,
                                           stamp_bottom_align + (stamp_row_count_down+1)*(stamp_height + stamp_space_horizontal) + c_bar_top_offset,
                                           cbar_width_hst_ha, cbar_height],)
            helper_func.create_cbar(ax_cbar=ax_hst_ha_cont_sub_cbar, cmap='Greys', norm=norm_ha_cont_sub_stamp,
                                    cbar_label=cbar_label, fontsize=fontsize_large, orientation='horizontal',
                                    extend='neither')

        nircam_stamp_index = 0
        if include_nircam:
            stamp_row_count_down -= 1
            ax_nircam_stamp_list = []
            for nircam_stamp_index, nircam_stamp_band in enumerate(nircam_stamp_band_list):

                ax_nircam_stamp_list.append(
                    figure.add_axes([stamp_left_align + nircam_stamp_index*(stamp_width + stamp_space_vertical),
                                     stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                     stamp_width, stamp_height],
                                    projection=cutout_dict_stamp['%s_img_cutout' % nircam_stamp_band].wcs))
                ax_nircam_stamp_list[nircam_stamp_index].imshow(cutout_dict_stamp['%s_img_cutout' % nircam_stamp_band].data,
                                                norm=norm_nircam_stamps, cmap='Greys')
                ax_nircam_stamp_list[nircam_stamp_index].set_title(nircam_stamp_band.upper(), fontsize=fontsize_large)
                if nircam_stamp_index == 0:
                    ra_tick_label, dec_tick_label = True, True
                else:
                    ra_tick_label, dec_tick_label = False, False
                self.arr_axis_params(ax=ax_nircam_stamp_list[nircam_stamp_index], ra_tick_label=ra_tick_label,
                                     dec_tick_label=dec_tick_label, ra_axis_label=' ', dec_axis_label=' ',
                                     fontsize=fontsize_small, labelsize=fontsize_small)

            ax_nircam_cbar = figure.add_axes([stamp_left_align + c_bar_left_offset,
                                           stamp_bottom_align + (stamp_row_count_down+1)*(stamp_height + stamp_space_horizontal) + c_bar_top_offset,
                                           cbar_width, cbar_height],)
            helper_func.create_cbar(ax_cbar=ax_nircam_cbar, cmap='Greys', norm=norm_nircam_stamps, cbar_label=cbar_label,
                                    fontsize=fontsize_large, ticks=None, labelpad=2, tick_width=2, orientation='horizontal',
                                    extend='neither')

        if include_miri:
            # stamp_row_count_down -= 1
            ax_miri_stamp_list = []
            for miri_stamp_index, miri_stamp_band in enumerate(miri_stamp_band_list):

                ax_miri_stamp_list.append(
                    figure.add_axes([stamp_left_align + (nircam_stamp_index + miri_stamp_index + 2)*(stamp_width + stamp_space_vertical),
                                     stamp_bottom_align + stamp_row_count_down*(stamp_height + stamp_space_horizontal),
                                     stamp_width, stamp_height],
                                    projection=cutout_dict_stamp['%s_img_cutout' % miri_stamp_band].wcs))
                ax_miri_stamp_list[miri_stamp_index].imshow(cutout_dict_stamp['%s_img_cutout' % miri_stamp_band].data,
                                                norm=norm_miri_stamps, cmap='Greys')
                ax_miri_stamp_list[miri_stamp_index].set_title(miri_stamp_band.upper(), fontsize=fontsize_large)
                if miri_stamp_index == 0:
                    ra_tick_label, dec_tick_label = True, True
                else:
                    ra_tick_label, dec_tick_label = False, False
                self.arr_axis_params(ax=ax_miri_stamp_list[miri_stamp_index], ra_tick_label=ra_tick_label,
                                     dec_tick_label=dec_tick_label, ra_axis_label=' ', dec_axis_label=' ',
                                     fontsize=fontsize_small, labelsize=fontsize_small)

            ax_miri_cbar = figure.add_axes([stamp_left_align + (nircam_stamp_index + 2)*(stamp_width + stamp_space_vertical) + c_bar_left_offset,
                                           stamp_bottom_align + (stamp_row_count_down+1)*(stamp_height + stamp_space_horizontal) + c_bar_top_offset,
                                           cbar_width, cbar_height],)
            helper_func.create_cbar(ax_cbar=ax_miri_cbar, cmap='Greys', norm=norm_miri_stamps, cbar_label=cbar_label,
                                    fontsize=fontsize_large, ticks=None, labelpad=2, tick_width=2, orientation='horizontal',
                    extend='neither')

        # plot muse fit
        # 'wavelength': wavelength, 'total_flux': total_flux, 'total_flux_err': total_flux_err,
        # 'best_fit': best_fit, 'gas_best_fit': gas_best_fit, 'continuum_best_fit': continuum_best_fit,
        # 'pp': pp, 'ages': ages, 'met': met, 'mass2light': mass2light}

        # add axis
        ax_muse_spec = figure.add_axes([0.065, 0.025, 0.94, 0.2])
        ax_muse_stamp = figure.add_axes([0.8, 0.13, 0.18, 0.18], projection=muse_spec_dict['cutout_muse'].wcs)

        # plot spectrum
        ax_muse_spec.errorbar(muse_fit_dict['wavelength'], muse_fit_dict['total_flux'] * 1e-4,
                              yerr=muse_fit_dict['total_flux_err'] * 1e-4, fmt='.', color='k', label='Data')
        ax_muse_spec.plot(muse_fit_dict['wavelength'], muse_fit_dict['best_fit'] * 1e-4, color='tab:blue', label='Best fit')
        ax_muse_spec.plot(muse_fit_dict['wavelength'], muse_fit_dict['continuum_best_fit'] * 1e-4, color='tab:orange',
                          label='Stellar continuum')
        # ax_muse_spec.plot(muse_fit_dict['wavelength'], muse_fit_dict['gas_best_fit'])
        # get limits
        ax_muse_spec.set_xlim(np.nanmin(muse_fit_dict['wavelength']) - 10,
                              np.nanmax(muse_fit_dict['wavelength']) + 10)
        min_cont = np.nanmin(muse_fit_dict['continuum_best_fit'] * 1e-4)
        max_cont = np.nanmax(muse_fit_dict['continuum_best_fit'] * 1e-4)
        ax_muse_spec.set_ylim(min_cont - 0.05*(max_cont - min_cont), max_cont + 0.05*(max_cont - min_cont))
        ax_muse_spec.legend(frameon=False, fontsize=fontsize_large)

        ax_muse_spec.tick_params(axis='both', which='both', width=2, direction='in', labelsize=fontsize_large)
        ax_muse_spec.set_xlabel(r'Wavelength [$\AA$]', fontsize=fontsize_large)

        ax_muse_spec.set_ylabel(r'Flux density [10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$]', fontsize=fontsize_large)



        norm_muse_stamps = helper_func.compute_cbar_norm(cutout_list=[muse_spec_dict['cutout_muse'].data],
                                                         log_scale=cbar_log_scale)
        ax_muse_stamp.imshow(muse_spec_dict['cutout_muse'].data, norm=norm_muse_stamps)
        self.arr_axis_params(ax=ax_muse_stamp, ra_tick_label=True, dec_tick_label=True,
                        ra_axis_label=' ', dec_axis_label=' ',
                        ra_minpad=0.8, dec_minpad=0.8, tick_color='white', label_color='k',
                        fontsize=fontsize_large, labelsize=fontsize_large, ra_tick_num=3, dec_tick_num=3)

        figure.savefig('plot_output/example_ultimate_zoom.png')
        figure.clf()

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

    def get_target_hst_overview_rgb_img(self, red_band, green_band, blue_band, ref_band='red',
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
