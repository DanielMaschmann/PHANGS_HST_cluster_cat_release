"""
Tool to visualize PHANGS imaging data
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.io import fits
from reproject import reproject_interp

import phot_data_access
import multicolorfits as mcf


class PhotVisualize(phot_data_access.PhotAccess):
    """
    Class to plot cutouts in multiple bands
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def plot_multi_band_artifacts(self, ra, dec, cluster_id, new_class, ml_class, cutout_size=2.0, circle_rad=0.2):
        """
        Function to create multi band panel for a single obejct

        Parameters
        ----------
        ra : float
            coordinate in deg without any units
        dec : float
            coordinate in deg without any units
        cluster_id : int
            id to identify the cluster
        new_class : str
            manually assigned object class
        ml_class : str
            ml assigned cluster class
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
        fontsize = 25
        axis_width = 0.27
        axis_hight = 0.27
        axis_space_x = -0.11
        start_x = -0.02
        hst_ax_y = 0.64
        nircam_ax_y = 0.345
        miri_ax_y = 0.05
        axis_dict = {}
        for band_index, hst_band in enumerate(self.get_hst_band_list()):
            axis_dict.update({
                'ax_%s' % hst_band: figure.add_axes([start_x + (axis_width + axis_space_x) * band_index, hst_ax_y,
                                                     axis_width, axis_hight],
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
            axis_dict.update({
                'ax_%s' % nircam_band: figure.add_axes([start_x + (axis_width + axis_space_x) * band_index, nircam_ax_y,
                                                        axis_width, axis_hight],
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
            axis_dict.update({
                'ax_%s' % miri_band: figure.add_axes([start_x + (axis_width + axis_space_x) * band_index, miri_ax_y,
                                                     axis_width, axis_hight],
                                                     projection=cutout_dict['%s_img_cutout' % miri_band].wcs)
            })

            if band_index == 0:
                self.arr_axis_params(ax=axis_dict['ax_%s' % miri_band], fontsize=fontsize, labelsize=fontsize,
                                     ra_tick_num=2)
            else:
                self.arr_axis_params(ax=axis_dict['ax_%s' % miri_band], dec_tick_label=False, dec_axis_label=' ',
                                     fontsize=fontsize, labelsize=fontsize, ra_tick_num=2)
            axis_dict['ax_%s' % miri_band].set_title(miri_band, fontsize=fontsize)

        # get rgb images
        axis_dict.update({
                'ax_hst_rgb': figure.add_axes([start_x+(axis_width+axis_space_x)*5, hst_ax_y, axis_width,
                                               axis_hight], projection=cutout_dict['%s_img_cutout' % 'F336W'].wcs)
            })

        hst_r_band = 'F555W'
        if 'F438W' in self.phangs_hst_obs_band_dict[self.target_name]['wfc3_uvis_observed_bands']:
            hst_g_band = 'F438W'
        else:
            hst_g_band = 'F435W'
        hst_b_band = 'F336W'

        hst_rgb = self.get_rgb_img(r_data=cutout_dict['%s_img_cutout' % hst_r_band].data,
                                   g_data=cutout_dict['%s_img_cutout' % hst_g_band].data,
                                   b_data=cutout_dict['%s_img_cutout' % hst_b_band].data)
        axis_dict['ax_hst_rgb'].imshow(hst_rgb)
        if nircam_obs_flag | miri_obs_flag:
            ra_axis_label = ' '
        else:
            ra_axis_label = 'R.A. (2000.0)'
        self.arr_axis_params(ax=axis_dict['ax_hst_rgb'], ra_tick_label=np.invert((nircam_obs_flag | miri_obs_flag)),
                             dec_tick_label=False, ra_axis_label=ra_axis_label, dec_axis_label=' ',
                             fontsize=fontsize, labelsize=fontsize)
        axis_dict['ax_hst_rgb'].set_title(hst_r_band + '+' + hst_g_band + '+' + hst_b_band, fontsize=fontsize)

        if nircam_obs_flag:
            axis_dict.update({
                    'ax_nircam_rgb': figure.add_axes([start_x+(axis_width+axis_space_x)*5, nircam_ax_y, axis_width,
                                                      axis_hight],
                                                     projection=cutout_dict['F200W_img_cutout'].wcs)
                })
            b_data = cutout_dict['F200W_img_cutout'].data
            g_data = self.reproject_image(data=cutout_dict['F300M_img_cutout'].data,
                                          wcs=cutout_dict['F300M_img_cutout'].wcs,
                                          new_wcs=cutout_dict['F200W_img_cutout'].wcs,
                                          new_shape=cutout_dict['F200W_img_cutout'].data.shape)
            r_data = self.reproject_image(data=cutout_dict['F360M_img_cutout'].data,
                                          wcs=cutout_dict['F360M_img_cutout'].wcs,
                                          new_wcs=cutout_dict['F200W_img_cutout'].wcs,
                                          new_shape=cutout_dict['F200W_img_cutout'].data.shape)
            nircam_rgb = self.get_rgb_img(r_data=b_data, g_data=g_data, b_data=r_data)
            axis_dict['ax_nircam_rgb'].imshow(nircam_rgb)
            if miri_obs_flag:
                ra_axis_label = ' '
            else:
                ra_axis_label = 'R.A. (2000.0)'
            self.arr_axis_params(ax=axis_dict['ax_nircam_rgb'], ra_tick_label=np.invert(miri_obs_flag),
                                 dec_tick_label=False, ra_axis_label=ra_axis_label, dec_axis_label=' ',
                                 fontsize=fontsize, labelsize=fontsize, ra_tick_num=2)
        axis_dict['ax_nircam_rgb'].set_title('F200W+F300M+F360M', fontsize=fontsize)

        if miri_obs_flag:
            axis_dict.update({
                    'ax_miri_rgb': figure.add_axes([start_x+(axis_width+axis_space_x)*5, miri_ax_y, axis_width,
                                                    axis_hight], projection=cutout_dict['F770W_img_cutout'].wcs)
                })

            b_data = cutout_dict['F770W_img_cutout'].data
            g_data = self.reproject_image(data=cutout_dict['F1000W_img_cutout'].data,
                                          wcs=cutout_dict['F1000W_img_cutout'].wcs,
                                          new_wcs=cutout_dict['F770W_img_cutout'].wcs,
                                          new_shape=cutout_dict['F770W_img_cutout'].data.shape)
            r_data = self.reproject_image(data=cutout_dict['F1130W_img_cutout'].data,
                                          wcs=cutout_dict['F1130W_img_cutout'].wcs,
                                          new_wcs=cutout_dict['F770W_img_cutout'].wcs,
                                          new_shape=cutout_dict['F770W_img_cutout'].data.shape)
            miri_rgb = self.get_rgb_img(r_data=b_data, g_data=g_data, b_data=r_data)
            axis_dict['ax_miri_rgb'].imshow(miri_rgb)
            self.arr_axis_params(ax=axis_dict['ax_miri_rgb'], dec_tick_label=False, dec_axis_label=' ',
                                 fontsize=fontsize, labelsize=fontsize, ra_tick_num=2)
        axis_dict['ax_miri_rgb'].set_title('F770W+F1000W+F1130W', fontsize=fontsize)

        # plot each band
        for band in self.get_hst_nircam_miri_band_list():
            m, s = (np.nanmean(cutout_dict['%s_img_cutout' % band].data),
                    np.nanstd(cutout_dict['%s_img_cutout' % band].data))
            axis_dict['ax_%s' % band].imshow(cutout_dict['%s_img_cutout' % band].data,
                                             cmap='Greys', vmin=m-s, vmax=m+5*s)
            pos = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
            self.plot_circle_on_wcs_img(ax=axis_dict['ax_%s' % band], pos=pos, rad=circle_rad, color='r')

        # add descriptive string
        string = 'ID_PHANGS_CLUSTERS_v1p2 = %i, ML_CLASS = %s, NEW_CLASS = %s' % (cluster_id, ml_class, new_class)
        plt.figtext(0.05, 0.97, string, fontsize=fontsize)

        return figure

    @staticmethod
    def reproject_image(data, wcs, new_wcs, new_shape):
        """
        Function to reproject an image to the resolution of another image while using the WCS as reference

        Parameters
        ----------
        data : ``numpy.ndarray``
        wcs : ``astropy.wcs.WCS``
        new_wcs : ``astropy.wcs.WCS``
        new_shape : tuple

        Returns
        -------
        reprojected_image : ``numpy.ndarray``
            reprojected image only the data not the WCS
        """
        hdu = fits.PrimaryHDU(data=data, header=wcs.to_header())
        return reproject_interp(hdu, new_wcs, shape_out=new_shape, return_footprint=False)

    @staticmethod
    def get_rgb_img(r_data, g_data, b_data, r_color='#FF4433', g_color='#0FFF50', b_color='#1F51FF',
                    r_min_max=None,
                    g_min_max=None,
                    b_min_max=None,
                    r_gamma=2.2, g_gamma=2.2, b_gamma=2.2,
                    r_gamma_corr=2.2, g_gamma_corr=2.2, b_gamma_corr=2.2, combined_gamma=2.2):
        """
        Function to create an RGB image

        Parameters
        ----------
        r_data, g_data, b_data : ``numpy.ndarray``, array
            color images. Must be all same shape
        r_color, g_color, b_color: str
            hex code for color
        r_min_max, g_min_max, b_min_max : tuple or None
            denotes the percentages till where the data is used
        r_gamma, g_gamma, b_gamma : float
            gamma factor for each individual color band
        r_gamma_corr, g_gamma_corr, b_gamma_corr : float
            gamma correction factor for each grey scale image
        combined_gamma : float
            gamma factor of resulting rgb image

        Returns
        -------
        rgb image : ``numpy.ndarray``
            of shape (N,N, 3)
        """
        if r_min_max is None:
            r_min_max = [0.01, 99.9]
        if g_min_max is None:
            g_min_max = [0.01, 99.9]
        if b_min_max is None:
            b_min_max = [0.01, 99.9]

        grey_r = mcf.greyRGBize_image(r_data, rescalefn='asinh', scaletype='perc', min_max=r_min_max, gamma=r_gamma)
        grey_g = mcf.greyRGBize_image(g_data, rescalefn='asinh', scaletype='perc', min_max=g_min_max, gamma=g_gamma)
        grey_b = mcf.greyRGBize_image(b_data, rescalefn='asinh', scaletype='perc', min_max=b_min_max, gamma=b_gamma)
        r = mcf.colorize_image(grey_r, r_color, colorintype='hex', gammacorr_color=r_gamma_corr)
        g = mcf.colorize_image(grey_g, g_color, colorintype='hex', gammacorr_color=g_gamma_corr)
        b = mcf.colorize_image(grey_b, b_color, colorintype='hex', gammacorr_color=b_gamma_corr)
        return mcf.combine_multicolor([r, g, b], gamma=combined_gamma)

    @staticmethod
    def plot_circle_on_wcs_img(ax, pos, rad, color, linestyle='-', linewidth=2., alpha=1., fill=False):
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
        linestyle : str
            matplotlib line style
        linewidth : float
        alpha : float
            matplotlib alpha factor
        fill : bool
            flag whether circle should be filled or not

        Returns
        -------
        None
        """

        if fill:
            facecolor = color
        else:
            facecolor = 'none'

        if isinstance(pos, list):
            if not isinstance(rad, list):
                rad = [rad] * len(pos)
            if not isinstance(color, list):
                color = [color] * len(pos)
            if not isinstance(linestyle, list):
                linestyle = [linestyle] * len(pos)
            if not isinstance(linewidth, list):
                linewidth = [linewidth] * len(pos)
            if not isinstance(alpha, list):
                alpha = [alpha] * len(pos)
            for pos_i, rad_i, color_i, linestyle_i, linewidth_i, alpha_i in zip(pos, rad, color, linestyle, linewidth,
                                                                                alpha):
                circle = SphericalCircle(pos_i, rad_i * u.arcsec, edgecolor=color_i, facecolor=facecolor,
                                         linewidth=linewidth_i,
                                         linestyle=linestyle_i, alpha=alpha_i, transform=ax.get_transform('fk5'))
                ax.add_patch(circle)
        else:
            circle = SphericalCircle(pos, rad * u.arcsec, edgecolor=color, facecolor=facecolor, linewidth=linewidth,
                                     linestyle=linestyle, alpha=alpha, transform=ax.get_transform('fk5'))
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
