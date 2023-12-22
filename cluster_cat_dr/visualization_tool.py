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
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground

from cluster_cat_dr import phot_data_access, helper_func
import multicolorfits as mcf


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
        axis_hight = 0.27
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
            if cutout_dict['%s_img_cutout' % nircam_band].data is not None:
                axis_dict.update({
                    'ax_%s' % nircam_band: figure.add_axes([start_x + (axis_width + axis_space_x) * band_index,
                                                            nircam_ax_y, axis_width, axis_hight],
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

        # now plot all bands
        for band in self.get_hst_nircam_miri_band_list():
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
                                               axis_hight], projection=cutout_dict['%s_img_cutout' % 'F336W'].wcs)
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
                                                          axis_hight],
                                                         projection=cutout_dict['%s_img_cutout' % nircam_blue_band].wcs)
                    })

                blue_data_nircam = cutout_dict['%s_img_cutout' % nircam_blue_band].data
                wcs_nircam = cutout_dict['%s_img_cutout' % nircam_blue_band].wcs
                green_data_nircam = self.reproject_image(data=cutout_dict['%s_img_cutout' % nircam_green_band].data,
                                                         wcs=cutout_dict['%s_img_cutout' % nircam_green_band].wcs,
                                                         new_wcs=wcs_nircam,
                                                         new_shape=blue_data_nircam.shape)
                red_data_nircam = self.reproject_image(data=cutout_dict['%s_img_cutout' % nircam_red_band].data,
                                                       wcs=cutout_dict['%s_img_cutout' % nircam_red_band].wcs,
                                                       new_wcs=wcs_nircam,
                                                       new_shape=blue_data_nircam.shape)

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
                                                        axis_hight],
                                                       projection=cutout_dict['%s_img_cutout' % miri_blue_band].wcs)
                    })

                blue_data_miri = cutout_dict['%s_img_cutout' % miri_blue_band].data
                wcs_miri = cutout_dict['%s_img_cutout' % miri_blue_band].wcs
                green_data_miri = self.reproject_image(data=cutout_dict['%s_img_cutout' % miri_green_band].data,
                                                       wcs=cutout_dict['%s_img_cutout' % miri_green_band].wcs,
                                                       new_wcs=wcs_miri,
                                                       new_shape=blue_data_miri.shape)
                red_data_miri = self.reproject_image(data=cutout_dict['%s_img_cutout' % miri_red_band].data,
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

        grey_r = mcf.greyRGBize_image(data_r, rescalefn='asinh', scaletype=scaletype_r, min_max=min_max_r, gamma=gamma_r)
        grey_g = mcf.greyRGBize_image(data_g, rescalefn='asinh', scaletype=scaletype_g, min_max=min_max_g, gamma=gamma_g)
        grey_b = mcf.greyRGBize_image(data_b, rescalefn='asinh', scaletype=scaletype_b, min_max=min_max_b, gamma=gamma_b)
        r = mcf.colorize_image(grey_r, color_r, colorintype='hex', gammacorr_color=gamma_corr_r)
        g = mcf.colorize_image(grey_g, color_g, colorintype='hex', gammacorr_color=gamma_corr_g)
        b = mcf.colorize_image(grey_b, color_b, colorintype='hex', gammacorr_color=gamma_corr_b)
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
