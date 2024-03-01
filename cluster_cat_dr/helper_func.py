"""
This script gathers function to support the HST catalog release
"""

import os
from pathlib import Path
import warnings

from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from astropy import constants as const
speed_of_light_kmps = const.c.to('km/s').value

from reproject import reproject_interp
from TardisPipeline.readData.MUSE_WFM import get_MUSE_polyFWHM
import ppxf.ppxf_util as util

import numpy as np


def identify_file_in_folder(folder_path, str_in_file_name_1, str_in_file_name_2=None):
    """
    Identify a file inside a folder that contains a specific string.

    Parameters
    ----------
    folder_path : Path or str
    str_in_file_name_1 : str
    str_in_file_name_2 : str

    Returns
    -------
    file_name : Path
    """

    if str_in_file_name_2 is None:
        str_in_file_name_2 = str_in_file_name_1

    if isinstance(folder_path, str):
        folder_path = Path(folder_path)
    identified_files_1 = list(filter(lambda x: str_in_file_name_1 in x, os.listdir(folder_path)))

    identified_files_2 = list(filter(lambda x: str_in_file_name_2 in x, os.listdir(folder_path)))

    if not identified_files_1 and not identified_files_2:
        raise FileNotFoundError('The data file containing the string %s or %s does not exist.' %
                                (str_in_file_name_1, str_in_file_name_2))
    elif len(identified_files_1) > 1:
        raise FileExistsError('There are more than one data files containing the string %s .' % str_in_file_name_1)
    elif len(identified_files_2) > 1:
        raise FileExistsError('There are more than one data files containing the string %s .' % str_in_file_name_2)
    else:
        if not identified_files_2:
            return folder_path / str(identified_files_1[0])
        if not identified_files_1:
            return folder_path / str(identified_files_2[0])
        if identified_files_1 and identified_files_2:
            return folder_path / str(identified_files_1[0])


def load_img(file_name, hdu_number=0):
    """function to open hdu using astropy.

    Parameters
    ----------
    file_name : str or Path
        file name to open
    hdu_number : int or str
        hdu number which should be opened. can be also a string such as 'SCI' for JWST images

    Returns
    -------
    array-like, astropy.header and astropy.wcs object
    """
    # get hdu
    hdu = fits.open(file_name)
    # get header
    header = hdu[hdu_number].header
    # get WCS
    wcs = WCS(header)
    # update the header
    header.update(wcs.to_header())
    # reload the WCS and header
    header = hdu[hdu_number].header
    wcs = WCS(header)
    # load data
    data = hdu[hdu_number].data
    # close hdu again
    hdu.close()
    return data, header, wcs


def get_img_cutout(img, wcs, coord, cutout_size):
    """function to cut out a region of a larger image with an WCS.
    Parameters
    ----------
    img : ndarray
        (Ny, Nx) image
    wcs : astropy.wcs.WCS()
        astropy world coordinate system object describing the parameter image
    coord : astropy.coordinates.SkyCoord
        astropy coordinate object to point to the selected area which to cutout
    cutout_size : float or tuple
        Units in arcsec. Cutout size of a box cutout. If float it will be used for both box length.

    Returns
    -------
    cutout : astropy.nddata.Cutout2D object
        cutout object of the initial image
    """
    if isinstance(cutout_size, tuple):
        size = cutout_size * u.arcsec
    elif isinstance(cutout_size, float) | isinstance(cutout_size, int):
        size = (cutout_size, cutout_size) * u.arcsec
    else:
        raise KeyError('cutout_size must be float or tuple')

    # check if cutout is inside the image
    pix_pos = wcs.world_to_pixel(coord)
    if (pix_pos[0] > 0) & (pix_pos[0] < img.shape[1]) & (pix_pos[1] > 0) & (pix_pos[1] < img.shape[0]):
        return Cutout2D(data=img, position=coord, size=size, wcs=wcs)
    else:
        warnings.warn("The selected cutout is outside the original dataset. The data and WCS will be None",
                      DeprecationWarning)
        cut_out = type('', (), {})()
        cut_out.data = None
        cut_out.wcs = None
        return cut_out


def transform_world2pix_scale(length_in_arcsec, wcs, dim=0):
    """ Function to get the pixel length of a length in arcseconds
    Parameters
    ----------
    length_in_arcsec : float
        length
    wcs : ``astropy.wcs.WCS``
        astropy world coordinate system object describing the parameter image
    dim : int, 0 or 1
        specifys the dimension 0 for ra and 1 for dec. This should be however always the same values...

    Returns
    -------
    length_in_pixel : tuple
        length in pixel along ra and dec
    """

    return (length_in_arcsec*u.arcsec).to(u.deg) / wcs.proj_plane_pixel_scales()[dim]


def points_in_hull(p, hull, tol=1e-12):
    return np.all(hull.equations[:,:-1] @ p.T + np.repeat(hull.equations[:,-1][None,:], len(p), axis=0).T <= tol, 0)


def construct_wcs(ra_min, ra_max, dec_min, dec_max, img_shape, quadratic_image=True):
    """Function to generate a WCS from scratch by only using a box of coordinates and pixel sizes.
    Parameters
    ----------
    ra_min, ra_max, dec_min, dec_max,  : float
        outer coordinates of the new frame.
    img_shape : tuple
        number of pixels
    quadratic_image : bool
        flag whether the resulting WCS is quadratic or not

    Returns
    -------
    wcs : astropy.wcs.WCS()
        new WCS system centered on the coordinates
    """
    # get length of image
    pos_coord_lower_left = SkyCoord(ra=ra_min*u.deg, dec=dec_min*u.deg)
    pos_coord_lower_right = SkyCoord(ra=ra_max*u.deg, dec=dec_min*u.deg)
    pos_coord_upper_left = SkyCoord(ra=ra_min*u.deg, dec=dec_max*u.deg)
    # now get the size of the image
    ra_width = (pos_coord_lower_left.separation(pos_coord_lower_right)).degree
    dec_width = (pos_coord_lower_left.separation(pos_coord_upper_left)).degree

    # if we want to have a quadratic image we use the largest width
    if quadratic_image:
        ra_image_width = np.max([ra_width, dec_width])
        dec_image_width = np.max([ra_width, dec_width])
    else:
        ra_image_width = ra_width
        dec_image_width = dec_width

    # get central coordinates
    ra_center = (ra_min + ra_max) / 2
    dec_center = (dec_min + dec_max) / 2

    # now create a WCS for this histogram
    new_wcs = WCS(naxis=2)
    # what is the center pixel of the XY grid.
    new_wcs.wcs.crpix = [img_shape[0]/2, img_shape[1]/2]
    # what is the galactic coordinate of that pixel.
    new_wcs.wcs.crval = [ra_center, dec_center]
    # what is the pixel scale in lon, lat.
    new_wcs.wcs.cdelt = np.array([-ra_image_width / img_shape[0], dec_image_width / img_shape[1]])
    # you would have to determine if this is in fact a tangential projection.
    new_wcs.wcs.ctype = ["RA---AIR", "DEC--AIR"]

    return new_wcs


def reproject_image(data, wcs, new_wcs, new_shape):
    """function to reproject an image with na existing WCS to a new WCS
    Parameters
    ----------
    data : ndarray
    wcs : astropy.wcs.WCS()
    new_wcs : astropy.wcs.WCS()
    new_shape : tuple

    Returns
    -------
    new_data : ndarray
        new data reprojected to the new wcs
    """
    hdu = fits.PrimaryHDU(data=data, header=wcs.to_header())
    return reproject_interp(hdu, new_wcs, shape_out=new_shape, return_footprint=False)


def draw_box(ax, wcs, coord, box_size, color='k', line_width=2, line_style='-'):
    """
    function to draw a box around a coordinate on an axis with a WCS projection

    Parameters
    ----------
    ax : ``astropy.visualization.wcsaxes.core.WCSAxes``
    wcs : ``astropy.wcs.WCS``
    coord : ``astropy.coordinates.SkyCoord``
    box_size : tuple of float
        box size in arcsec
    color : str
    line_width : float
    line_style : str

    Returns
    -------
    None
    """
    if isinstance(box_size, tuple):
        box_size = box_size * u.arcsec
    elif isinstance(box_size, float) | isinstance(box_size, int):
        box_size = (box_size, box_size) * u.arcsec
    else:
        raise KeyError('cutout_size must be float or tuple')

    top_left_pix = wcs.world_to_pixel(SkyCoord(ra=coord.ra + (box_size[1] / 2)/np.cos(coord.dec.degree*np.pi/180),
                                               dec=coord.dec + (box_size[0] / 2)))
    top_right_pix = wcs.world_to_pixel(SkyCoord(ra=coord.ra - (box_size[1] / 2)/np.cos(coord.dec.degree*np.pi/180),
                                                dec=coord.dec + (box_size[0] / 2)))
    bottom_left_pix = wcs.world_to_pixel(SkyCoord(ra=coord.ra + (box_size[1] / 2)/np.cos(coord.dec.degree*np.pi/180),
                                                  dec=coord.dec - (box_size[0] / 2)))
    bottom_right_pix = wcs.world_to_pixel(SkyCoord(ra=coord.ra - (box_size[1] / 2)/np.cos(coord.dec.degree*np.pi/180),
                                                   dec=coord.dec - (box_size[0] / 2)))

    ax.plot([top_left_pix[0], top_right_pix[0]], [top_left_pix[1], top_right_pix[1]], color=color,
            linewidth=line_width, linestyle=line_style)
    ax.plot([bottom_left_pix[0], bottom_right_pix[0]], [bottom_left_pix[1], bottom_right_pix[1]], color=color,
            linewidth=line_width, linestyle=line_style)
    ax.plot([top_left_pix[0], bottom_left_pix[0]], [top_left_pix[1], bottom_left_pix[1]], color=color,
            linewidth=line_width, linestyle=line_style)
    ax.plot([top_right_pix[0], bottom_right_pix[0]], [top_right_pix[1], bottom_right_pix[1]], color=color,
            linewidth=line_width, linestyle=line_style)


def plot_coord_circle(ax, pos, rad, color, line_style='-', line_width=3, alpha=1., fill=False):
    """
    function to draw circles around a coordinate on an axis with a WCS projection

    Parameters
    ----------
    ax : ``astropy.visualization.wcsaxes.core.WCSAxes``
    pos : ``astropy.coordinates.SkyCoord``
    rad : float
        circle_radius in arcsec
    color : str
    line_style: str
    line_width: float
    alpha: float
    fill: bool



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
        for pos_i, rad_i, color_i, line_style_i, line_width_i, alpha_i in zip(pos, rad, color, line_style, line_width,
                                                                              alpha):
            circle = SphericalCircle(pos_i, rad_i * u.arcsec, edgecolor=color_i, facecolor=face_color,
                                     linewidth=line_width_i,
                                     linestyle=line_style_i, alpha=alpha_i, transform=ax.get_transform('fk5'))
            ax.add_patch(circle)
    else:
        circle = SphericalCircle(pos, rad * u.arcsec, edgecolor=color, facecolor=face_color, linewidth=line_width,
                                 linestyle=line_style, alpha=alpha, transform=ax.get_transform('fk5'))
        ax.add_patch(circle)


def load_muse_data(muse_cube_path):
    # get MUSE data
    muse_hdu = fits.open(muse_cube_path)
    # get header
    hdr = muse_hdu['DATA'].header
    # get wavelength
    wave_muse = hdr['CRVAL3'] + np.arange(hdr['NAXIS3']) * hdr['CD3_3']
    # get data and variance cube
    data_cube_muse = muse_hdu['DATA'].data
    var_cube_muse = muse_hdu['STAT'].data
    # get WCS
    wcs_muse = WCS(hdr).celestial

    muse_data_dict = {'wave_muse': wave_muse, 'data_cube_muse': data_cube_muse,
                      'var_cube_muse': var_cube_muse, 'wcs_muse': wcs_muse}

    return muse_data_dict


def extract_muse_spec_circ_app(muse_data_dict, ra, dec, circ_rad):
    # get select spectra from coordinates
    obj_coords_world = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    obj_coords_muse_pix = muse_data_dict['wcs_muse'].world_to_pixel(obj_coords_world)
    selection_radius_pix = transform_world2pix_scale(length_in_arcsec=circ_rad, wcs=muse_data_dict['wcs_muse'])

    x_lin_muse = np.linspace(1, muse_data_dict['data_cube_muse'].shape[2],
                             muse_data_dict['data_cube_muse'].shape[2])
    y_lin_muse = np.linspace(1, muse_data_dict['data_cube_muse'].shape[1],
                             muse_data_dict['data_cube_muse'].shape[1])
    x_data_muse, y_data_muse = np.meshgrid(x_lin_muse, y_lin_muse)
    mask_spectrum = (np.sqrt((x_data_muse - obj_coords_muse_pix[0]) ** 2 +
                             (y_data_muse - obj_coords_muse_pix[1]) ** 2) < selection_radius_pix)

    spec_flux = np.sum(muse_data_dict['data_cube_muse'][:, mask_spectrum], axis=1)
    spec_flux_err = np.sqrt(np.sum(muse_data_dict['var_cube_muse'][:, mask_spectrum], axis=1))

    lsf = get_MUSE_polyFWHM(muse_data_dict['wave_muse'], version="udf10")
    lam_range = [np.min(muse_data_dict['wave_muse'][np.invert(np.isnan(spec_flux))]),
                 np.max(muse_data_dict['wave_muse'][np.invert(np.isnan(spec_flux))])]
    lam = muse_data_dict['wave_muse']

    mask_wave_range = (lam > lam_range[0]) & (lam < lam_range[1])
    spec_flux = spec_flux[mask_wave_range]
    spec_flux_err = spec_flux_err[mask_wave_range]
    lam = lam[mask_wave_range]
    lsf = lsf[mask_wave_range]
    good_pixel_mask = np.invert(np.isnan(spec_flux) + np.isinf(spec_flux))

    return {'lam_range': lam_range, 'spec_flux': spec_flux, 'spec_flux_err': spec_flux_err, 'lam': lam,
            'lsf': lsf, 'good_pixel_mask': good_pixel_mask}


def fit_ppxf2spec(spec_dict, redshift, sps_name='fsps', age_range=None, metal_range=None):
    """

    Parameters
    ----------
    spec_dict : dict
    sps_name : str
        can be fsps, galaxev or emiles



    Returns
    -------
    dict
    """
    from os import path
    import ppxf.sps_util as lib
    from urllib import request
    from ppxf.ppxf import ppxf
    import matplotlib.pyplot as plt

    # plt.plot(spec_dict['lam'], spec_dict['spec_flux'])
    # plt.show()
    # rebin spectra
    # good_spec = np.invert(np.isnan(spec_dict['spec_flux']) + np.isnan(spec_dict['spec_flux_err']) + np.isinf(spec_dict['spec_flux']) + np.isinf(spec_dict['spec_flux_err']))
    # print(spec_dict['lam_range'])
    # spec_dict['spec_flux'] = spec_dict['spec_flux'][good_spec]
    # spec_dict['lam'] = spec_dict['lam'][good_spec]
    # spec_dict['spec_flux_err'] = spec_dict['spec_flux_err'][good_spec]
    # spec_dict['lam_range'] = [np.min(spec_dict['lam'][np.invert(np.isnan(spec_dict['spec_flux']))]),
    #                           np.max(spec_dict['lam'][np.invert(np.isnan(spec_dict['spec_flux']))])]
    # spec_dict['lsf'] = spec_dict['lsf'][good_spec]
    # print(spec_dict['lam_range'])
    # print(util.vac_to_air(spec_dict['lam']))
    # spec_dict['lam'] = util.vac_to_air(spec_dict['lam'])

    velscale = speed_of_light_kmps * np.diff(np.log(spec_dict['lam'][-2:]))[0]  # Smallest velocity step
    # velscale = speed_of_light_kmps*np.log(spec_dict['lam'][1]/spec_dict['lam'][0])

    spectra_muse, ln_lam_gal, velscale = util.log_rebin(lam=spec_dict['lam_range'], spec=spec_dict['spec_flux'],
                                                        velscale=velscale)
    spectra_muse_err, ln_lam_gal, velscale = util.log_rebin(lam=spec_dict['lam_range'],
                                                            spec=spec_dict['spec_flux_err'], velscale=velscale)

    # print(sum(np.isnan(spec_dict['spec_flux'])))
    # print(sum(np.isnan(spectra_muse)))
    #
    # plt.plot(ln_lam_gal, spectra_muse_err)
    # plt.show()

    lsf_dict = {"lam": spec_dict['lam'], "fwhm": spec_dict['lsf']}
    # get new wavelength array
    lam_gal = np.exp(ln_lam_gal)
    # goodpixels = util.determine_goodpixels(ln_lam=ln_lam_gal, lam_range_temp=spec_dict['lam_range'], z=redshift)
    goodpixels = None
    # goodpixels = (np.isnan(spectra_muse) + np.isnan(spectra_muse_err) + np.isinf(spectra_muse) + np.isinf(spectra_muse_err))
    # print(sum(np.invert(np.isnan(spectra_muse) + np.isnan(spectra_muse_err) + np.isinf(spectra_muse) + np.isinf(spectra_muse_err))))
    # print(sum(((spectra_muse > 0) & (spectra_muse < 100000000000000))))

    # get stellar library
    ppxf_dir = path.dirname(path.realpath(lib.__file__))
    basename = f"spectra_{sps_name}_9.0.npz"
    filename = path.join(ppxf_dir, 'sps_models', basename)
    if not path.isfile(filename):
        url = "https://raw.githubusercontent.com/micappe/ppxf_data/main/" + basename
        request.urlretrieve(url, filename)

    sps = lib.sps_lib(filename=filename, velscale=velscale, fwhm_gal=lsf_dict, norm_range=[5070, 5950],
                      wave_range=None,
                      age_range=age_range, metal_range=metal_range)
    reg_dim = sps.templates.shape[1:]  # shape of (n_ages, n_metal)
    stars_templates = sps.templates.reshape(sps.templates.shape[0], -1)

    gas_templates, gas_names, line_wave = util.emission_lines(ln_lam_temp=sps.ln_lam_temp,
                                                              lam_range_gal=spec_dict['lam_range'],
                                                              FWHM_gal=get_MUSE_polyFWHM)

    templates = np.column_stack([stars_templates, gas_templates])

    n_star_temps = stars_templates.shape[1]
    component = [0] * n_star_temps
    for line_name in gas_names:
        if '[' in line_name:
            component += [2]
        else:
            component += [1]

    gas_component = np.array(component) > 0  # gas_component=True for gas templates

    moments = [4, 4, 4]

    vel = speed_of_light_kmps * np.log(1 + redshift)   # eq.(8) of Cappellari (2017)
    start_gas = [vel, 150., 0, 0]     # starting guess
    start_star = [vel, 150., 0, 0]
    print(start_gas)
    start = [start_star, start_gas, start_gas]

    pp = ppxf(templates=templates, galaxy=spectra_muse, noise=spectra_muse_err, velscale=velscale, start=start,
              moments=moments, degree=-1, mdegree=4, lam=lam_gal, lam_temp=sps.lam_temp, #regul=1/rms,
              reg_dim=reg_dim, component=component, gas_component=gas_component, #reddening=0,
              gas_names=gas_names, goodpixels=goodpixels)

    light_weights = pp.weights[~gas_component]      # Exclude weights of the gas templates
    light_weights = light_weights.reshape(reg_dim)  # Reshape to (n_ages, n_metal)
    light_weights /= light_weights.sum()            # Normalize to light fractions

    # light_weights = pp.weights[~gas_component]      # Exclude weights of the gas templates
    # light_weights = light_weights.reshape(reg_dim)

    ages, met = sps.mean_age_metal(light_weights)
    mass2light = sps.mass_to_light(light_weights, redshift=redshift)

    return {'pp': pp, 'ages': ages, 'met': met, 'mass2light': mass2light}



    # wavelength = pp.lam
    # total_flux = pp.galaxy
    # total_flux_err = pp.noise
    #
    # best_fit = pp.bestfit
    # gas_best_fit = pp.gas_bestfit
    # continuum_best_fit = best_fit - gas_best_fit
    #
    # plt.errorbar(wavelength, total_flux, yerr=total_flux_err)
    # plt.plot(wavelength, continuum_best_fit + gas_best_fit)
    # plt.plot(wavelength, gas_best_fit)
    # plt.show()
    #
    #
    #
    #
    # plt.figure(figsize=(17, 6))
    # plt.subplot(111)
    # pp.plot()
    # plt.show()
    #
    # plt.figure(figsize=(9, 3))
    # sps.plot(light_weights)
    # plt.title("Light Weights Fractions");
    # plt.show()
    #
    # exit()

