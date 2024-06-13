"""
This script gathers function to support the HST catalog release
"""

import os
from pathlib import Path
import warnings

from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize, LogNorm

from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.stats import SigmaClip
from astropy.visualization.wcsaxes import SphericalCircle
from astropy import constants as const
speed_of_light_kmps = const.c.to('km/s').value

from reproject import reproject_interp
from TardisPipeline.readData.MUSE_WFM import get_MUSE_polyFWHM
import ppxf.ppxf_util as util
from scipy import odr
import sep

import numpy as np

from os import path
import ppxf.sps_util as lib
from urllib import request
from ppxf.ppxf import ppxf

from astropy.convolution import Gaussian2DKernel, convolve
from scipy.stats import gaussian_kde
from matplotlib import cm
import photometry_tools
import dust_tools.extinction_tools


def download_file(file_path, url, unpack=False, reload=False):
    """

    Parameters
    ----------
    file_path : str or ``pathlib.Path``
    url : str
    unpack : bool
        In case the downloaded file is zipped, this function can unpack it and remove the downloaded file,
        leaving only the extracted file
    reload : bool
        If the file is corrupted, this removes the file and reloads it

    Returns
    -------

    """
    if reload:
        # if reload file the file will be removed to re download it
        os.remove(file_path)
    # check if file already exists
    if os.path.isfile(file_path):
        print(file_path, 'already exists')
        return True
    else:
        from urllib3 import PoolManager
        # download file
        http = PoolManager()
        r = http.request('GET', url, preload_content=False)

        if unpack:
            with open(file_path.with_suffix(".gz"), 'wb') as out:
                while True:
                    data = r.read()
                    if not data:
                        break
                    out.write(data)
            r.release_conn()
            # uncompress file
            from gzip import GzipFile
            # read compressed file
            compressed_file = GzipFile(file_path.with_suffix(".gz"), 'rb')
            s = compressed_file.read()
            compressed_file.close()
            # save compressed file
            uncompressed_file = open(file_path, 'wb')
            uncompressed_file.write(s)
            uncompressed_file.close()
            # delete compressed file
            os.remove(file_path.with_suffix(".gz"))
        else:
            with open(file_path, 'wb') as out:
                while True:
                    data = r.read()
                    if not data:
                        break
                    out.write(data)
            r.release_conn()


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


def plot_coord_croshair(ax, pos, wcs, rad, hair_length, color, line_style='-', line_width=3, alpha=1.):
    """
    function to draw croshair around a coordinate on an axis with a WCS projection

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

    Returns
    -------
    None
    """

    pos_pix = wcs.world_to_pixel(pos)
    horizontal_rad = transform_world2pix_scale(length_in_arcsec=rad, wcs=wcs, dim=0)
    vertical_rad = transform_world2pix_scale(length_in_arcsec=rad, wcs=wcs, dim=1)
    horizontal_hair_length = transform_world2pix_scale(length_in_arcsec=hair_length, wcs=wcs, dim=0)
    vertical_hair_length = transform_world2pix_scale(length_in_arcsec=hair_length, wcs=wcs, dim=1)

    ax.plot([pos_pix[0] + horizontal_rad, pos_pix[0] + horizontal_rad + horizontal_hair_length],
            [pos_pix[1], pos_pix[1]], color=color, linestyle=line_style, linewidth=line_width, alpha=alpha)
    ax.plot([pos_pix[0] - horizontal_rad, pos_pix[0] - horizontal_rad - horizontal_hair_length],
            [pos_pix[1], pos_pix[1]], color=color, linestyle=line_style, linewidth=line_width, alpha=alpha)
    ax.plot([pos_pix[0], pos_pix[0]], [pos_pix[1] + vertical_rad, pos_pix[1] + vertical_rad + vertical_hair_length],
        color=color, linestyle=line_style, linewidth=line_width, alpha=alpha)
    ax.plot([pos_pix[0], pos_pix[0]], [pos_pix[1] - vertical_rad, pos_pix[1] - vertical_rad - vertical_hair_length],
        color=color, linestyle=line_style, linewidth=line_width, alpha=alpha)


def extract_flux_from_circ_aperture(data, wcs, pos, aperture_rad, data_err=None):
    """

    Parameters
    ----------
    data : ``numpy.ndarray``
    wcs : ``astropy.wcs.WCS``
    pos : ``astropy.coordinates.SkyCoord``
    aperture_rad : float
    data_err : ``numpy.ndarray``

    Returns
    -------
    flux : float
    flux_err : float
    """
    # estimate background
    bkg = sep.Background(np.array(data, dtype=float))
    # get radius in pixel scale
    pix_radius = transform_world2pix_scale(length_in_arcsec=aperture_rad, wcs=wcs, dim=1)
    # pix_radius_old = (wcs.world_to_pixel(pos)[0] -
    #               wcs.world_to_pixel(SkyCoord(ra=pos.ra + aperture_rad * u.arcsec, dec=pos.dec))[0])
    # print(pix_radius)
    # print(pix_radius_old)
    # exit()
    # get the coordinates in pixel scale
    pixel_coords = wcs.world_to_pixel(pos)

    data = np.array(data.byteswap().newbyteorder(), dtype=float)
    if data_err is None:
        bkg_rms = bkg.rms()
        data_err = np.array(bkg_rms.byteswap().newbyteorder(), dtype=float)
    else:
        data_err = np.array(data_err.byteswap().newbyteorder(), dtype=float)

    flux, flux_err, flag = sep.sum_circle(data=data - bkg.globalback, x=np.array([float(pixel_coords[0])]),
                                          y=np.array([float(pixel_coords[1])]), r=np.array([float(pix_radius)]),
                                          err=data_err)

    return float(flux), float(flux_err)

def load_muse_cube(muse_cube_path):
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

    muse_hdu.close()

    muse_data_dict = {'wave_muse': wave_muse, 'data_cube_muse': data_cube_muse,
                      'var_cube_muse': var_cube_muse, 'wcs_muse': wcs_muse, 'hdr_muse': hdr}

    return muse_data_dict


def load_muse_dap_map(muse_dap_map_path, map='HA6562_FLUX'):
    # get MUSE data
    muse_hdu = fits.open(muse_dap_map_path)
    muse_map_data = muse_hdu[map].data
    muse_map_wcs = WCS(muse_hdu[map].header)
    muse_hdu.close()
    return {'muse_map_data': muse_map_data, 'muse_map_wcs': muse_map_wcs}


def extract_muse_spec_circ_app(muse_data_dict, ra, dec, circ_rad, cutout_size, wave_range=None):
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
    if wave_range is None:
        lam_range = [np.min(muse_data_dict['wave_muse'][np.invert(np.isnan(spec_flux))]),
                     np.max(muse_data_dict['wave_muse'][np.invert(np.isnan(spec_flux))])]
    else:
        lam_range = wave_range
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

    import matplotlib.pyplot as plt
    spec_dict['spec_flux'] *= 1e-20
    spec_dict['spec_flux_err'] *= 1e-20

    velscale = speed_of_light_kmps * np.diff(np.log(spec_dict['lam'][-2:]))[0]  # Smallest velocity step
    # print('velscale ', velscale)
    # velscale = speed_of_light_kmps*np.log(spec_dict['lam'][1]/spec_dict['lam'][0])
    # print('velscale ', velscale)
    # velscale = 10
    spectra_muse, ln_lam_gal, velscale = util.log_rebin(lam=spec_dict['lam'], spec=spec_dict['spec_flux'],
                                                        velscale=velscale)
    spectra_muse_err, ln_lam_gal, velscale = util.log_rebin(lam=spec_dict['lam'],
                                                            spec=spec_dict['spec_flux_err'], velscale=velscale)

    lsf_dict = {"lam": spec_dict['lam'], "fwhm": spec_dict['lsf']}
    # get new wavelength array
    lam_gal = np.exp(ln_lam_gal)

    # get stellar library
    ppxf_dir = path.dirname(path.realpath(lib.__file__))
    basename = f"spectra_{sps_name}_9.0.npz"
    filename = path.join(ppxf_dir, 'sps_models', basename)
    if not path.isfile(filename):
        url = "https://raw.githubusercontent.com/micappe/ppxf_data/main/" + basename
        request.urlretrieve(url, filename)

    sps = lib.sps_lib(filename=filename, velscale=velscale, fwhm_gal=lsf_dict, norm_range=[5070, 5950],
                      wave_range=None, age_range=age_range, metal_range=metal_range)
    reg_dim = sps.templates.shape[1:]  # shape of (n_ages, n_metal)
    stars_templates = sps.templates.reshape(sps.templates.shape[0], -1)

    gas_templates, gas_names, line_wave = util.emission_lines(ln_lam_temp=sps.ln_lam_temp,
                                                              lam_range_gal=spec_dict['lam_range'],
                                                              FWHM_gal=get_MUSE_polyFWHM,
                                                              limit_doublets=False)

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
    start = [start_star, start_gas, start_gas]

    # mask bad values
    mask = np.invert(np.isnan(spectra_muse_err))
    spectra_muse_err[np.isnan(spectra_muse_err)] = np.nanmean(spectra_muse_err)
    spectra_muse[np.isnan(spectra_muse)] = 0


    pp = ppxf(templates=templates, galaxy=spectra_muse, noise=spectra_muse_err, velscale=velscale, start=start,
              moments=moments, degree=-1, mdegree=4, lam=lam_gal, lam_temp=sps.lam_temp,
              reg_dim=reg_dim, component=component, gas_component=gas_component,
              reddening=2.5, gas_reddening=0.0, gas_names=gas_names, mask=mask)

    light_weights = pp.weights[~gas_component]      # Exclude weights of the gas templates
    light_weights = light_weights.reshape(reg_dim)  # Reshape to (n_ages, n_metal)
    light_weights /= light_weights.sum()            # Normalize to light fractions

    ages, met = sps.mean_age_metal(light_weights)
    mass2light = sps.mass_to_light(light_weights, redshift=redshift)

    wavelength = pp.lam
    total_flux = pp.galaxy
    total_flux_err = pp.noise

    best_fit = pp.bestfit
    gas_best_fit = pp.gas_bestfit
    continuum_best_fit = best_fit - gas_best_fit

    # get velocity of balmer component
    sol_kin_comp = pp.sol[0]
    balmer_kin_comp = pp.sol[1]
    forbidden_kin_comp = pp.sol[2]

    h_beta_rest_air = 4861.333
    h_alpha_rest_air = 6562.819

    balmer_redshift = np.exp(balmer_kin_comp[0] / speed_of_light_kmps) - 1

    observed_h_beta = h_beta_rest_air * (1 + balmer_redshift)
    observed_h_alpha = h_alpha_rest_air * (1 + balmer_redshift)

    observed_sigma_h_alpha = (balmer_kin_comp[1] / speed_of_light_kmps) * h_alpha_rest_air
    observed_sigma_h_alpha = np.sqrt(observed_sigma_h_alpha ** 2 + get_MUSE_polyFWHM(observed_h_alpha))
    observed_sigma_h_beta = (balmer_kin_comp[1] / speed_of_light_kmps) * h_beta_rest_air
    observed_sigma_h_beta = np.sqrt(observed_sigma_h_beta ** 2 + get_MUSE_polyFWHM(observed_h_beta))

    mask_ha = (wavelength > (observed_h_alpha - 3*observed_sigma_h_alpha)) & (wavelength < (observed_h_alpha + 3*observed_sigma_h_alpha))
    mask_hb = (wavelength > (observed_h_beta - 3*observed_sigma_h_beta)) & (wavelength < (observed_h_beta + 3*observed_sigma_h_beta))

    ha_line_comp = (total_flux - continuum_best_fit)[mask_ha]
    ha_cont_comp = continuum_best_fit[mask_ha]
    ha_wave_comp = wavelength[mask_ha]
    delta_lambda_ha = np.mean((ha_wave_comp[1:] - ha_wave_comp[:-1]) / 2)
    ha_ew = np.sum(((ha_cont_comp - ha_line_comp) / ha_cont_comp) * delta_lambda_ha)

    hb_line_comp = (total_flux - continuum_best_fit)[mask_hb]
    hb_cont_comp = continuum_best_fit[mask_hb]
    hb_wave_comp = wavelength[mask_hb]
    delta_lambda_hb = np.mean((hb_wave_comp[1:] - hb_wave_comp[:-1]) / 2)
    hb_ew = np.sum(((hb_cont_comp - hb_line_comp) / hb_cont_comp) * delta_lambda_hb)

    # gas_phase_metallicity
    flux_ha = pp.gas_flux[pp.gas_names == 'Halpha']
    flux_hb = pp.gas_flux[pp.gas_names == 'Hbeta']
    flux_nii = pp.gas_flux[pp.gas_names == '[OIII]5007_d']
    flux_oiii = pp.gas_flux[pp.gas_names == '[NII]6583_d']

    o3n2 = np.log10((flux_oiii / flux_hb) / (flux_nii / flux_ha))
    gas_phase_met = 8.73 - 0.32 * o3n2[0]
    # plt.plot(hb_wave_comp, hb_line_comp)
    # plt.plot(hb_wave_comp, hb_cont_comp)
    # plt.show()
    # exit()
    #
    # # exit()
    # plt.errorbar(wavelength, total_flux, yerr=total_flux_err)
    # plt.plot(wavelength, continuum_best_fit)
    # plt.scatter(wavelength[left_idx_ha[0][0]], continuum_best_fit[left_idx_ha[0][0]])
    # plt.scatter(wavelength[right_idx_ha[0][0]], continuum_best_fit[right_idx_ha[0][0]])
    # plt.plot(wavelength, continuum_best_fit + gas_best_fit)
    # plt.plot(wavelength, gas_best_fit)
    # plt.plot([observed_nii_1, observed_nii_1], [np.min(total_flux), np.max(total_flux)])
    # plt.plot([observed_h_alpha, observed_h_alpha], [np.min(total_flux), np.max(total_flux)])
    # plt.plot([observed_nii_2, observed_nii_2], [np.min(total_flux), np.max(total_flux)])
    # plt.show()
    #

    return {
        'wavelength': wavelength, 'total_flux': total_flux, 'total_flux_err': total_flux_err,
        'best_fit': best_fit, 'gas_best_fit': gas_best_fit, 'continuum_best_fit': continuum_best_fit,
        'ages': ages, 'met': met, 'mass2light': mass2light,
        'star_red': pp.dust[0]['sol'][0], 'gas_red': pp.dust[1]['sol'][0],
        'sol_kin_comp': sol_kin_comp, 'balmer_kin_comp': balmer_kin_comp, 'forbidden_kin_comp': forbidden_kin_comp,
        'ha_ew': ha_ew, 'hb_ew': hb_ew, 'gas_phase_met': gas_phase_met
    }

    #
    #
    #
    # plt.figure(figsize=(17, 6))
    # plt.subplot(111)
    # pp.plot()
    # plt.show()
    #

    #
    # exit()


def fit_tardis2spec(spec_dict, velocity, hdr, sps_name='fsps', age_range=None, metal_range=None, name='explore1'):
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
    # import ppxf.sps_util as lib
    # from urllib import request
    # from ppxf.ppxf import ppxf

    import matplotlib.pyplot as plt

    from TardisPipeline.utilities import util_ppxf, util_ppxf_stellarpops, util_sfh_quantities,util_ppxf_emlines
    import TardisPipeline as tardis_module
    codedir = os.path.dirname(os.path.realpath(tardis_module.__file__))

    import ppxf.ppxf_util as util
    from astropy.io import fits, ascii
    from astropy import constants as const
    from astropy.table import Table
    import extinction

    # tardis_path = '/home/egorov/Soft/ifu-pipeline/TardisPipeline/' # change to directory where you have installed DAP
    ncpu = 20  # how many cpu would you like to use? (20-30 is fine for our server, but use no more than 8 for laptop)
    # print(codedir+'/Templates/spectralTemplates/eMILES-noyoung/')
    # exit()
    configs = {#'SSP_LIB': os.path.join(codedir, 'Templates/spectralTemplates/eMILES-noyoung/'),
               #'SSP_LIB_SFH': os.path.join(codedir, 'Templates/spectralTemplates/eMILES-noyoung/'),
               'SSP_LIB': codedir+'/Templates/spectralTemplates/CB07_chabrier-young-selection-MetalPoorRemoved/',  # stellar library to use
               'SSP_LIB_SFH': codedir+'/Templates/spectralTemplates/CB07_chabrier-young-selection-MetalPoorRemoved/',  # stellar library to use
        # 'SSP_LIB': codedir+'/Templates/spectralTemplates/eMILES-noyoung/',  # stellar library to use
               'NORM_TEMP': 'LIGHT', 'REDSHIFT': velocity, 'MOM': 4, 'MC_PPXF': 0, 'PARALLEL': 1,
        'ADEG': 12,
        'ADEG_SFH': 12,
               'MDEG': 0,
               'MDEG_SFH': 0,
        'MDEG_EMS': 24,
        'NCPU': ncpu,
        'ROOTNAME': name,
        'SPECTRUM_SIZE': abs(hdr['CD1_1']) * 3600.,  # spaxel size in arcsec
               # 'EMI_FILE': os.path.join(codedir, '/Templates/configurationTemplates/emission_lines.setup'),
        'MC_PPXF_SFH': 10,
                'EMI_FILE': codedir+'/Templates/configurationTemplates/emission_lines.setup',  # set of emission lines to fit
        'SKY_LINES_RANGES': codedir+'/Templates/configurationTemplates/sky_lines_ranges.setup',
               'OUTDIR': 'data_output/',
        'MASK_WIDTH': 150,
    'GAS_MOMENTS': 4}


    velscale = speed_of_light_kmps * np.diff(np.log(spec_dict['lam'][-2:]))[0]  # Smallest velocity step
    log_spec, logLam, velscale = util.log_rebin(lam=spec_dict['lam_range'], spec=spec_dict['spec_flux'],
                                                        velscale=velscale)
    c1=fits.Column(name='LOGLAM',array=logLam,format='D')
    c2 = fits.Column(name='LOGSPEC', array=log_spec,format='D')
    t=fits.BinTableHDU.from_columns([c1,c2])
    t.writeto('{}{}-ppxf_obsspec.fits'.format(configs['OUTDIR'],name),overwrite=True)
    log_err, _, _ = util.log_rebin(spec_dict['lam_range'], spec_dict['spec_flux_err'], velscale=velscale)
    ww = ~np.isfinite(log_spec) | ~np.isfinite(log_err) | (log_err <= 0)
    log_err[ww] = 9999
    log_spec[ww] = 0.
    # # the DAP fitting routines expect log_spec and log_err to be 2D arrays containing N spectra,
    # # here we add a dummy dimension since we are fitting only one spectrum
    # # to fit more than one spectrum at the same time these lines can be easily adapted
    log_err = np.expand_dims(log_err, axis=1)
    log_spec = np.expand_dims(log_spec, axis=1)

    # define the LSF of the MUSE data
    LSF = get_MUSE_polyFWHM(np.exp(logLam), version="udf10")

    # define the velocity scale in kms
    velscale = (logLam[1] - logLam[0]) * speed_of_light_kmps

# this is the stellar kinematics ppxf wrapper function
    ppxf_result = util_ppxf.runModule_PPXF(configs=configs, #tasks='',
                                           logLam=logLam,
                                           log_spec=log_spec, log_error=log_err,
                                           LSF=LSF)#, velscale=velscale)
    util_ppxf_emlines.runModule_PPXF_emlines(configs=configs, #tasks='',
                                             logLam=logLam,
                                             log_spec=log_spec, log_error=log_err,
                                             LSF=LSF, ppxf_results=ppxf_result)

    # exit()
    util_ppxf_stellarpops.runModule_PPXF_stellarpops(configs, logLam, log_spec, log_err, LSF, np.arange(1), ppxf_result)
    masses_density, mass_density_err, ages_mw, ages_mw_err, z_mw, z_mw_err, ages_lw, ages_lw_err, z_lw, z_lw_err = util_sfh_quantities.compute_sfh_relevant_quantities(configs)
    print(masses_density, mass_density_err, ages_mw, ages_mw_err, z_mw, z_mw_err, ages_lw, ages_lw_err, z_lw, z_lw_err )


    # read the output file which contains the best-fit from the emission lines fitting stage
    ppxf_bestfit_gas = configs['OUTDIR']+configs['ROOTNAME']+'_ppxf-bestfit-emlines.fits'
    hdu3 = fits.open(ppxf_bestfit_gas)
    bestfit_gas = hdu3['FIT'].data["BESTFIT"][0]
    mask = (hdu3['FIT'].data['BESTFIT'][0]==0)
    gas_templ = hdu3['FIT'].data["GAS_BESTFIT"][0]

    ppxf_bestfit = configs['OUTDIR']+configs['ROOTNAME']+'_ppxf-bestfit.fits'
    hdu_best_fit = fits.open(ppxf_bestfit)
    cont_fit = hdu_best_fit['FIT'].data["BESTFIT"][0]



    # # reddening = ppxf_sfh_data['REDDENING']
    # hdu_best_fit_sfh = fits.open('data_output/explore1_ppxf-bestfit.fits')
    # print(hdu_best_fit_sfh.info())
    # print(hdu_best_fit_sfh[1].data.names)
    #
    # print(hdu_best_fit_sfh['FIT'].data['BESTFIT'])
    # print(hdu_best_fit_sfh['FIT'].data['BESTFIT'].shape)
    # print(logLam.shape)
    # print(spec_dict['lam'].shape)
    # # exit()
    # # hdu_best_fit = fits.open('data_output/explore1_templates_SFH_info.fits')
    # # print(hdu_best_fit.info())
    # # print(hdu_best_fit[1].data.names)
    # # print(hdu_best_fit[1].data['Age'])

    plt.plot(spec_dict['lam'], spec_dict['spec_flux'])
    plt.plot(np.exp(logLam), cont_fit)
    plt.plot(np.exp(logLam), gas_templ)
    plt.plot(np.exp(logLam), cont_fit+gas_templ)
    plt.show()


    exit()
    # this the ppxf wrapper function to simulataneously fit the continuum plus emission lines
    # util_ppxf_emlines.runModule_PPXF_emlines(configs,# '',
    #                                          logLam, log_spec,
    #                                          log_err, LSF, #velscale,
    #                                          np.arange(1), ppxf_result)
    util_ppxf_emlines.runModule_PPXF_emlines(configs=configs, #tasks='',
                                             logLam=logLam,
                                             log_spec=log_spec, log_error=log_err,
                                             LSF=LSF, ppxf_results=ppxf_result)

    emlines = configs['OUTDIR'] + configs['ROOTNAME'] + '_emlines.fits'
    with fits.open(emlines) as hdu_emis:
        ems = Table(hdu_emis['EMLDATA_DATA'].data)

    # This is to include SFH results, NOT TESTED!
    with fits.open(configs['OUTDIR'] + configs['ROOTNAME'] + '_ppxf_SFH.fits') as hdu_ppxf_sfh:
        ppxf_sfh_data = hdu_ppxf_sfh[1].data
        masses_density, mass_density_err, ages_mw, ages_mw_err, z_mw, z_mw_err, ages_lw, ages_lw_err, z_lw, z_lw_err = util_sfh_quantities.compute_sfh_relevant_quantities(configs)
        reddening = ppxf_sfh_data['REDDENING']
        st_props = masses_density, mass_density_err, ages_mw, ages_mw_err, z_mw, z_mw_err, ages_lw, ages_lw_err, z_lw, z_lw_err, reddening

    exit()

    return ems, st_props




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


def compute_cbar_norm(vmin_vmax=None, cutout_list=None, log_scale=False):
    """
    Computing the color bar scale for a single or multiple cutouts.

    Parameters
    ----------
    vmin_vmax : tuple
    cutout_list : list
        This list should include all cutouts
    log_scale : bool

    Returns
    -------
    norm : ``matplotlib.colors.Normalize``  or ``matplotlib.colors.LogNorm``
    """
    if (vmin_vmax is None) & (cutout_list is None):
        raise KeyError('either vmin_vmax or cutout_list must be not None')

    # get maximal value
    # vmin_vmax

    if vmin_vmax is None:
        vmin = None
        vmax = None
        for cutout in cutout_list:
            sigma_clip = SigmaClip(sigma=3)
            mask_zeros = np.invert(cutout == 0)
            if len(sigma_clip(cutout[mask_zeros])) == 0:
                return None
            min = np.nanmin(sigma_clip(cutout[mask_zeros]))
            max = np.nanmax(sigma_clip(cutout[mask_zeros]))
            if vmin is None:
               vmin = min
            if vmax is None:
               vmax = max
            if min < vmin:
                vmin = min
            if max > vmax:
                vmax = max

        # list_of_means = [np.nanmean(cutout) for cutout in cutout_list]
        # list_of_stds = [np.nanstd(cutout) for cutout in cutout_list]
        # mean, std = (np.nanmean(list_of_means), np.nanstd(list_of_stds))
        #
        # vmin = mean - 5 * std
        # vmax = mean + 20 * std


    else:
        vmin, vmax = vmin_vmax[0], vmin_vmax[1]
    if log_scale:

        if vmax < 0:
            vmax = 0.000001
        if vmin < 0:
            vmin = vmax / 100
        norm = LogNorm(vmin, vmax)
    else:
        norm = Normalize(vmin, vmax)
    return norm


def create_cbar(ax_cbar, cmap, norm, cbar_label, fontsize, ticks=None, labelpad=2, tick_width=2, orientation='vertical',
                extend='neither'):
    """

    Parameters
    ----------
    ax_cbar : ``matplotlib.pylab.axis``
    cmap : str
        same as name parameter of ``matplotlib.colors.Colormap.name``
    norm : ``matplotlib.colors.Normalize``  or ``matplotlib.colors.LogNorm``
    cbar_label : str
    fontsize : int or float
    ticks : list
    labelpad : int or float
    tick_width : int or float
    orientation : str
        default is `vertical`
    extend : str
        default is 'neither'
        can be 'neither', 'min' , 'max' or 'both'
    """
    ColorbarBase(ax_cbar, orientation=orientation, cmap=cmap, norm=norm, extend=extend, ticks=ticks)
    if orientation == 'vertical':
        ax_cbar.set_ylabel(cbar_label, labelpad=labelpad, fontsize=fontsize)
        ax_cbar.tick_params(axis='both', which='both', width=tick_width, direction='in', top=True, labelbottom=False,
                            labeltop=True, labelsize=fontsize)
    elif orientation == 'horizontal':
        # ax_cbar.set_xlabel(cbar_label, labelpad=labelpad, fontsize=fontsize)
        ax_cbar.tick_params(width=tick_width, direction='in', top=True, labeltop=True, bottom=False, labelbottom=False,
                            labelsize=fontsize)
        # also put the minor ticks to the top
        ax_cbar.tick_params(which='minor', width=tick_width, direction='in',
                            top=True, labeltop=True, bottom=False, labelbottom=False,
                            labelsize=fontsize/1.5)
        ax_cbar.set_title(cbar_label, fontsize=fontsize)


def lin_func(p, x):
    gradient, intersect = p
    return gradient*x + intersect


def fit_line(x_data, y_data, x_data_err, y_data_err):

    # Create a model for fitting.
    lin_model = odr.Model(lin_func)

    # Create a RealData object using our initiated data from above.
    data = odr.RealData(x_data, y_data, sx=x_data_err, sy=y_data_err)

    # Set up ODR with the model and data.
    odr_object = odr.ODR(data, lin_model, beta0=[0., 1.])

    # Run the regression.
    out = odr_object.run()

    # Use the in-built pprint method to give us results.
    # out.pprint()

    gradient, intersect = out.beta
    gradient_err, intersect_err = out.sd_beta

    # calculate sigma around fit
    sigma = np.std(y_data - lin_func(p=(gradient, intersect), x=x_data))

    return {
        'gradient': gradient,
        'intersect': intersect,
        'gradient_err': gradient_err,
        'intersect_err': intersect_err,
        'sigma': sigma
    }


def conv_mjy2ab_mag(flux):

    """conversion of mJy to AB mag.
    See definition on Wikipedia : https://en.wikipedia.org/wiki/AB_magnitude """

    return -2.5 * np.log10(flux * 1e-3) + 8.90


def conv_ab_mag2mjy(mag):

    """conversion of AB mag to mJy.
    See definition on Wikipedia : https://en.wikipedia.org/wiki/AB_magnitude """

    return 1e3 * 10 ** ((8.5 - mag)/2.5)


def conv_mjy2vega(flux, ab_zp, vega_zp):

    """This function (non-sophisticated as of now)
     assumes the flux are given in units of milli-Janskies"""

    """First, convert mjy to f_nu"""
    conv_f_nu = flux*np.power(10.0, -26)
    """Convert f_nu in ergs s^-1 cm^-2 Hz^-1 to AB mag"""
    ABmag = -2.5*np.log10(conv_f_nu) - 48.60
    """Convert AB mag to Vega mag"""
    vega_mag = ABmag + (vega_zp - ab_zp)

    return vega_mag


def density_with_points(ax, x, y, binx=None, biny=None, threshold=1, kernel_std=2.0, save=False, save_name='',
                        cmap='inferno', scatter_size=10, scatter_alpha=0.3, invert_y_axis=False):

    if binx is None:
        binx = np.linspace(-1.5, 2.5, 190)
    if biny is None:
        biny = np.linspace(-2.0, 2.0, 190)

    good = np.invert(((np.isnan(x) | np.isnan(y)) | (np.isinf(x) | np.isinf(y))))

    hist, xedges, yedges = np.histogram2d(x[good], y[good], bins=(binx, biny))

    if save:
        np.save('data_output/binx.npy', binx)
        np.save('data_output/biny.npy', biny)
        np.save('data_output/hist_%s_un_smoothed.npy' % save_name, hist)

    kernel = Gaussian2DKernel(x_stddev=kernel_std)
    hist = convolve(hist, kernel)

    if save:
        np.save('data_output/hist_%s_smoothed.npy' % save_name, hist)

    over_dense_regions = hist > threshold
    mask_high_dens = np.zeros(len(x), dtype=bool)

    for x_index in range(len(xedges)-1):
        for y_index in range(len(yedges)-1):
            if over_dense_regions[x_index, y_index]:
                mask = (x > xedges[x_index]) & (x < xedges[x_index + 1]) & (y > yedges[y_index]) & (y < yedges[y_index + 1])
                mask_high_dens += mask
    print(sum(mask_high_dens) / len(mask_high_dens))
    hist[hist <= threshold] = np.nan

    cmap = cm.get_cmap(cmap)

    scatter_color = cmap(0)

    ax.imshow(hist.T, origin='lower', extent=(binx.min(), binx.max(), biny.min(), biny.max()), cmap=cmap,
              interpolation='nearest', aspect='auto')
    ax.scatter(x[~mask_high_dens], y[~mask_high_dens], color=scatter_color, marker='.', s=scatter_size, alpha=scatter_alpha)
    if invert_y_axis:
        ax.set_ylim(ax.get_ylim()[::-1])





def plot_reddening_vect(ax, x_color_1='v', x_color_2='i',  y_color_1='u', y_color_2='b',
                        x_color_int=0, y_color_int=0, av_val=1,
                        linewidth=2, line_color='k',
                        text=False, fontsize=20, text_color='k', x_text_offset=0.1, y_text_offset=-0.3):

    catalog_access = photometry_tools.data_access.CatalogAccess()

    nuv_wave = catalog_access.hst_wfc3_uvis1_bands_mean_wave['F275W']*1e-4
    u_wave = catalog_access.hst_wfc3_uvis1_bands_mean_wave['F336W']*1e-4
    b_wave = catalog_access.hst_wfc3_uvis1_bands_mean_wave['F438W']*1e-4
    v_wave = catalog_access.hst_wfc3_uvis1_bands_mean_wave['F555W']*1e-4
    i_wave = catalog_access.hst_wfc3_uvis1_bands_mean_wave['F814W']*1e-4

    x_wave_1 = locals()[x_color_1 + '_wave']
    x_wave_2 = locals()[x_color_2 + '_wave']
    y_wave_1 = locals()[y_color_1 + '_wave']
    y_wave_2 = locals()[y_color_2 + '_wave']

    color_ext_x = dust_tools.extinction_tools.ExtinctionTools.color_ext_ccm89_av(wave1=x_wave_1, wave2=x_wave_2, av=av_val)
    color_ext_y = dust_tools.extinction_tools.ExtinctionTools.color_ext_ccm89_av(wave1=y_wave_1, wave2=y_wave_2, av=av_val)

    slope_av_vector = ((y_color_int + color_ext_y) - y_color_int) / ((x_color_int + color_ext_x) - x_color_int)

    angle_av_vector = np.arctan(color_ext_y/color_ext_x) * 180/np.pi

    ax.annotate('', xy=(x_color_int + color_ext_x, y_color_int + color_ext_y), xycoords='data',
                xytext=(x_color_int, y_color_int), fontsize=fontsize,
                textcoords='data', arrowprops=dict(arrowstyle='-|>', color=line_color, lw=linewidth, ls='-'))

    if text:
        if isinstance(av_val, int):
            arrow_text = r'A$_{\rm V}$=%i mag' % av_val
        else:
            arrow_text = r'A$_{\rm V}$=%.1f mag' % av_val
        ax.text(x_color_int + x_text_offset, y_color_int + y_text_offset, arrow_text,
                horizontalalignment='left', verticalalignment='bottom',
                transform_rotates_text=True, rotation_mode='anchor',
                rotation=angle_av_vector, fontsize=fontsize, color=text_color)
