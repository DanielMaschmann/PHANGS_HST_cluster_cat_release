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
