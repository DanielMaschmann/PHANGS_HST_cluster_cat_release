"""
Construct a data access structure for HST and JWST imaging data
"""

from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from scipy.constants import c as speed_of_light

from cluster_cat_dr import helper_func
from cluster_cat_dr import phangs_info


class PhotAccess(phangs_info.PhangsObsInfo, phangs_info.PhysParams):
    """
    Access class to organize data structure of HST, NIRCAM and MIRI imaging data
    """

    def __init__(self, target_name=None, hst_data_path=None, nircam_data_path=None, miri_data_path=None,
                 hst_data_ver='v1.0', nircam_data_ver='v0p4p2', miri_data_ver='v0p5'):
        """
        Parameters
        ----------
        target_name : str
            Default None. Target name
        hst_data_path : str
            Default None. Path to HST imaging data
        nircam_data_path : str
            Default None. Path to NIRCAM imaging data
        miri_data_path : str
            Default None. Path to MIRI imaging data
        hst_data_ver : str
            Default None. denotes the internal data release version of the HST imaging data.
        nircam_data_ver : str
            Default None. denotes the internal data release version of the NIRCAM imaging data.
        miri_data_ver : str
            Default None. denotes the internal data release version of the MIRI imaging data.

        """
        super().__init__()

        self.target_name = target_name
        if (((self.target_name not in self.phangs_galaxy_list) &
            (self.target_name not in self.phangs_hst_obs_target_list) &
            (self.target_name not in self.phangs_nircam_obs_target_list) &
            (self.target_name not in self.phangs_miri_obs_target_list)) &
                (self.target_name is not None)):
            raise AttributeError('The target %s is not in the PHANGS photometric sample or has not been added to '
                                 'the current package version' % self.target_name)

        if hst_data_path is None:
            self.hst_data_path = None
        else:
            self.hst_data_path = Path(hst_data_path)
        if nircam_data_path is None:
            self.nircam_data_path = None
        else:
            self.nircam_data_path = Path(nircam_data_path)
        if miri_data_path is None:
            self.miri_data_path = None
        else:
            self.miri_data_path = Path(miri_data_path)

        self.hst_data_ver = hst_data_ver
        self.nircam_data_ver = nircam_data_ver
        self.miri_data_ver = miri_data_ver

        # loaded data dictionaries
        self.hst_bands_data = {}
        self.nircam_bands_data = {}
        self.miri_bands_data = {}

    def get_hst_img_file_name(self, band, hst_data_folder=None, file_name=None):
        """

        Parameters
        ----------
        band : str
        hst_data_folder : str
        file_name : str

        Returns
        -------
        data_file_path : Path
        """
        if (band not in self.hst_acs_wfc1_bands) & (band not in self.hst_wfc3_uvis2_bands):
            raise AttributeError('The band <%s> is not in the list of possible HST bands.' % band)

        if hst_data_folder is None:
            hst_data_folder = (self.hst_data_path / self.hst_ver_folder_names[self.hst_data_ver] /
                               self.phangs_hst_obs_band_dict[self.target_name]['folder_name'])
        if self.hst_data_ver == 'v1.0':
            hst_data_ver_str = 'v1'
        else:
            hst_data_ver_str = self.hst_data_ver
        ending_of_band_file_1 = '%s_%s_exp-drc-sci.fits' % (band.lower(), hst_data_ver_str)
        ending_of_band_file_2 = '%s_exp_drc_sci.fits' % (band.lower())

        if file_name is None:
            return helper_func.identify_file_in_folder(folder_path=hst_data_folder,
                                                       str_in_file_name_1=ending_of_band_file_1,
                                                       str_in_file_name_2=ending_of_band_file_2)
        else:
            return Path(hst_data_folder) / Path(file_name)

    def get_hst_exp_time_file_name(self, band, hst_data_folder=None, file_name=None):
        """

        Parameters
        ----------
        band : str
        hst_data_folder : str
        file_name : str

        Returns
        -------
        data_file_path : Path
        """
        if (band not in self.hst_acs_wfc1_bands) & (band not in self.hst_wfc3_uvis2_bands):
            raise AttributeError('The band <%s> is not in the list of possible HST bands.' % band)

        if hst_data_folder is None:
            hst_data_folder = (self.hst_data_path / self.hst_ver_folder_names[self.hst_data_ver] /
                               self.phangs_hst_obs_band_dict[self.target_name]['folder_name'])
        if self.hst_data_ver == 'v1.0':
            hst_data_ver_str = 'v1'
        else:
            hst_data_ver_str = self.hst_data_ver
        ending_of_band_file_1 = '%s_%s_exp-drc-wht.fits' % (band.lower(), hst_data_ver_str)
        ending_of_band_file_2 = '%s_exp_drc_wht.fits' % (band.lower())

        if file_name is None:
            return helper_func.identify_file_in_folder(folder_path=hst_data_folder,
                                                       str_in_file_name_1=ending_of_band_file_1,
                                                       str_in_file_name_2=ending_of_band_file_2)
        else:
            return Path(hst_data_folder) / Path(file_name)

    def get_hst_err_file_name(self, band, hst_data_folder=None, file_name=None):
        """

        Parameters
        ----------
        band : str
        hst_data_folder : str
        file_name : str

        Returns
        -------
        data_file_path : Path
        """
        if (band not in self.hst_acs_wfc1_bands) & (band not in self.hst_wfc3_uvis2_bands):
            raise AttributeError('The band <%s> is not in the list of possible HST bands.' % band)
        if hst_data_folder is None:
            hst_data_folder = (self.hst_data_path / self.hst_ver_folder_names[self.hst_data_ver] /
                               self.phangs_hst_obs_band_dict[self.target_name]['folder_name'])

        if self.hst_data_ver == 'v1.0':
            hst_data_ver_str = 'v1'
        else:
            hst_data_ver_str = self.hst_data_ver
        ending_of_band_file = '%s_%s_err-drc-wht.fits' % (band.lower(), hst_data_ver_str)

        if file_name is None:
            return helper_func.identify_file_in_folder(folder_path=hst_data_folder,
                                                       str_in_file_name_1=ending_of_band_file)
        else:
            return Path(hst_data_folder) / Path(file_name)

    def get_nircam_img_file_name(self, band, nircam_data_folder=None, file_name=None):
        """

        Parameters
        ----------
        band : str
        nircam_data_folder : str
        file_name : str

        Returns
        -------
        data_file_path : Path
        """
        if band not in self.nircam_bands:
            raise AttributeError('The band <%s> is not in the list of possible NIRCAM bands.' % band)

        if nircam_data_folder is None:
            nircam_data_folder = (self.nircam_data_path / self.nircam_ver_folder_names[self.nircam_data_ver] /
                                  self.nircam_targets[self.target_name]['folder_name'])

        ending_of_band_file = 'nircam_lv3_%s_i2d_align.fits' % band.lower()
        if file_name is None:
            return helper_func.identify_file_in_folder(folder_path=nircam_data_folder,
                                                       str_in_file_name_1=ending_of_band_file)
        else:
            return Path(nircam_data_folder) / Path(file_name)

    def get_miri_img_file_name(self, band, miri_data_folder=None, file_name=None):
        """

        Parameters
        ----------
        band : str
        miri_data_folder : str
        file_name : str


        Returns
        -------
        data_file_path : Path
        """
        if band not in self.miri_bands:
            raise AttributeError('The band <%s> is not in the list of possible MIRI bands.' % band)

        if miri_data_folder is None:
            miri_data_folder = (self.miri_data_path / self.miri_ver_folder_names[self.miri_data_ver] /
                                self.nircam_targets[self.target_name]['folder_name'])

        ending_of_band_file = '%s_miri_lv3_%s_i2d_align.fits' % (self.target_name, band.lower())

        if file_name is None:
            return helper_func.identify_file_in_folder(folder_path=miri_data_folder,
                                                       str_in_file_name_1=ending_of_band_file)
        else:
            return Path(miri_data_folder) / Path(file_name)

    def get_miri_err_file_name(self, band, miri_data_folder=None, file_name=None):
        """

        Parameters
        ----------
        band : str
        miri_data_folder : str
        file_name : str

        Returns
        -------
        data_file_path : Path
        """
        if band not in self.miri_bands:
            raise AttributeError('The band <%s> is not in the list of possible MIRI bands.' % band)

        if miri_data_folder is None:
            miri_data_folder = self.miri_data_path / self.miri_ver_folder_names[self.miri_data_ver]
        ending_of_band_file = '%s_miri_%s_noisemap.fits' % (self.target_name, band.lower())
        if file_name is None:
            return helper_func.identify_file_in_folder(folder_path=miri_data_folder,
                                                       str_in_file_name_1=ending_of_band_file)
        else:
            return Path(miri_data_folder) / Path(file_name)

    def load_hst_band(self, band, load_err=True, flux_unit='Jy', hst_data_folder=None, img_file_name=None,
                      err_file_name=None):
        """

        Parameters
        ----------
        band : str
        load_err : bool
        flux_unit : str
        hst_data_folder : str
        img_file_name : str
        err_file_name : str
        """
        # load the band observations
        img_file_name = self.get_hst_img_file_name(band=band, hst_data_folder=hst_data_folder, file_name=img_file_name)
        img_data, img_header, img_wcs = helper_func.load_img(file_name=img_file_name)

        # convert the flux unit
        if 'PHOTFNU' in img_header:
            conversion_factor = img_header['PHOTFNU']
        elif 'PHOTFLAM' in img_header:
            # wavelength in angstrom
            pivot_wavelength = img_header['PHOTPLAM']
            # inverse sensitivity, ergs/cm2/Ang/electron
            sensitivity = img_header['PHOTFLAM']
            # speed of light in Angstrom/s
            c = speed_of_light * 1e10
            # change the conversion facto to get erg s−1 cm−2 Hz−1
            f_nu = sensitivity * pivot_wavelength ** 2 / c
            # change to get Jy
            conversion_factor = f_nu * 1e23
        else:
            raise KeyError('there is no PHOTFNU or PHOTFLAM in the header')

        pixel_area_size_sr = img_wcs.proj_plane_pixel_area().value * self.sr_per_square_deg
        # rescale data image
        if flux_unit == 'Jy':
            # rescale to Jy
            conversion_factor = conversion_factor
        elif flux_unit == 'mJy':
            # rescale to Jy
            conversion_factor *= 1e3
        elif flux_unit == 'MJy/sr':
            # get the size of one pixel in sr with the factor 1e6 for the conversion of Jy to MJy later
            # change to MJy/sr
            conversion_factor /= (pixel_area_size_sr * 1e6)
        else:
            raise KeyError('flux_unit ', flux_unit, ' not understand!')

        img_data *= conversion_factor
        self.hst_bands_data.update({'%s_data_img' % band: img_data, '%s_header_img' % band: img_header,
                                    '%s_wcs_img' % band: img_wcs, '%s_unit_img' % band: flux_unit,
                                    '%s_pixel_area_size_sr_img' % band: pixel_area_size_sr})
        if load_err:
            err_file_name = self.get_hst_err_file_name(band=band, hst_data_folder=hst_data_folder,
                                                       file_name=err_file_name)
            err_data, err_header, err_wcs = helper_func.load_img(file_name=err_file_name)
            err_data = 1 / np.sqrt(err_data)
            err_data *= conversion_factor

            self.hst_bands_data.update({'%s_data_err' % band: err_data, '%s_header_err' % band: err_header,
                                        '%s_wcs_err' % band: err_wcs, '%s_unit_err' % band: flux_unit,
                                        '%s_pixel_area_size_sr_err' % band: pixel_area_size_sr})

    def load_nircam_band(self, band, load_err=True, flux_unit='Jy', nircam_data_folder=None, img_file_name=None):
        """

        Parameters
        ----------
        band : str
        load_err : bool
        flux_unit : str
        nircam_data_folder : str
        img_file_name : str

        """
        # load the band observations
        file_name = self.get_nircam_img_file_name(band=band, nircam_data_folder=nircam_data_folder,
                                                  file_name=img_file_name)
        img_data, img_header, img_wcs = helper_func.load_img(file_name=file_name, hdu_number='SCI')
        pixel_area_size_sr = img_wcs.proj_plane_pixel_area().value * self.sr_per_square_deg
        # rescale data image
        if flux_unit == 'Jy':
            # rescale to Jy
            conversion_factor = pixel_area_size_sr * 1e6

        elif flux_unit == 'mJy':
            # rescale to Jy
            conversion_factor = pixel_area_size_sr * 1e9
        elif flux_unit == 'MJy/sr':
            conversion_factor = 1
        else:
            raise KeyError('flux_unit ', flux_unit, ' not understand')

        img_data *= conversion_factor
        self.nircam_bands_data.update({'%s_data_img' % band: img_data, '%s_header_img' % band: img_header,
                                       '%s_wcs_img' % band: img_wcs, '%s_unit_img' % band: flux_unit,
                                       '%s_pixel_area_size_sr_img' % band: pixel_area_size_sr})
        if load_err:
            err_data, err_header, err_wcs = helper_func.load_img(file_name=file_name, hdu_number='ERR')
            err_data *= conversion_factor
            self.nircam_bands_data.update({'%s_data_err' % band: err_data, '%s_header_err' % band: err_header,
                                           '%s_wcs_err' % band: img_wcs, '%s_unit_err' % band: flux_unit,
                                           '%s_pixel_area_size_sr_err' % band: pixel_area_size_sr})

    def load_miri_band(self, band, load_err=True, flux_unit='Jy', miri_data_folder=None, img_file_name=None,
                       err_file_name=None):
        """

        Parameters
        ----------
        band : str
        load_err : bool
        flux_unit : str
        miri_data_folder : str
        img_file_name : str
        err_file_name : str
        """
        # load the band observations
        file_name = self.get_miri_img_file_name(band=band, miri_data_folder=miri_data_folder, file_name=img_file_name)
        # if self.miri_data_ver == 'v0p6':
        hdu_number = 'SCI'
        # else:
        #     hdu_number = 0
        img_data, img_header, img_wcs = helper_func.load_img(file_name=file_name, hdu_number=hdu_number)
        pixel_area_size_sr = img_wcs.proj_plane_pixel_area().value * self.sr_per_square_deg
        # rescale data image
        if flux_unit == 'Jy':
            # rescale to Jy
            conversion_factor = pixel_area_size_sr * 1e6

        elif flux_unit == 'mJy':
            # rescale to Jy
            conversion_factor = pixel_area_size_sr * 1e9
        elif flux_unit == 'MJy/sr':
            conversion_factor = 1
        else:
            raise KeyError('flux_unit ', flux_unit, ' not understand')

        img_data *= conversion_factor
        self.miri_bands_data.update({'%s_data_img' % band: img_data, '%s_header_img' % band: img_header,
                                     '%s_wcs_img' % band: img_wcs, '%s_unit_img' % band: flux_unit,
                                     '%s_pixel_area_size_sr_img' % band: pixel_area_size_sr})
        if load_err:
            err_data, err_header, err_wcs = helper_func.load_img(file_name=file_name, hdu_number='ERR')
            err_data *= conversion_factor
            self.miri_bands_data.update({'%s_data_err' % band: err_data, '%s_header_err' % band: err_header,
                                         '%s_wcs_err' % band: img_wcs, '%s_unit_err' % band: flux_unit,
                                         '%s_pixel_area_size_sr_err' % band: pixel_area_size_sr})

    def get_hst_band_list(self):
        """
        gets list of bands of HST
        Returns
        -------
        band_list : list
        """
        band_list = []
        if self.target_name in self.phangs_hst_obs_target_list:
            for band in list(set(self.hst_acs_wfc1_bands + self.hst_wfc3_uvis2_bands)):
                if band in (self.phangs_hst_obs_band_dict[self.target_name]['acs_wfc1_observed_bands'] +
                            self.phangs_hst_obs_band_dict[self.target_name]['wfc3_uvis_observed_bands']):
                    band_list.append(band)
        return self.sort_band_list(band_list=band_list)

    def get_nircam_band_list(self):
        """
        gets list of bands of NIRCAM
        Returns
        -------
        band_list : list
        """
        band_list = []
        if self.target_name in self.phangs_nircam_obs_target_list:
            for band in self.nircam_bands:
                if band in self.nircam_targets[self.target_name]['observed_bands']:
                    band_list.append(band)
            band_list = self.sort_band_list(band_list=band_list)
        return band_list

    def get_miri_band_list(self):
        """
        gets list of bands of NIRCAM
        Returns
        -------
        band_list : list
        """
        band_list = []
        if self.target_name in self.phangs_miri_obs_target_list:
            for band in self.miri_bands:
                if band in self.miri_targets[self.target_name]['observed_bands']:
                    band_list.append(band)
            band_list = self.sort_band_list(band_list=band_list)
        return band_list

    def get_hst_nircam_miri_band_list(self):
        """
        wrapper to load all bands observed with HST, NIRCAM and MIRI for one target
        Returns
        -------
        band_list : list
        """
        band_list = self.get_hst_band_list() + self.get_nircam_band_list() + self.get_miri_band_list()

        return self.sort_band_list(band_list=band_list)

    def load_hst_nircam_miri_bands(self, band_list=None, flux_unit='Jy', folder_name_list=None,
                                   img_file_name_list=None, err_file_name_list=None,
                                   load_err=True):
        """
        wrapper to load all available HST, NIRCAM and MIRI observations into the constructor

        Parameters
        ----------
        band_list : list or str
        flux_unit : str
        folder_name_list : list
        img_file_name_list : list
        err_file_name_list : list
        load_err : bool
        """
        # geta list with all observed bands in order of wavelength
        if band_list is None:
            band_list = self.get_hst_nircam_miri_band_list()
        elif isinstance(band_list, str):
            band_list = [band_list]

        if folder_name_list is None:
            folder_name_list = [None] * len(band_list)
        elif isinstance(folder_name_list, str):
            folder_name_list = [folder_name_list]
        if img_file_name_list is None:
            img_file_name_list = [None] * len(band_list)
        elif isinstance(img_file_name_list, str):
            img_file_name_list = [img_file_name_list]
        if err_file_name_list is None:
            err_file_name_list = [None] * len(band_list)
        elif isinstance(err_file_name_list, str):
            err_file_name_list = [err_file_name_list]

        # load bands
        for band, folder_name, img_file_name, err_file_name in zip(band_list, folder_name_list, img_file_name_list,
                                                                   err_file_name_list):
            if band in list(set(self.hst_acs_wfc1_bands + self.hst_wfc3_uvis2_bands)):
                self.load_hst_band(band=band, flux_unit=flux_unit, hst_data_folder=folder_name,
                                   img_file_name=img_file_name, err_file_name=err_file_name, load_err=load_err)

            elif band in self.nircam_bands:
                self.load_nircam_band(band=band, flux_unit=flux_unit, nircam_data_folder=folder_name,
                                      img_file_name=img_file_name, load_err=load_err)
            elif band in self.miri_bands:
                self.load_miri_band(band=band, flux_unit=flux_unit, miri_data_folder=folder_name,
                                    img_file_name=img_file_name, err_file_name=err_file_name, load_err=load_err)
            else:
                raise KeyError('Band is not found in possible band lists')

    def change_hst_nircam_miri_band_units(self, band_list=None, new_unit='MJy/sr'):
        """

        Parameters
        ----------
        band_list : list
        new_unit : str
        """
        if band_list is None:
            band_list = []
            for band in list(set(self.hst_acs_wfc1_bands + self.hst_wfc3_uvis2_bands)):
                if band in (self.phangs_hst_obs_band_dict[self.target_name]['acs_wfc1_observed_bands'] +
                            self.phangs_hst_obs_band_dict[self.target_name]['wfc3_uvis_observed_bands']):
                    band_list.append(band)
            for band in self.nircam_bands:
                if band in self.nircam_targets[self.target_name]['observed_bands']:
                    band_list.append(band)
            for band in self.miri_bands:
                if band in self.miri_targets[self.target_name]['observed_bands']:
                    band_list.append(band)

        for band in band_list:
            self.change_band_unit(band=band, new_unit=new_unit)

    def change_band_unit(self, band, new_unit='MJy/sr'):
        """
        will change loaded data to the needed unit
        Parameters
        ----------
        band : str
        new_unit : str
        """
        if band in list(set(self.hst_acs_wfc1_bands + self.hst_wfc3_uvis2_bands)):

            old_unit = self.hst_bands_data['%s_unit_img' % band]
            conversion_factor = 1
            # change to Jy
            if 'm' in old_unit:
                conversion_factor *= 1e-3
            if 'M' in old_unit:
                conversion_factor *= 1e6
            if '/sr' in old_unit:
                conversion_factor *= self.hst_bands_data['%s_pixel_area_size_sr_img' % band]

            # change to the new unit
            if 'm' in new_unit:
                conversion_factor *= 1e3
            if 'M' in new_unit:
                conversion_factor *= 1e-6
            if '/sr' in new_unit:
                conversion_factor /= self.hst_bands_data['%s_pixel_area_size_sr_img' % band]

            self.hst_bands_data['%s_data_img' % band] *= conversion_factor
            self.hst_bands_data['%s_unit_img' % band] = new_unit
            if '%s_data_err' % band in self.hst_bands_data.keys():
                self.hst_bands_data['%s_data_err' % band] *= conversion_factor
                self.hst_bands_data['%s_unit_err' % band] = new_unit

        if band in self.nircam_bands:
            old_unit = self.nircam_bands_data['%s_unit_img' % band]

            conversion_factor = 1
            # change to Jy
            if 'm' in old_unit:
                conversion_factor *= 1e-3
            if 'M' in old_unit:
                conversion_factor *= 1e6
            if '/sr' in old_unit:
                conversion_factor *= self.nircam_bands_data['%s_pixel_area_size_sr_img' % band]

            # change to the new unit
            if 'm' in new_unit:
                conversion_factor *= 1e3
            if 'M' in new_unit:
                conversion_factor *= 1e-6
            if '/sr' in new_unit:
                conversion_factor /= self.nircam_bands_data['%s_pixel_area_size_sr_img' % band]

            self.nircam_bands_data['%s_data_img' % band] *= conversion_factor
            self.nircam_bands_data['%s_unit_img' % band] = new_unit

            if '%s_data_err' % band in self.nircam_bands_data.keys():
                self.nircam_bands_data['%s_data_err' % band] *= conversion_factor
                self.nircam_bands_data['%s_unit_err' % band] = new_unit

        if band in self.miri_bands:
            old_unit = self.miri_bands_data['%s_unit_img' % band]

            conversion_factor = 1
            # change to Jy
            if 'm' in old_unit:
                conversion_factor *= 1e-3
            if 'M' in old_unit:
                conversion_factor *= 1e6
            if '/sr' in old_unit:
                conversion_factor *= self.miri_bands_data['%s_pixel_area_size_sr_img' % band]

            # change to the new unit
            if 'm' in new_unit:
                conversion_factor *= 1e3
            if 'M' in new_unit:
                conversion_factor *= 1e-6
            if '/sr' in new_unit:
                conversion_factor /= self.miri_bands_data['%s_pixel_area_size_sr_img' % band]

            self.miri_bands_data['%s_data_img' % band] *= conversion_factor
            self.miri_bands_data['%s_unit_img' % band] = new_unit
            if '%s_data_err' % band in self.miri_bands_data.keys():
                self.miri_bands_data['%s_data_err' % band] *= conversion_factor
                self.miri_bands_data['%s_unit_err' % band] = new_unit

    def get_band_cutout_dict(self, ra_cutout, dec_cutout, cutout_size, include_err=False, band_list=None):
        """

        Parameters
        ----------
        ra_cutout : float
        dec_cutout : float
        cutout_size : float, tuple or list
            Units in arcsec. Cutout size of a box cutout. If float it will be used for both box length.
        include_err : bool
        band_list : list

        Returns
        -------
        cutout_dict : dict
        each element in dictionary is of type astropy.nddata.Cutout2D object
        """
        # geta list with all observed bands in order of wavelength
        if band_list is None:
            band_list = self.get_hst_nircam_miri_band_list()

        if not isinstance(cutout_size, list):
            cutout_size = [cutout_size] * len(band_list)

        cutout_pos = SkyCoord(ra=ra_cutout, dec=dec_cutout, unit=(u.degree, u.degree), frame='fk5')
        cutout_dict = {'cutout_pos': cutout_pos}
        cutout_dict.update({'band_list': band_list})

        for band, band_index in zip(band_list, range(len(band_list))):
            if band in list(set(self.hst_acs_wfc1_bands + self.hst_wfc3_uvis2_bands)):
                cutout_dict.update({
                    '%s_img_cutout' % band:
                        helper_func.get_img_cutout(img=self.hst_bands_data['%s_data_img' % band],
                                                   wcs=self.hst_bands_data['%s_wcs_img' % band],
                                                   coord=cutout_pos, cutout_size=cutout_size[band_index])})
                if include_err:
                    cutout_dict.update({
                        '%s_err_cutout' % band:
                            helper_func.get_img_cutout(img=self.hst_bands_data['%s_data_err' % band],
                                                       wcs=self.hst_bands_data['%s_wcs_err' % band],
                                                       coord=cutout_pos, cutout_size=cutout_size[band_index])})

            elif band in self.nircam_bands:
                cutout_dict.update({
                    '%s_img_cutout' % band:
                        helper_func.get_img_cutout(img=self.nircam_bands_data['%s_data_img' % band],
                                                   wcs=self.nircam_bands_data['%s_wcs_img' % band],
                                                   coord=cutout_pos, cutout_size=cutout_size[band_index])})
                if include_err:
                    cutout_dict.update({
                        '%s_err_cutout' % band:
                            helper_func.get_img_cutout(img=self.nircam_bands_data['%s_data_err' % band],
                                                       wcs=self.nircam_bands_data['%s_wcs_err' % band],
                                                       coord=cutout_pos, cutout_size=cutout_size[band_index])})

            elif band in self.miri_bands:
                cutout_dict.update({
                    '%s_img_cutout' % band:
                        helper_func.get_img_cutout(img=self.miri_bands_data['%s_data_img' % band],
                                                   wcs=self.miri_bands_data['%s_wcs_img' % band],
                                                   coord=cutout_pos, cutout_size=cutout_size[band_index])})
                if include_err:
                    cutout_dict.update({
                        '%s_err_cutout' % band:
                            helper_func.get_img_cutout(img=self.miri_bands_data['%s_data_err' % band],
                                                       wcs=self.miri_bands_data['%s_wcs_err' % band],
                                                       coord=cutout_pos, cutout_size=cutout_size[band_index])})

        return cutout_dict

    def get_band_wave(self, band, unit='mu'):
        """
        Returns mean wavelength of a specific band
        Parameters
        ----------
        unit : str
        band : str

        Returns
        -------
        wavelength : float
        """
        if band in self.phangs_hst_obs_band_dict[self.target_name]['acs_wfc1_observed_bands']:
            wave = self.hst_acs_wfc1_bands_mean_wave[band]
        elif band in self.phangs_hst_obs_band_dict[self.target_name]['wfc3_uvis_observed_bands']:
            wave = self.hst_wfc3_uvis1_bands_mean_wave[band]
        elif band in self.nircam_targets[self.target_name]['observed_bands']:
            wave = self.nircam_bands_mean_wave[band]
        elif band in self.miri_targets[self.target_name]['observed_bands']:
            wave = self.miri_bands_mean_wave[band]
        else:
            raise KeyError(band, 'not understand')

        if unit == 'angstrom':
            return wave
        if unit == 'nano':
            return wave * 1e-1
        elif unit == 'mu':
            return wave * 1e-4
        else:
            raise KeyError('return unit not understand')

    def sort_band_list(self, band_list):
        """
        sorts a band list with increasing wavelength
        Parameters
        ----------
        band_list : list

        Returns
        -------
        sorted_band_list : list
        """
        wave_list = []
        for band in band_list:
            wave_list.append(self.get_band_wave(band=band))

        # sort wavelength bands
        sort = np.argsort(wave_list)
        return list(np.array(band_list)[sort])

    def get_hst_median_exp_time(self, band):

        exp_file_name = self.get_hst_exp_time_file_name(band=band)
        # print(exp_file_name)
        data, header, wcs = helper_func.load_img(file_name=exp_file_name)
        return np.nanmedian(data[data != 0])
