"""
class to gather all information need to access PHANGS observational data products
"""


class PhangsObsInfo:
    """
    all info are gathered in dictionaries as attributes
    """

    def __init__(self):
        super().__init__()

        # List of all Phangs galaxies.
        # This list is biased towards photometry observations and some ALMA observed galaxies are not incuded
        self.phangs_galaxy_list = ['ic1954', 'ic5332', 'ngc0628', 'ngc0685', 'ngc1087', 'ngc1097', 'ngc1300', 'ngc1317',
                                   'ngc1365', 'ngc1385', 'ngc1433', 'ngc1512', 'ngc1559', 'ngc1566', 'ngc1672',
                                   'ngc1792', 'ngc2775', 'ngc2835', 'ngc2903', 'ngc3351', 'ngc3621', 'ngc3627',
                                   'ngc4254', 'ngc4298', 'ngc4303', 'ngc4321', 'ngc4535', 'ngc4536', 'ngc4548',
                                   'ngc4569', 'ngc4571', 'ngc4654', 'ngc4689', 'ngc4826', 'ngc5068', 'ngc5248',
                                   'ngc6744', 'ngc7496']

        # Phangs target with existing HST observations
        self.phangs_hst_obs_target_list = ['ic1954', 'ic5332', 'ngc0628', 'ngc0628e', 'ngc0628c', 'ngc0685', 'ngc1087',
                                           'ngc1097', 'ngc1300', 'ngc1317', 'ngc1365', 'ngc1385', 'ngc1433', 'ngc1512',
                                           'ngc1559', 'ngc1566', 'ngc1672', 'ngc1792', 'ngc2775', 'ngc2835', 'ngc2903',
                                           'ngc3351', 'ngc3621', 'ngc3627', 'ngc4254', 'ngc4298', 'ngc4303', 'ngc4321',
                                           'ngc4535', 'ngc4536', 'ngc4548', 'ngc4569', 'ngc4571', 'ngc4654', 'ngc4689',
                                           'ngc4826', 'ngc5068', 'ngc5248', 'ngc6744', 'ngc7496']
        # Phangs target with existing NIRCAM observations
        self.phangs_nircam_obs_target_list = ['ngc0628', 'ngc1087', 'ngc1300', 'ngc1365', 'ngc1385', 'ngc1433',
                                              'ngc1512', 'ngc1566', 'ngc1672', 'ngc2835', 'ngc3351', 'ngc3627', 'ngc4254',
                                              'ngc4303', 'ngc4321', 'ngc4535', 'ngc5068', 'ngc7496']
        # Phangs target with existing MIRI observations
        self.phangs_miri_obs_target_list = ['ic5332', 'ngc0628', 'ngc1087', 'ngc1300', 'ngc1365', 'ngc1385', 'ngc1433',
                                            'ngc1512', 'ngc1566', 'ngc1672', 'ngc2835', 'ngc3351', 'ngc3627', 'ngc4254',
                                            'ngc4303', 'ngc4321', 'ngc4535', 'ngc5068', 'ngc7496']

        # # targets for which a cluster catalog was created
        # self.phangs_hst_cluster_cat_target_list = ['ic1954', 'ic5332', 'ngc0628e', 'ngc0628c', 'ngc0685', 'ngc1087',
        #                                            'ngc1097', 'ngc1300', 'ngc1317', 'ngc1365', 'ngc1385', 'ngc1433',
        #                                            'ngc1512', 'ngc1559', 'ngc1566', 'ngc1672', 'ngc1792', 'ngc2775',
        #                                            'ngc2835', 'ngc2903', 'ngc3351', 'ngc3621', 'ngc3627', 'ngc4254',
        #                                            'ngc4298', 'ngc4303', 'ngc4321', 'ngc4535', 'ngc4536', 'ngc4548',
        #                                            'ngc4569', 'ngc4571', 'ngc4654', 'ngc4689', 'ngc4826', 'ngc5068',
        #                                            'ngc5248', 'ngc6744', 'ngc7496']
        # targets for which a cluster catalog was created
        self.phangs_hst_cluster_cat_target_list = ['ngc1365']

        self.hst_ver_folder_names = {'v1.0': 'v1.0'}
        self.nircam_ver_folder_names = {'v0p9': 'v0p9', 'v0p9p1': 'v0p9p1'}
        self.miri_ver_folder_names = {'v0p9': 'v0p9', 'v0p9p1': 'v0p9p1'}

        # specification of observed bands for each HST target
        self.phangs_hst_obs_band_dict = {
            'ngc0628':
                {'folder_name': 'ngc628mosaic',
                 'acs_wfc1_observed_bands': ['F435W', 'F814W'],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F555W']},
            'ngc0628e':
                {'folder_name': 'ngc628mosaic',
                 'acs_wfc1_observed_bands': ['F435W', 'F814W'],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F555W']},
            'ngc0628c':
                {'folder_name': 'ngc628mosaic',
                 'acs_wfc1_observed_bands': ['F435W', 'F555W', 'F814W'],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W']},
            'ngc0685':
                {'folder_name': 'ngc685',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc1087':
                {'folder_name': 'ngc1087',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc1097':
                {'folder_name': 'ngc1097mosaic',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc1300':
                {'folder_name': 'ngc1300mosaic',
                 'acs_wfc1_observed_bands': ['F435W', 'F555W', 'F814W'],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W']},
            'ngc1317':
                {'folder_name': 'ngc1317',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ic1954':
                {'folder_name': 'ic1954',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc1365':
                {'folder_name': 'ngc1365',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc1385':
                {'folder_name': 'ngc1385',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc1433':
                {'folder_name': 'ngc1433',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc1512':
                {'folder_name': 'ngc1512mosaic',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc1559':
                {'folder_name': 'ngc1559',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc1566':
                {'folder_name': 'ngc1566',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc1672':
                {'folder_name': 'ngc1672mosaic',
                 'acs_wfc1_observed_bands': ['F555W', 'F435W', 'F814W'],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F555W']},
            'ngc1792':
                {'folder_name': 'ngc1792',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc2775':
                {'folder_name': 'ngc2775',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc2835':
                {'folder_name': 'ngc2835',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc2903':
                {'folder_name': 'ngc2903mosaic',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc3351':
                {'folder_name': 'ngc3351mosaic',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc3621':
                {'folder_name': 'ngc3621mosaic',
                 'acs_wfc1_observed_bands': ['F435W', 'F555W', 'F814W'],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W']},
            'ngc3627':
                {'folder_name': 'ngc3627mosaic',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc4254':
                {'folder_name': 'ngc4254mosaic',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc4298':
                {'folder_name': 'ngc4298',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc4303':
                {'folder_name': 'ngc4303',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc4321':
                {'folder_name': 'ngc4321mosaic',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc4535':
                {'folder_name': 'ngc4535',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc4536':
                {'folder_name': 'ngc4536mosaic',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc4548':
                {'folder_name': 'ngc4548',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc4569':
                {'folder_name': 'ngc4569',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc4571':
                {'folder_name': 'ngc4571',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc4654':
                {'folder_name': 'ngc4654',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc4689':
                {'folder_name': 'ngc4689',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc4826':
                {'folder_name': 'ngc4826',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc5068':
                {'folder_name': 'ngc5068mosaic',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc5248':
                {'folder_name': 'ngc5248',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc6744':
                {'folder_name': 'ngc6744mosaic',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ngc7496':
                {'folder_name': 'ngc7496',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
            'ic5332':
                {'folder_name': 'ic5332',
                 'acs_wfc1_observed_bands': [],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F438W', 'F555W', 'F814W']},
        }
        # specification of observed bands for each NIRCAM target
        self.nircam_targets = {
            'ngc0628': {'folder_name': 'ngc0628', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc1087': {'folder_name': 'ngc1087', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc1300': {'folder_name': 'ngc1300', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc1365': {'folder_name': 'ngc1365', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc1385': {'folder_name': 'ngc1385', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc1433': {'folder_name': 'ngc1433', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc1512': {'folder_name': 'ngc1512', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc1566': {'folder_name': 'ngc1566', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc1672': {'folder_name': 'ngc1672', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc2835': {'folder_name': 'ngc2835', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc3351': {'folder_name': 'ngc3351', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc3627': {'folder_name': 'ngc3627', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc4254': {'folder_name': 'ngc4254', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc4303': {'folder_name': 'ngc4303', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc4321': {'folder_name': 'ngc4321', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc4535': {'folder_name': 'ngc4535', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc5068': {'folder_name': 'ngc5068', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
            'ngc7496': {'folder_name': 'ngc7496', 'observed_bands': ['F200W', 'F300M', 'F335M', 'F360M']},
        }
        # specification of observed bands for each MIRI target
        self.miri_targets = {
            'ic5332': {'folder_name': 'ic5332', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc0628': {'folder_name': 'ngc0628', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc1087': {'folder_name': 'ngc1087', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc1300': {'folder_name': 'ngc1300', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc1365': {'folder_name': 'ngc1365', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc1385': {'folder_name': 'ngc1385', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc1433': {'folder_name': 'ngc1433', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc1512': {'folder_name': 'ngc1512', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc1566': {'folder_name': 'ngc1566', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc1672': {'folder_name': 'ngc1672', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc2835': {'folder_name': 'ngc2835', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc3351': {'folder_name': 'ngc3351', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc3627': {'folder_name': 'ngc3627', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc4254': {'folder_name': 'ngc4254', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc4303': {'folder_name': 'ngc4303', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc4321': {'folder_name': 'ngc4321', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc4535': {'folder_name': 'ngc4535', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc5068': {'folder_name': 'ngc5068', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
            'ngc7496': {'folder_name': 'ngc7496', 'observed_bands': ['F770W', 'F1000W', 'F1130W', 'F2100W']},
        }


class PhysParams:
    """
    Class to gather all physical params
    """

    def __init__(self):
        super().__init__()

        """
        distances need to be done!!!!! See Lee et al 2022 Table 1
        """

        self.sr_per_square_deg = 0.00030461741978671  # steradians per square degree

        # zero point NIRCAM flux corrections for data from the pipeline version v0p4p2
        self.nircam_zero_point_flux_corr = {'F200W': 0.854, 'F300M': 0.997, 'F335M': 1.000, 'F360M': 1.009}

        # filter names from http://svo2.cab.inta-csic.es
        self.hst_acs_wfc1_bands = ['FR388N', 'FR423N', 'F435W', 'FR459M', 'FR462N', 'F475W', 'F502N', 'FR505N', 'F555W',
                                   'FR551N', 'F550M', 'FR601N', 'F606W', 'F625W', 'FR647M', 'FR656N', 'F658N', 'F660N',
                                   'FR716N', 'POL_UV', 'POL_V', 'G800L', 'F775W', 'FR782N', 'F814W', 'FR853N', 'F892N',
                                   'FR914M', 'F850LP', 'FR931N', 'FR1016N']
        self.hst_wfc3_uvis2_bands = ['F218W', 'FQ232N', 'F225W', 'FQ243N', 'F275W', 'F280N', 'F300X', 'F336W', 'F343N',
                                     'F373N', 'FQ378N', 'FQ387N', 'F390M', 'F390W', 'F395N', 'F410M', 'FQ422M', 'F438W',
                                     'FQ436N', 'FQ437N', 'G280', 'F467M', 'F469N', 'F475W', 'F487N', 'FQ492N', 'F502N',
                                     'F475X', 'FQ508N', 'F555W', 'F547M', 'FQ575N', 'F606W', 'F200LP', 'FQ619N',
                                     'F621M', 'F625W', 'F631N', 'FQ634N', 'F645N', 'F350LP', 'F656N', 'F657N', 'F658N',
                                     'F665N', 'FQ672N', 'FQ674N', 'F673N', 'F680N', 'F689M', 'FQ727N', 'FQ750N',
                                     'F763M', 'F600LP', 'F775W', 'F814W', 'F845M', 'FQ889N', 'FQ906N', 'F850LP',
                                     'FQ924N', 'FQ937N', 'F953N']
        self.nircam_bands = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W', 'F162M', 'F164N', 'F150W2', 'F182M', 'F187N',
                             'F200W', 'F210M', 'F212N', 'F250M', 'F277W', 'F300M', 'F323N', 'F322W2', 'F335M', 'F356W',
                             'F360M', 'F405N', 'F410M', 'F430M', 'F444W', 'F460M', 'F466N', 'F470N', 'F480M']
        self.miri_bands = ['F560W', 'F770W', 'F1000W', 'F1065C', 'F1140C', 'F1130W', 'F1280W', 'F1500W', 'F1550C',
                           'F1800W', 'F2100W', 'F2300C', 'F2550W']

        # band wavelength taken from
        # http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=HST&gname2=ACS_WFC&asttype=
        self.hst_acs_wfc1_bands_mean_wave = {
            'FR388N': 3881.71,
            'FR423N': 4230.39,
            'F435W': 4360.06,
            'FR459M': 4592.76,
            'FR462N': 4620.13,
            'F475W': 4802.31,
            'F502N': 5023.13,
            'FR505N': 5050.39,
            'F555W': 5397.60,
            'FR551N': 5510.30,
            'F550M': 5588.24,
            'FR601N': 6010.50,
            'F606W': 6035.73,
            'F625W': 6352.46,
            'FR647M': 6476.15,
            'FR656N': 6560.36,
            'F658N': 6584.10,
            'F660N': 6599.50,
            'FR716N': 7160.04,
            'POL_UV': 7294.11,
            'POL_V': 7523.49,
            'G800L': 7704.08,
            'F775W': 7730.77,
            'FR782N': 7819.44,
            'F814W': 8129.21,
            'FR853N': 8528.80,
            'F892N': 8915.37,
            'FR914M': 9079.84,
            'F850LP': 9080.26,
            'FR931N': 9306.31,
            'FR1016N': 10150.22,
        }
        self.hst_wfc3_uvis1_bands_mean_wave = {
            'F218W': 2231.14,
            'FQ232N': 2327.12,
            'F225W': 2377.24,
            'FQ243N': 2420.59,
            'F275W': 2718.36,
            'F280N': 2796.98,
            'F300X': 2867.82,
            'F336W': 3365.86,
            'F343N': 3438.50,
            'F373N': 3730.19,
            'FQ378N': 3792.78,
            'FQ387N': 3873.61,
            'F390M': 3898.62,
            'F390W': 3952.50,
            'F395N': 3955.38,
            'F410M': 4109.81,
            'FQ422M': 4219.70,
            'F438W': 4338.57,
            'FQ436N': 4367.41,
            'FQ437N': 4371.30,
            'G280': 4628.43,
            'F467M': 4683.55,
            'F469N': 4688.29,
            'F475W': 4827.71,
            'F487N': 4871.54,
            'FQ492N': 4933.83,
            'F502N': 5009.93,
            'F475X': 5076.23,
            'FQ508N': 5091.59,
            'F555W': 5388.55,
            'F547M': 5459.04,
            'FQ575N': 5756.92,
            'F606W': 5999.27,
            'F200LP': 6043.00,
            'FQ619N': 6198.49,
            'F621M': 6227.39,
            'F625W': 6291.29,
            'F631N': 6304.27,
            'FQ634N': 6349.37,
            'F645N': 6453.59,
            'F350LP': 6508.00,
            'F656N': 6561.54,
            'F657N': 6566.93,
            'F658N': 6585.64,
            'F665N': 6656.23,
            'FQ672N': 6717.13,
            'FQ674N': 6730.58,
            'F673N': 6766.27,
            'F680N': 6880.13,
            'F689M': 6885.92,
            'FQ727N': 7275.84,
            'FQ750N': 7502.54,
            'F763M': 7623.09,
            'F600LP': 7656.67,
            'F775W': 7683.41,
            'F814W': 8117.36,
            'F845M': 8449.34,
            'FQ889N': 8892.56,
            'FQ906N': 9058.19,
            'F850LP': 9207.49,
            'FQ924N': 9247.91,
            'FQ937N': 9372.90,
            'F953N': 9531.11,
        }
        self.nircam_bands_mean_wave = {
            'F070W': 7088.30,
            'F090W': 9083.40,
            'F115W': 11623.89,
            'F140M': 14074.46,
            'F150W': 15104.23,
            'F162M': 16296.59,
            'F164N': 16445.95,
            'F150W2': 17865.58,
            'F182M': 18494.30,
            'F187N': 18739.65,
            'F200W': 20028.15,
            'F210M': 20982.22,
            'F212N': 21213.97,
            'F250M': 25049.39,
            'F277W': 27844.64,
            'F300M': 29940.44,
            'F323N': 32369.29,
            'F322W2': 33334.98,
            'F335M': 33675.24,
            'F356W': 35934.49,
            'F360M': 36298.10,
            'F405N': 40517.39,
            'F410M': 40886.55,
            'F430M': 42829.39,
            'F444W': 44393.50,
            'F460M': 46315.57,
            'F466N': 46545.31,
            'F470N': 47078.82,
            'F480M': 48213.27
        }
        self.miri_bands_mean_wave = {
            'F560W': 56651.28,
            'F770W': 77111.39,
            'F1000W': 99981.09,
            'F1065C': 105681.52,
            'F1140C': 113156.52,
            'F1130W': 113159.44,
            'F1280W': 128738.34,
            'F1500W': 151469.08,
            'F1550C': 155219.65,
            'F1800W': 180508.31,
            'F2100W': 209373.20,
            'F2300C': 227630.49,
            'F2550W': 254994.19
        }