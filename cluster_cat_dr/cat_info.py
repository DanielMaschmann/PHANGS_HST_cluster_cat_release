"""
Script to gather all information going into data-description tables, documentation & released catalogs
"""
import astropy.units as u
import numpy as np


class CatalogInfo:
    """
    Class to gather all catalog information as attributes
    """
    def __init__(self):

        # PHANGS-HST galaxy sample list
        self.phangs_hst_target_list = ['ic1954', 'ic5332', 'ngc628e', 'ngc628c', 'ngc685', 'ngc1087', 'ngc1097',
                                       'ngc1300', 'ngc1317', 'ngc1365', 'ngc1385', 'ngc1433', 'ngc1512', 'ngc1559',
                                       'ngc1566', 'ngc1672', 'ngc1792', 'ngc2775', 'ngc2835', 'ngc2903', 'ngc3351',
                                       'ngc3621', 'ngc3627', 'ngc4254', 'ngc4298', 'ngc4303', 'ngc4321', 'ngc4535',
                                       'ngc4536', 'ngc4548', 'ngc4569', 'ngc4571', 'ngc4654', 'ngc4689', 'ngc4826',
                                       'ngc5068', 'ngc5248', 'ngc6744', 'ngc7496']
        # specification of observed bands for each HST target
        self.phangs_hst_obs_band_dict = {
            'ngc628':
                {'folder_name': 'ngc628mosaic',
                 'acs_wfc1_observed_bands': ['F435W', 'F814W'],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F555W']},
            'ngc628e':
                {'folder_name': 'ngc628mosaic',
                 'acs_wfc1_observed_bands': ['F435W', 'F814W'],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W', 'F555W']},
            'ngc628c':
                {'folder_name': 'ngc628mosaic',
                 'acs_wfc1_observed_bands': ['F435W', 'F555W', 'F814W'],
                 'wfc3_uvis_observed_bands': ['F275W', 'F336W']},
            'ngc685':
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

        # Catalog info is a dictionary providing the final column name for the fits files of the final data release.
        # The keys are sorted for the column names used in the final internal release.
        # 'col_name' provides the name we will use in the articles, documentations and DR files.
        # 'units' provides the units in astropy.units format for example seconds would be units.s
        # 'units_str' provides the unit in string format which we will provide in the documentation.
        # 'doc_comment' is the description for the documentation.
        # 'tab_comment' is the description going into the table in the article.
        self.cat_info = {

            # cluster identifier
            'INDEX': {'col_name': 'index',
                      'unit': None,
                      'unit_str': 'int',
                      'doc_comment': 'A running index from 1 to N, where N is the total number of objects in the '
                                     'catalog. Objects sorted in order of increasing Y pixel value on image.',
                      'tab_comment': 'Running index from 1 to N for each individual target'},
            'ID_PHANGS_CLUSTERS_v1p2': {'col_name': 'phangs_cluster_id',
                                        'unit': None,
                                        'unit_str': 'int',
                                        'doc_comment': 'Running PHANGS cluster ID. These values correspond to the '
                                                       'candidate catalogs produced by the PHANGS-HST team. Objects '
                                                       'sorted in order of increasing Y pixel value on image.',
                                        'tab_comment': 'Running PHANGS cluster ID for candidate catalog '
                                                       'cross-identification.'
                                        },
            'ID_PHANGS_ALLSOURCES': {'col_name': 'phangs_all_cluster_id',
                                     'unit': None,
                                     'unit_str': 'int',
                                     'doc_comment': 'Running PHANGS cluster ID. These values correspond to the '
                                                    'initial detection catalogs produced by the PHANGS-HST team. '
                                                    'Objects sorted in order of increasing Y pixel value on image.',
                                     'tab_comment': 'Running PHANGS cluster ID for initial source detection catalog '
                                                    'cross-identification.'
                                     },

            # coordinates of each cluster
            'PHANGS_X': {'col_name': 'x_pix_coord',
                         'unit': None,
                         'unit_str': 'pix',
                         'doc_comment': 'X pixel coordinate on HST image. Scale = 0.03962 arcsec/pixel. '
                                        'Numeration goes from 0 to n-1, where n is the number of Pixels on the X-axis. '
                                        'The value of 1 must be added to some other systems which star counting at 1 '
                                        '(e.g., IRAF, ds9).',
                         'tab_comment': 'X coordinates on HST X-pixel grid (0...n-1). Scale = 0.03962 arcsec/pixel.'
                         },
            'PHANGS_Y': {'col_name': 'y_pix_coord',
                         'unit': None,
                         'unit_str': 'pix',
                         'doc_comment': 'Y pixel coordinate on HST image. Scale = 0.03962 arcsec/pixel. '
                                        'Numeration goes from 0 to n-1, where n is the number of Pixels on the Y-axis. '
                                        'The value of 1 must be added to some other systems which star counting at 1 '
                                        '(e.g., IRAF, ds9).',
                         'tab_comment': 'Y coordinates on HST Y-pixel grid (0...n-1). Scale = 0.03962 arcsec/pixel.'
                         },

            'PHANGS_RA': {'col_name': 'ra',
                          'unit': u.deg,
                          'unit_str': 'deg',
                          'doc_comment': 'J2000 Right ascension, ICRS frame, calibrated against selected Gaia sources.',
                          'tab_comment': 'J2000 Right ascension, ICRS frame, calibrated against selected Gaia sources.'
                          },
            'PHANGS_DEC': {'col_name': 'dec',
                           'unit': u.deg,
                           'unit_str': 'deg',
                           'doc_comment': 'J2000 Declination, ICRS frame, calibrated against selected Gaia sources.',
                           'tab_comment': 'J2000 Declination, ICRS frame, calibrated against selected Gaia sources.'
                           },

            # Morphological Classification
            'PHANGS_CLUSTER_CLASS_HUMAN': {'col_name': 'cluster_class_human',
                                           'unit': None,
                                           'unit_str': 'int',
                                           'doc_comment': 'Human visual classification determined by coauthor '
                                                          'Brad Whitmore (BCW). Further details in '
                                                          'Whitmore et al. (2021) (2021MNRAS.506.5294W).'
                                                          'The classification is encoded in integers: '
                                                          '1 and 2 for class 1 and 2, respectively and '
                                                          '3 for class 3 compact associations. ',
                                           'tab_comment': 'Cluster class assigned through visual inspection. '
                                                          'Integers 1, 2 and 3 denotes class 1, 2 and '
                                                          'class 3 compact associations, respectively.'
                                           },
            'PHANGS_CLUSTER_CLASS_ML_VGG': {'col_name': 'cluster_class_ml',
                                            'unit': None,
                                            'unit_str': 'int',
                                            'doc_comment': 'Classification determined by VGG neural network model of '
                                                           'Hannon et al. (2023, in review) (2023MNRAS.tmp.2242H); '
                                                           'This model was trained using human classification of '
                                                           'PHANGS-HST cluster classes provided in '
                                                           '`cluster_class_human\'. '
                                                           'Further details are provided in '
                                                           'Wei et al. (2020) (2020MNRAS.493.3178W), '
                                                           'Whitmore et al. (2021) (2021MNRAS.506.5294W) and '
                                                           'Thilker et al. (2022) (2022MNRAS.509.4094T). '
                                                           'The cluster class was determined from 10 randomly '
                                                           'initialized models. The values are the same as for '
                                                           '`cluster_class_human\': Integers 1, 2 and 3 denotes '
                                                           'class 1, 2 and class 3 compact associations, respectively.',
                                            'tab_comment': 'Cluster class determined by VGG neural network. '
                                                           'Integer 1, 2 and 3 denotes class 1, 2 and '
                                                           'class 3 compact associations, respectively.'
                                            },
            'PHANGS_CLUSTER_CLASS_ML_VGG_QUAL': {'col_name': 'cluster_class_ml_qual',
                                                 'unit': None,
                                                 'unit_str': 'float',
                                                 'doc_comment': 'Quality value providing accuracy of the VGG neural '
                                                                'network classification provided in column '
                                                                '`cluster_class_ml\'. '
                                                                'The provided value lies between 0.3 and 1 and denotes '
                                                                'the frequency of the classification among the 10 '
                                                                'randomly initialized models.',
                                                 'tab_comment': 'Quality value for `cluster_class_ml\' '
                                                                'with values between 0.3 and 1, providing the '
                                                                'frequency of the mode among the 10 randomly '
                                                                'initialized models.'
                                                 },

            # Photometry
            # Vega magnitudes
            'PHANGS_F275W_VEGA_TOT': {'col_name': 'f275w_vega_tot',
                                      'unit': u.mag,
                                      'unit_str': 'mag',
                                      'doc_comment': 'WFC3 f275w (NUV-band) total vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. For more details on the '
                                                     'aperture correction see '
                                                     'Deger et al. (2022) (2022MNRAS.510...32D). Set to -9999 if '
                                                     'source is not covered by HST filter. See also the '
                                                     '`no_coverage_flag\' column.',
                                      'tab_comment': 'WFC3 f275w (NUV-band) total vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. Set to -9999 if source is not '
                                                     'covered by HST filter.'
                                      },
            'PHANGS_F275W_VEGA_TOT_ERR': {'col_name': 'f275w_vega_tot_err',
                                          'unit': u.mag,
                                          'unit_str': 'mag',
                                          'doc_comment': 'Uncertainty of `f275w_vega_tot\'',
                                          'tab_comment': 'Uncertainty of `f275w_vega_tot\''
                                          },
            'PHANGS_F336W_VEGA_TOT': {'col_name': 'f336w_vega_tot',
                                      'unit': u.mag,
                                      'unit_str': 'mag',
                                      'doc_comment': 'WFC3 f336w (U-band) total vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. For more details on the '
                                                     'aperture correction see '
                                                     'Deger et al. (2022) (2022MNRAS.510...32D). Set to -9999 if '
                                                     'source is not covered by HST filter. See also the '
                                                     '`no_coverage_flag\' column.',
                                      'tab_comment': None
                                      },
            'PHANGS_F336W_VEGA_TOT_ERR': {'col_name': 'f336w_vega_tot_err',
                                          'unit': u.mag,
                                          'unit_str': 'mag',
                                          'doc_comment': 'Uncertainty of `f336w_vega_tot\'',
                                          'tab_comment': None
                                          },
            'PHANGS_F435W_VEGA_TOT': {'col_name': 'f435w_vega_tot',
                                      'unit': u.mag,
                                      'unit_str': 'mag',
                                      'doc_comment': 'Note that for Targets observed with the UVIS detector, '
                                                     'the filter name of the B-band is f438w. '
                                                     'WFC3 f435w (B-band) total vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. For more details on the '
                                                     'aperture correction see '
                                                     'Deger et al. (2022) (2022MNRAS.510...32D). Set to -9999 if '
                                                     'source is not covered by HST filter. See also the '
                                                     '`no_coverage_flag\' column.',
                                      'tab_comment': None
                                      },
            'PHANGS_F435W_VEGA_TOT_ERR': {'col_name': 'f435w_vega_tot_err',
                                          'unit': u.mag,
                                          'unit_str': 'mag',
                                          'doc_comment': 'Uncertainty of `f435w_vega_tot\'',
                                          'tab_comment': None
                                          },
            'PHANGS_F438W_VEGA_TOT': {'col_name': 'f438w_vega_tot',
                                      'unit': u.mag,
                                      'unit_str': 'mag',
                                      'doc_comment': 'Note that for Targets observed with the WFC detector, '
                                                     'the filter name of the B-band is f435w. '
                                                     'WFC3 f438w (B-band) total vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. For more details on the '
                                                     'aperture correction see '
                                                     'Deger et al. (2022) (2022MNRAS.510...32D). Set to -9999 if '
                                                     'source is not covered by HST filter. See also the '
                                                     '`no_coverage_flag\' column.',
                                      'tab_comment': None
                                      },
            'PHANGS_F438W_VEGA_TOT_ERR': {'col_name': 'f438w_vega_tot_err',
                                          'unit': u.mag,
                                          'unit_str': 'mag',
                                          'doc_comment': 'Uncertainty of `f438w_vega_tot\'',
                                          'tab_comment': None
                                          },
            'PHANGS_F555W_VEGA_TOT': {'col_name': 'f555w_vega_tot',
                                      'unit': u.mag,
                                      'unit_str': 'mag',
                                      'doc_comment': 'WFC3 f555w (V-band) total vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. For more details on the '
                                                     'aperture correction see '
                                                     'Deger et al. (2022) (2022MNRAS.510...32D). Set to -9999 if '
                                                     'source is not covered by HST filter. See also the '
                                                     '`no_coverage_flag\' column.',
                                      'tab_comment': None
                                      },
            'PHANGS_F555W_VEGA_TOT_ERR': {'col_name': 'f555w_vega_tot_err',
                                          'unit': u.mag,
                                          'unit_str': 'mag',
                                          'doc_comment': 'Uncertainty of `f555w_vega_tot\'',
                                          'tab_comment': None
                                          },
            'PHANGS_F814W_VEGA_TOT': {'col_name': 'f814w_vega_tot',
                                      'unit': u.mag,
                                      'unit_str': 'mag',
                                      'doc_comment': 'WFC3 f814w (I-band) total vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. For more details on the '
                                                     'aperture correction see '
                                                     'Deger et al. (2022) (2022MNRAS.510...32D). Set to -9999 if '
                                                     'source is not covered by HST filter. See also the '
                                                     '`no_coverage_flag\' column.',
                                      'tab_comment': None
                                      },
            'PHANGS_F814W_VEGA_TOT_ERR': {'col_name': 'f814w_vega_tot_err',
                                          'unit': u.mag,
                                          'unit_str': 'mag',
                                          'doc_comment': 'Uncertainty of `f814w_vega_tot\'',
                                          'tab_comment': None
                                          },
            # flux in mJy
            'PHANGS_F275W_mJy_TOT': {'col_name': 'f275w_mJy_tot',
                                     'unit': u.mJy,
                                     'unit_str': 'mJy',
                                     'doc_comment': 'WFC3 f275w (NUV-band) total flux in mJy, MW foreground reddening '
                                                    'and aperture corrected. For more details on the aperture '
                                                    'correction see Deger et al. (2022) (2022MNRAS.510...32D). '
                                                    'Set to -9999 if source is not covered by HST filter. See also the '
                                                    '`no_coverage_flag\' column.',
                                     'tab_comment': 'WFC3 f275w (NUV-band) total flux in mJy, MW foreground '
                                                    'reddening and aperture corrected. Set to -9999 if source is not '
                                                    'covered by HST filter.'
                                     },
            'PHANGS_F275W_mJy_TOT_ERR': {'col_name': 'f275w_mJy_tot_err',
                                         'unit': u.mJy,
                                         'unit_str': 'mJy',
                                         'doc_comment': 'Uncertainty of `f275w_mJy_tot\'',
                                         'tab_comment': 'Uncertainty of `f275w_mJy_tot\''
                                         },
            'PHANGS_F336W_mJy_TOT': {'col_name': 'f336w_mJy_tot',
                                     'unit': u.mJy,
                                     'unit_str': 'mJy',
                                     'doc_comment': 'WFC3 f336w (U-band) total flux in mJy, MW foreground reddening '
                                                    'and aperture corrected. For more details on the aperture '
                                                    'correction see Deger et al. (2022) (2022MNRAS.510...32D). '
                                                    'Set to -9999 if source is not covered by HST filter. See also the '
                                                    '`no_coverage_flag\' column.',
                                     'tab_comment': None
                                     },
            'PHANGS_F336W_mJy_TOT_ERR': {'col_name': 'f336w_mJy_tot_err',
                                         'unit': u.mJy,
                                         'unit_str': 'mJy',
                                         'doc_comment': 'Uncertainty of `f336w_mJy_tot\'',
                                         'tab_comment': None
                                         },
            'PHANGS_F435W_mJy_TOT': {'col_name': 'f435w_mJy_tot',
                                     'unit': u.mJy,
                                     'unit_str': 'mJy',
                                     'doc_comment': 'Note that for Targets observed with the UVIS detector, '
                                                    'the filter name of the B-band is f438w. '
                                                    'WFC3 f435w (B-band) total flux in mJy, MW foreground reddening '
                                                    'and aperture corrected. For more details on the aperture '
                                                    'correction see Deger et al. (2022) (2022MNRAS.510...32D). '
                                                    'Set to -9999 if source is not covered by HST filter. See also the '
                                                    '`no_coverage_flag\' column.',
                                     'tab_comment': None
                                     },
            'PHANGS_F435W_mJy_TOT_ERR': {'col_name': 'f435w_mJy_tot_err',
                                         'unit': u.mJy,
                                         'unit_str': 'mJy',
                                         'doc_comment': 'Uncertainty of `f435w_mJy_tot\'',
                                         'tab_comment': None
                                         },
            'PHANGS_F438W_mJy_TOT': {'col_name': 'f438w_mJy_tot',
                                     'unit': u.mJy,
                                     'unit_str': 'mJy',
                                     'doc_comment': 'Note that for Targets observed with the WFC detector, '
                                                    'the filter name of the B-band is f435w. '
                                                    'WFC3 f438w (B-band) total flux in mJy, MW foreground reddening '
                                                    'and aperture corrected. For more details on the aperture '
                                                    'correction see Deger et al. (2022) (2022MNRAS.510...32D). '
                                                    'Set to -9999 if source is not covered by HST filter. See also the '
                                                    '`no_coverage_flag\' column.',
                                     'tab_comment': None
                                     },
            'PHANGS_F438W_mJy_TOT_ERR': {'col_name': 'f438w_mJy_tot_err',
                                         'unit': u.mJy,
                                         'unit_str': 'mJy',
                                         'doc_comment': 'Uncertainty of `f438w_mJy_tot\'',
                                         'tab_comment': None
                                         },
            'PHANGS_F555W_mJy_TOT': {'col_name': 'f555w_mJy_tot',
                                     'unit': u.mJy,
                                     'unit_str': 'mJy',
                                     'doc_comment': 'WFC3 f555w (V-band) total flux in mJy, MW foreground reddening '
                                                    'and aperture corrected. For more details on the aperture '
                                                    'correction see Deger et al. (2022) (2022MNRAS.510...32D). '
                                                    'Set to -9999 if source is not covered by HST filter. See also the '
                                                    '`no_coverage_flag\' column.',
                                     'tab_comment': None
                                     },
            'PHANGS_F555W_mJy_TOT_ERR': {'col_name': 'f555w_mJy_tot_err',
                                         'unit': u.mJy,
                                         'unit_str': 'mJy',
                                         'doc_comment': 'Uncertainty of `f555w_mJy_tot\'',
                                         'tab_comment': None
                                         },
            'PHANGS_F814W_mJy_TOT': {'col_name': 'f814w_mJy_tot',
                                     'unit': u.mJy,
                                     'unit_str': 'mJy',
                                     'doc_comment': 'WFC3 f814w (I-band) total flux in mJy, MW foreground reddening '
                                                    'and aperture corrected. For more details on the aperture '
                                                    'correction see Deger et al. (2022) (2022MNRAS.510...32D). '
                                                    'Set to -9999 if source is not covered by HST filter. See also the '
                                                    '`no_coverage_flag\' column.',
                                     'tab_comment': None
                                     },
            'PHANGS_F814W_mJy_TOT_ERR': {'col_name': 'f814w_mJy_tot_err',
                                         'unit': u.mJy,
                                         'unit_str': 'mJy',
                                         'doc_comment': 'Uncertainty of `f814w_mJy_tot\'',
                                         'tab_comment': None
                                         },
            # detection and coverage flags
            'PHANGS_NON_DETECTION_FLAG': {'col_name': 'non_detection_flag',
                                          'unit': None,
                                          'unit_str': 'int',
                                          'doc_comment': 'An integer denoting the number of bands in which the '
                                                         'photometry for the object was below the requested '
                                                         'signal-to-noise ratio (S/N=1). A value of 0 '
                                                         'indicates all five bands had detections. A value of 1 '
                                                         'indicates the object was detected in four bands but one '
                                                         'other band had only an upper limit, and a value of 2 '
                                                         'indicates the object was detected in three bands but two '
                                                         'other bands had only upper limits. By design, this flag '
                                                         'cannot be higher than 2.',
                                          'tab_comment': 'Integer denoting the number of bands in which the '
                                                         'photometry for the object was below the requested '
                                                         'signal-to-noise ratio (S/N=1). 0 indicates all five bands '
                                                         'had detections. A value of 1 and 2 means the object was '
                                                         'detected in four and three bands, resectively. By design, '
                                                         'this flag cannot be higher than 2.'
                                          },
            'PHANGS_NO_COVERAGE_FLAG': {'col_name': 'no_coverage_flag',
                                        'unit': None,
                                        'unit_str': 'int',
                                        'doc_comment': 'An integer denoting the number of bands for which the object '
                                                       'was outside the footprint of the sky coverage. The specific '
                                                       'bands for which there was no available image data can be '
                                                       'identified as photometry columns set to -9999.',
                                        'tab_comment': 'Integer denoting the number of bands with no coverage for '
                                                       'object. The specific bands can be identified as photometry '
                                                       'columns are set to -9999.'
                                        },
            # Concentration index
            'PHANGS_CI': {'col_name': 'ci',
                          'unit': None,
                          'unit_str': 'float',
                          'doc_comment': 'Concentration index: difference in magnitudes measured in 1 pix and 3 pix '
                                         'radii apertures (r=0.03962", r=3.*0.03962"=0.11886").',
                          'tab_comment': 'Concentration index: difference in magnitudes measured in 1 pix and 3 pix '
                                         'radii apertures.',
                          },

            # Fit results
            'PHANGS_AGE_MINCHISQ': {'col_name': 'age',
                                    'unit': u.Myr,
                                    'unit_str': 'Myr',
                                    'doc_comment': 'Cluster age corresponding to SED fit with minimum reduced chisq',
                                    'tab_comment': 'Cluster age corresponding to SED fit with minimum reduced chisq',
                                    },
            'PHANGS_AGE_MINCHISQ_ERR': {'col_name': 'age_err',
                                        'unit': u.Myr,
                                        'unit_str': 'Myr',
                                        'doc_comment': 'Uncertainty of `age\'',
                                        'tab_comment': 'Uncertainty of `age\'',
                                        },
            'PHANGS_MASS_MINCHISQ': {'col_name': 'mass',
                                     'unit': u.Myr,
                                     'unit_str': 'Myr',
                                     'doc_comment': 'Cluster mass corresponding to SED fit with minimum reduced chisq',
                                     'tab_comment': 'Cluster mass corresponding to SED fit with minimum reduced chisq',
                                     },
            'PHANGS_MASS_MINCHISQ_ERR': {'col_name': 'mass_err',
                                         'unit': u.Myr,
                                         'unit_str': 'Myr',
                                         'doc_comment': 'Uncertainty of `mass\'',
                                         'tab_comment': 'Uncertainty of `mass\'',
                                         },
            'PHANGS_EBV_MINCHISQ': {'col_name': 'ebv',
                                    'unit': u.M_sun,
                                    'unit_str': r'M$_{\odot}$',
                                    'doc_comment': 'Cluster reddening E(B-V) corresponding to SED fit with minimum '
                                                   'reduced chisq',
                                    'tab_comment': 'Cluster reddening E(B-V) corresponding to SED fit with minimum '
                                                   'reduced chisq',
                                    },
            'PHANGS_EBV_MINCHISQ_ERR': {'col_name': 'ebv_err',
                                        'unit': u.mag,
                                        'unit_str': 'mag',
                                        'doc_comment': 'Uncertainty of `ebv\'',
                                        'tab_comment': 'Uncertainty of `ebv\'',
                                        },
            'PHANGS_REDUCED_MINCHISQ': {'col_name': 'chisq',
                                        'unit': None,
                                        'unit_str': 'float',
                                        'doc_comment': 'The reduced chisq value from the SED fit that computed the '
                                                       'age, mass, and reddening of the cluster',
                                        'tab_comment': 'The reduced chisq value from the SED fit',
                                        },
            'PHANGS_SEDFIX_CATEGORY_DR4': {'col_name': 'sed_fix_category',
                                           'unit': None,
                                           'unit_str': 'str',
                                           'doc_comment': 'Category used to assign parameter grid for SED fit based on '
                                                          'color-color topology and H-alpha surface brightness. '
                                                          'Possible values can be `YRO\', `UNCHANGED\', or `OGC\'',
                                           'tab_comment': 'Category used to assign parameter grid for SED fit. '
                                                          'Possible values can be `YRO\', `UNCHANGED\', or `OGC\''
                                           },

        }

        self.identifier_columns = ['INDEX', 'ID_PHANGS_CLUSTERS_v1p2',
                                   # 'ID_PHANGS_ALLSOURCES',
                                   'PHANGS_X', 'PHANGS_Y', 'PHANGS_RA', 'PHANGS_DEC']

        self.classification_columns = ['PHANGS_CLUSTER_CLASS_HUMAN', 'PHANGS_CLUSTER_CLASS_ML_VGG',
                                       'PHANGS_CLUSTER_CLASS_ML_VGG_QUAL']

        self.example_photometry_columns = ['PHANGS_F275W_VEGA_TOT', 'PHANGS_F275W_VEGA_TOT_ERR',
                                           'PHANGS_F275W_mJy_TOT', 'PHANGS_F275W_mJy_TOT_ERR']

        self.detect_shape_columns = ['PHANGS_NON_DETECTION_FLAG', 'PHANGS_NO_COVERAGE_FLAG', 'PHANGS_CI']

        # list of columns entering table 1 (Observational properties)
        self.tab1_columns = (self.identifier_columns + self.classification_columns + self.example_photometry_columns +
                             self.detect_shape_columns)

        self.tab2_columns = ['INDEX', 'ID_PHANGS_CLUSTERS', 'ID_PHANGS_ALLSOURCES',
                             'PHANGS_AGE_MINCHISQ', 'PHANGS_AGE_MINCHISQ_ERR',
                             'PHANGS_MASS_MINCHISQ', 'PHANGS_MASS_MINCHISQ_ERR',
                             'PHANGS_EBV_MINCHISQ', 'PHANGS_EBV_MINCHISQ_ERR',
                             'PHANGS_REDUCED_MINCHISQ', 'PHANGS_SEDFIX_CATEGORY_DR4']

    def get_obs_table_column_list(self, target):
        """
        returns list with all column names for the observational table

        Parameters
        ----------
        target : str
            Target name

        Returns
        -------
        list : list
        """

        # get identifier and classification coulmns
        columns_list = self.identifier_columns + self.classification_columns

        # to add the photometry bands we need to access the individual bands for each target and put them in order
        band_list = []
        for band in list(set(self.hst_acs_wfc1_bands + self.hst_wfc3_uvis2_bands)):
            if band in (self.phangs_hst_obs_band_dict[target]['acs_wfc1_observed_bands'] +
                        self.phangs_hst_obs_band_dict[target]['wfc3_uvis_observed_bands']):
                band_list.append(band)
        band_list = self.sort_band_name_list(band_list=band_list)

        # create the photometry coulmn names
        vega_column_names = []
        mjy_column_names = []
        for band in band_list:
            vega_column_names.append('PHANGS_%s_VEGA_TOT' % band)
            vega_column_names.append('PHANGS_%s_VEGA_TOT_ERR' % band)
            mjy_column_names.append('PHANGS_%s_mJy_TOT' % band)
            mjy_column_names.append('PHANGS_%s_mJy_TOT_ERR' % band)
        columns_list += vega_column_names
        columns_list += mjy_column_names
        columns_list += self.detect_shape_columns

        return columns_list

    @staticmethod
    def sort_band_name_list(band_list):
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
            wave_list.append(int(band[1:-1]))

        # sort wavelength bands
        sort = np.argsort(wave_list)
        return list(np.array(band_list)[sort])

    @staticmethod
    def split_line(line_string, max_length):
        """
        Function to split a string into 2 up to the last spacing before string is larger than the maximal length
        Parameters
        ----------
        line_string : str
            String of the line which needs to be split
        max_length : int
            length of characters after which line needs to be split

        Returns
        -------
        shortened_str : str
            String of original line up to maximal character number
        rest_str : str
            Rest of the string beginning at the split of maximal character number
        """

        # check if line string is longer than maximal character length
        if len(line_string) < max_length:
            raise KeyError('line_string is already short enough.')

        # get all the indices of spacings in the string
        positions_of_spaces = [pos for pos, char in enumerate(line_string) if char == ' ']
        # get all the space position which are smaller than the maximal character length
        smaller_space_positions = np.where(np.array(positions_of_spaces) < max_length)
        # split the string
        shortened_str = line_string[:positions_of_spaces[smaller_space_positions[0][-1]]]
        rest_str = line_string[positions_of_spaces[smaller_space_positions[0][-1]]:]
        return shortened_str, rest_str

    def print_content_table(self, table_name='tab1', max_length_description=90):
        """
        Function to print table name
        ----------
        table_name : str
            String to provide the table_name can be tab1 for the observation paper or tab2 for the SED paper
        max_length_description : int
            length of characters in the description line. Some descriptions are too long and need to be spit into two

        Returns
        -------
        None
        """

        # print table header
        print("")
        print('\\hline\\hline')
        print('\\multicolumn{1}{c}{Column name} & '
              '\\multicolumn{1}{c}{Unit} & '
              '\\multicolumn{1}{c}{Description} \\\\ ')
        print('\\hline')

        # print table content
        # get the according set of column names
        table_keys = getattr(self, table_name + '_columns')
        # for loop to go through all column names
        for key in table_keys:
            name_str = self.cat_info[key]['col_name'].replace('_', '\\_')
            unit_str = self.cat_info[key]['unit_str']
            content_str = self.cat_info[key]['tab_comment'].replace('_', '\\_')
            # check if the description part is not too long
            if len(content_str) < max_length_description:
                # in case of an adequate description length we can simply print the table line
                line_str = (name_str + ' & ' + unit_str + ' & ' + content_str + ' \\\\ ')
                print(line_str)
            else:
                # in order to shorten the description line we need to split the description string
                shortened_str, current_str = self.split_line(line_string=content_str, max_length=max_length_description)
                # print first line of the string
                line_str_1 = (name_str + ' & ' + unit_str + ' & ' + shortened_str + ' \\\\ ')
                print(line_str_1)
                # use while loop to keep on dividing the rest of the description string till it is short enough
                while_flag = True
                while while_flag:
                    if len(current_str) > max_length_description:
                        shortened_str, current_str = self.split_line(line_string=current_str,
                                                                     max_length=max_length_description)
                    else:
                        shortened_str = current_str
                        # after this condition is fulfilled the while loop will stop after the current iteration
                        while_flag = False
                    # print out the next line
                    line_str_2 = (' &  &  \\,\\,\\,' + shortened_str + ' \\\\ ')
                    print(line_str_2)

        print('\\hline')
