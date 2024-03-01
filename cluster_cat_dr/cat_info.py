"""
Script to gather all information going into data-description tables, documentation & released catalogs
"""
import astropy.units as u
import numpy as np
from phangs_info import PhangsObsInfo, PhysParams


class CatalogInfo(PhangsObsInfo, PhysParams):
    """
    Class to gather all catalog information as attributes
    """
    def __init__(self):
        super().__init__()

        # Catalog info is a dictionary providing the final column name for the fits files of the final data release.
        # The keys are sorted for the column names used in the final internal release.
        # 'col_name' provides the name we will use in the articles, documentations and DR files.
        # 'units' provides the units in astropy.units format for example seconds would be units.s
        # 'units_str' provides the unit in string format which we will provide in the documentation.
        # 'doc_comment' is the description for the documentation.
        # 'tab_comment' is the description going into the table in the article.
        self.cat_info = {

            # cluster identifier
            'INDEX': {'col_name': 'INDEX',
                      # 'col_name': 'index',
                      'unit': None,
                      'unit_str': 'int',
                      'doc_comment': 'A running index from 1 to N, where N is the total number of objects in the '
                                     'catalog. Objects are sorted in order of increasing Y pixel value on image.',
                      'tab_comment': 'Running index from 1 to N for each individual target'},

            'id_phangs_cluster': {  # 'col_name': 'id_phangs_cluster',
                                  'col_name': 'ID_PHANGS_CLUSTER',
                                  'unit': None,
                                  'unit_str': 'int',
                                  'doc_comment': 'PHANGS cluster ID for each individual object classified as '
                                                 'class 1,2 or 3, ordered by '
                                                 'increasing Y pixel coordinate',
                                  'tab_comment': 'PHANGS cluster ID for each individual object classified as '
                                                 'class 1,2 or 3, ordered by '
                                                 'increasing Y pixel coordinate'},


            'ID_PHANGS_CLUSTERS_v1p2': {# 'col_name': 'ID_PHANGS_CLUSTERS_v1p2',
                                        # 'col_name': 'id_phangs_candidate',
                                        'col_name': 'ID_PHANGS_CANDIDATE',
                                        'unit': None,
                                        'unit_str': 'int',
                                        'doc_comment': 'Running PHANGS candidate ID. These values correspond to the '
                                                       'candidate catalogs produced by the PHANGS-HST team. Objects '
                                                       'are sorted in order of increasing Y pixel value on image.',
                                        'tab_comment': 'ID in the PHANGS-HST candidate catalog for each '
                                                       'individual target, for cross-identification.'
                                        },
            'ID_PHANGS_ALLSOURCES_v1p2': {# 'col_name': 'ID_PHANGS_ALLSOURCES',
                                     # 'col_name': 'id_phangs_allsource',
                                     'col_name': 'ID_PHANGS_ALLSOURCES',
                                     'unit': None,
                                     'unit_str': 'int',
                                     'doc_comment': 'Running PHANGS source ID. These values correspond to the '
                                                    'initial detection catalogs produced by the PHANGS-HST team. '
                                                    'Objects are sorted in order of increasing Y pixel value on image.',
                                     'tab_comment': 'ID in the initial PHANGS-HST “all-source” detection catalog '
                                                    'for each individual target, for cross-identification.'
                                     },

            # coordinates of each cluster
            'PHANGS_X': {'col_name': 'PHANGS_X',
                         # 'col_name': 'x_pix_coord',
                         'unit': None,
                         'unit_str': 'pix',
                         'doc_comment': 'X pixel coordinate on HST image. Scale = 0.03962 arcsec/pixel. '
                                        'Numeration goes from 0 to n-1, where n is the number of Pixels on the X-axis. '
                                        'The value of 1 must be added to some other systems which star counting at 1 '
                                        '(e.g., IRAF, ds9).',
                         'tab_comment': 'X coordinates on HST X-pixel grid (0...n-1). Scale = 0.03962 arcsec/pixel.'
                         },
            'PHANGS_Y': {'col_name': 'PHANGS_Y',
                         # 'col_name': 'y_pix_coord',
                         'unit': None,
                         'unit_str': 'pix',
                         'doc_comment': 'Y pixel coordinate on HST image. Scale = 0.03962 arcsec/pixel. '
                                        'Numeration goes from 0 to n-1, where n is the number of Pixels on the Y-axis. '
                                        'The value of 1 must be added to some other systems which star counting at 1 '
                                        '(e.g., IRAF, ds9).',
                         'tab_comment': 'Y coordinates on HST Y-pixel grid (0...n-1). Scale = 0.03962 arcsec/pixel.'
                         },

            'PHANGS_RA': {'col_name': 'PHANGS_RA',
                          # 'col_name': 'ra',
                          'unit': u.deg,
                          'unit_str': 'deg',
                          'doc_comment': 'J2000 Right ascension, ICRS frame, calibrated against selected Gaia sources.',
                          'tab_comment': 'J2000 Right ascension, ICRS frame, calibrated against selected Gaia sources.'
                          },
            'PHANGS_DEC': {'col_name': 'PHANGS_DEC',
                           # 'col_name': 'dec',
                           'unit': u.deg,
                           'unit_str': 'deg',
                           'doc_comment': 'J2000 Declination, ICRS frame, calibrated against selected Gaia sources.',
                           'tab_comment': 'J2000 Declination, ICRS frame, calibrated against selected Gaia sources.'
                           },

            # Morphological Classification
            'PHANGS_CLUSTER_CLASS_HUMAN': {'col_name': 'PHANGS_CLUSTER_CLASS_HUMAN',
                                           # 'col_name': 'class_human',
                                           'unit': None,
                                           'unit_str': 'int',
                                           'doc_comment': 'Human visual classification determined by coauthor '
                                                          'Brad Whitmore (BCW). Further details in '
                                                          'Whitmore et al. (2021) (2021MNRAS.506.5294W).'
                                                          'The classification is encoded in integers: '
                                                          '1 and 2 for class 1 and 2 compact clusters, respectively and '
                                                          '3 for class 3 compact associations. '
                                                          'Intengers > 3 are artefacts, see above.'
                                                          'Classification numbers for artefacts are described above',
                                           'tab_comment': 'Cluster class assigned through visual inspection. '
                                                          'Integers 1 and 2 represent C1 and C2 compact clusters. '
                                                          '3 stands for C3 compact associations. '
                                                          'Intengers > 3 are artefacts.'
                                           },
            'PHANGS_CLUSTER_CLASS_ML_VGG': {'col_name': 'PHANGS_CLUSTER_CLASS_ML_VGG',
                                            # 'col_name': 'class_ml_vgg',
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
                                                           'class 1, 2 and class 3 compact associations, respectively. '
                                                           'Intengers > 3 are artefacts, see above.'
                                                           ' Classification numbers for artefacts are described above',
                                            'tab_comment': 'Cluster class determined by VGG neural network. '
                                                           'Integers 1 and 2 represent C1 and C2 compact clusters. '
                                                           '3 stands for C3 compact associations. '
                                                           'Intengers > 3 are artefacts.'
                                            },
            'PHANGS_CLUSTER_CLASS_ML_VGG_QUAL': {'col_name': 'PHANGS_CLUSTER_CLASS_ML_VGG_QUAL',
                                                 # 'col_name': 'class_ml_vgg_qual',
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
            'class_ml_vgg_corr': {# 'col_name': 'class_ml_vgg_corr',
                                  'col_name': 'PHANGS_CLUSTER_CLASS_ML_VGG_CORR',
                                  'unit': None,
                                  'unit_str': 'int',
                                  'doc_comment': 'Same classification as in column `cluster_class_ml\' '
                                                 'but corrected for human identified artefacts. '
                                                 'Thus, all object classes identified as artefacts in column '
                                                 '`class_human\' are replaced.',
                                  'tab_comment': ' '
                                  },
            # Photometry
            # placeholder example
            'PHANGS_[BAND]_VEGA_TOT': {'col_name': 'PHANGS_[BAND]_VEGA',
                                      'unit': u.mag,
                                      'unit_str': 'mag',
                                      'doc_comment': '',
                                      'tab_comment': 'HST band apparent vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. Set to -9999 if source is not '
                                                     'covered by HST filter.'
                                      },
            'PHANGS_[BAND]_VEGA_TOT_ERR': {'col_name': 'PHANGS_[BAND]_VEGA_ERR',
                                          # 'col_name': 'f275w_vega_err',
                                          'unit': u.mag,
                                          'unit_str': 'mag',
                                          'doc_comment': 'Uncertainty of `[BAND]_VEGA\'',
                                          'tab_comment': 'Uncertainty of `[BAND]_VEGA\''
                                          },
            'PHANGS_[BAND]_mJy_TOT': {'col_name': 'PHANGS_[BAND]_mJy',
                                     # 'col_name': 'f275w_mJy',
                                     'unit': u.mJy,
                                     'unit_str': 'mJy',
                                     'doc_comment': '',
                                     'tab_comment': 'HST band flux in mJy, MW foreground '
                                                    'reddening and aperture corrected. Set to -9999 if source is not '
                                                    'covered by HST filter.'
                                     },
            'PHANGS_[BAND]_mJy_TOT_ERR': {'col_name': 'PHANGS_[BAND]_mJy_ERR',
                                         # 'col_name': 'f275w_mJy_err',
                                         'unit': u.mJy,
                                         'unit_str': 'mJy',
                                         'doc_comment': 'Uncertainty of `[BAND]_mJy\'',
                                         'tab_comment': 'Uncertainty of `[BAND]_mJy\''
                                         },

            # Vega magnitudes
            'PHANGS_F275W_VEGA_TOT': {'col_name': 'PHANGS_F275W_VEGA',
                                      # 'col_name': 'f275w_vega',
                                      'unit': u.mag,
                                      'unit_str': 'mag',
                                      'doc_comment': 'WFC3 f275w (NUV-band) apparent apparent vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. For more details on the '
                                                     'aperture correction see '
                                                     'Deger et al. (2022) (2022MNRAS.510...32D). Set to -9999 if '
                                                     'source is not covered by HST filter. See also the '
                                                     '`no_coverage_flag\' column.',
                                      'tab_comment': 'WFC3 f275w (NUV-band) apparent vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. Set to -9999 if source is not '
                                                     'covered by HST filter.'
                                      },
            'PHANGS_F275W_VEGA_TOT_ERR': {'col_name': 'PHANGS_F275W_VEGA_ERR',
                                          # 'col_name': 'f275w_vega_err',
                                          'unit': u.mag,
                                          'unit_str': 'mag',
                                          'doc_comment': 'Uncertainty of `F275W_VEGA\'',
                                          'tab_comment': 'Uncertainty of `F275W_VEGA\''
                                          },
            'PHANGS_F336W_VEGA_TOT': {'col_name': 'PHANGS_F336W_VEGA',
                                      # 'col_name': 'f336w_vega',
                                      'unit': u.mag,
                                      'unit_str': 'mag',
                                      'doc_comment': 'WFC3 f336w (U-band) apparent vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. For more details on the '
                                                     'aperture correction see '
                                                     'Deger et al. (2022) (2022MNRAS.510...32D). Set to -9999 if '
                                                     'source is not covered by HST filter. See also the '
                                                     '`no_coverage_flag\' column.',
                                      'tab_comment': None
                                      },
            'PHANGS_F336W_VEGA_TOT_ERR': {'col_name': 'PHANGS_F336W_VEGA_ERR',
                                          # 'col_name': 'f336w_vega_err',
                                          'unit': u.mag,
                                          'unit_str': 'mag',
                                          'doc_comment': 'Uncertainty of `F336W_VEGA\'',
                                          'tab_comment': None
                                          },
            'PHANGS_F435W_VEGA_TOT': {'col_name': 'PHANGS_F435W_VEGA',
                                      # 'col_name': 'f435w_vega',
                                      'unit': u.mag,
                                      'unit_str': 'mag',
                                      'doc_comment': 'ACS WFC3 f435w (B-band) apparent vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. '
                                                     'This column is only provided for targets observed with ACS WFC3 '
                                                     'detector. '
                                                     'For more details on the aperture correction see '
                                                     'Deger et al. (2022) (2022MNRAS.510...32D). Set to -9999 if '
                                                     'source is not covered by HST filter. See also the '
                                                     '`no_coverage_flag\' column.',
                                      'tab_comment': None
                                      },
            'PHANGS_F435W_VEGA_TOT_ERR': {'col_name': 'PHANGS_F435W_VEGA_ERR',
                                          # 'col_name': 'f435w_vega_err',
                                          'unit': u.mag,
                                          'unit_str': 'mag',
                                          'doc_comment': 'Uncertainty of `F435W_VEGA\'',
                                          'tab_comment': None
                                          },
            'PHANGS_F438W_VEGA_TOT': {'col_name': 'PHANGS_F438W_VEGA',
                                      # 'col_name': 'f438w_vega',
                                      'unit': u.mag,
                                      'unit_str': 'mag',
                                      'doc_comment': 'WFC3 f438w (B-band) apparent vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. '
                                                     'This column is only provided for targets observed with UVIS WFC3 '
                                                     'detector. '
                                                     'For more details on the '
                                                     'aperture correction see '
                                                     'Deger et al. (2022) (2022MNRAS.510...32D). Set to -9999 if '
                                                     'source is not covered by HST filter. See also the '
                                                     '`no_coverage_flag\' column.',
                                      'tab_comment': None
                                      },
            'PHANGS_F438W_VEGA_TOT_ERR': {'col_name': 'PHANGS_F438W_VEGA_ERR',
                                          # 'col_name': 'f438w_vega_err',
                                          'unit': u.mag,
                                          'unit_str': 'mag',
                                          'doc_comment': 'Uncertainty of `F438W_VEGA\'',
                                          'tab_comment': None
                                          },
            'PHANGS_F555W_VEGA_TOT': {'col_name': 'PHANGS_F555W_VEGA',
                                      # 'col_name': 'f555w_vega',
                                      'unit': u.mag,
                                      'unit_str': 'mag',
                                      'doc_comment': 'WFC3 f555w (V-band) apparent vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. For more details on the '
                                                     'aperture correction see '
                                                     'Deger et al. (2022) (2022MNRAS.510...32D). Set to -9999 if '
                                                     'source is not covered by HST filter. See also the '
                                                     '`no_coverage_flag\' column.',
                                      'tab_comment': None
                                      },
            'PHANGS_F555W_VEGA_TOT_ERR': {'col_name': 'PHANGS_F555W_VEGA_ERR',
                                          # 'col_name': 'f555w_vega_err',
                                          'unit': u.mag,
                                          'unit_str': 'mag',
                                          'doc_comment': 'Uncertainty of `F555W_VEGA\'',
                                          'tab_comment': None
                                          },
            'PHANGS_F814W_VEGA_TOT': {'col_name': 'PHANGS_F814W_VEGA',
                                      # 'col_name': 'f814w_vega',
                                      'unit': u.mag,
                                      'unit_str': 'mag',
                                      'doc_comment': 'WFC3 f814w (I-band) apparent vega magnitude, MW foreground '
                                                     'reddening and aperture corrected. For more details on the '
                                                     'aperture correction see '
                                                     'Deger et al. (2022) (2022MNRAS.510...32D). Set to -9999 if '
                                                     'source is not covered by HST filter. See also the '
                                                     '`no_coverage_flag\' column.',
                                      'tab_comment': None
                                      },
            'PHANGS_F814W_VEGA_TOT_ERR': {'col_name': 'PHANGS_F814W_VEGA_ERR',
                                          # 'col_name': 'f814w_vega_err',
                                          'unit': u.mag,
                                          'unit_str': 'mag',
                                          'doc_comment': 'Uncertainty of `F814W_VEGA\'',
                                          'tab_comment': None
                                          },
            # flux in mJy
            'PHANGS_F275W_mJy_TOT': {'col_name': 'PHANGS_F275W_mJy',
                                     # 'col_name': 'f275w_mJy',
                                     'unit': u.mJy,
                                     'unit_str': 'mJy',
                                     'doc_comment': 'WFC3 f275w (NUV-band) flux in mJy, MW foreground reddening '
                                                    'and aperture corrected. For more details on the aperture '
                                                    'correction see Deger et al. (2022) (2022MNRAS.510...32D). '
                                                    'Set to -9999 if source is not covered by HST filter. See also the '
                                                    '`no_coverage_flag\' column.',
                                     'tab_comment': 'WFC3 f275w (NUV-band) flux in mJy, MW foreground '
                                                    'reddening and aperture corrected. Set to -9999 if source is not '
                                                    'covered by HST filter.'
                                     },
            'PHANGS_F275W_mJy_TOT_ERR': {'col_name': 'PHANGS_F275W_mJy_ERR',
                                         # 'col_name': 'f275w_mJy_err',
                                         'unit': u.mJy,
                                         'unit_str': 'mJy',
                                         'doc_comment': 'Uncertainty of `F275W_mJy\'',
                                         'tab_comment': 'Uncertainty of `F275W_mJy\''
                                         },
            'PHANGS_F336W_mJy_TOT': {'col_name': 'PHANGS_F336W_mJy',
                                     # 'col_name': 'f336w_mJy',
                                     'unit': u.mJy,
                                     'unit_str': 'mJy',
                                     'doc_comment': 'WFC3 f336w (U-band) flux in mJy, MW foreground reddening '
                                                    'and aperture corrected. For more details on the aperture '
                                                    'correction see Deger et al. (2022) (2022MNRAS.510...32D). '
                                                    'Set to -9999 if source is not covered by HST filter. See also the '
                                                    '`no_coverage_flag\' column.',
                                     'tab_comment': None
                                     },
            'PHANGS_F336W_mJy_TOT_ERR': {'col_name': 'PHANGS_F336W_mJy_ERR',
                                         # 'col_name': 'f336w_mJy_err',
                                         'unit': u.mJy,
                                         'unit_str': 'mJy',
                                         'doc_comment': 'Uncertainty of `F336W_mJy\'',
                                         'tab_comment': None
                                         },
            'PHANGS_F435W_mJy_TOT': {'col_name': 'PHANGS_F435W_mJy',
                                     # 'col_name': 'f435w_mJy',
                                     'unit': u.mJy,
                                     'unit_str': 'mJy',
                                     'doc_comment': 'Note that for Targets observed with the UVIS detector, '
                                                    'the filter name of the B-band is f438w. '
                                                    'WFC3 f435w (B-band) flux in mJy, MW foreground reddening '
                                                    'and aperture corrected. For more details on the aperture '
                                                    'correction see Deger et al. (2022) (2022MNRAS.510...32D). '
                                                    'Set to -9999 if source is not covered by HST filter. See also the '
                                                    '`no_coverage_flag\' column.',
                                     'tab_comment': None
                                     },
            'PHANGS_F435W_mJy_TOT_ERR': {'col_name': 'PHANGS_F435W_mJy_ERR',
                                         # 'col_name': 'f435w_mJy_err',
                                         'unit': u.mJy,
                                         'unit_str': 'mJy',
                                         'doc_comment': 'Uncertainty of `F435W_mJy\'',
                                         'tab_comment': None
                                         },
            'PHANGS_F438W_mJy_TOT': {'col_name': 'PHANGS_F438W_mJy',
                                     # 'col_name': 'f438w_mJy',
                                     'unit': u.mJy,
                                     'unit_str': 'mJy',
                                     'doc_comment': 'Note that for Targets observed with the WFC detector, '
                                                    'the filter name of the B-band is f435w. '
                                                    'WFC3 f438w (B-band) flux in mJy, MW foreground reddening '
                                                    'and aperture corrected. For more details on the aperture '
                                                    'correction see Deger et al. (2022) (2022MNRAS.510...32D). '
                                                    'Set to -9999 if source is not covered by HST filter. See also the '
                                                    '`no_coverage_flag\' column.',
                                     'tab_comment': None
                                     },
            'PHANGS_F438W_mJy_TOT_ERR': {'col_name': 'PHANGS_F438W_mJy_ERR',
                                         # 'col_name': 'f438w_mJy_err',
                                         'unit': u.mJy,
                                         'unit_str': 'mJy',
                                         'doc_comment': 'Uncertainty of `F438W_mJy\'',
                                         'tab_comment': None
                                         },
            'PHANGS_F555W_mJy_TOT': {'col_name': 'PHANGS_F555W_mJy',
                                     # 'col_name': 'f555w_mJy',
                                     'unit': u.mJy,
                                     'unit_str': 'mJy',
                                     'doc_comment': 'WFC3 f555w (V-band) flux in mJy, MW foreground reddening '
                                                    'and aperture corrected. For more details on the aperture '
                                                    'correction see Deger et al. (2022) (2022MNRAS.510...32D). '
                                                    'Set to -9999 if source is not covered by HST filter. See also the '
                                                    '`no_coverage_flag\' column.',
                                     'tab_comment': None
                                     },
            'PHANGS_F555W_mJy_TOT_ERR': {'col_name': 'PHANGS_F555W_mJy_ERR',
                                         # 'col_name': 'f555w_mJy_err',
                                         'unit': u.mJy,
                                         'unit_str': 'mJy',
                                         'doc_comment': 'Uncertainty of `F555W_mJy\'',
                                         'tab_comment': None
                                         },
            'PHANGS_F814W_mJy_TOT': {'col_name': 'PHANGS_F814W_mJy',
                                     # 'col_name': 'f814w_mJy',
                                     'unit': u.mJy,
                                     'unit_str': 'mJy',
                                     'doc_comment': 'WFC3 f814w (I-band) flux in mJy, MW foreground reddening '
                                                    'and aperture corrected. For more details on the aperture '
                                                    'correction see Deger et al. (2022) (2022MNRAS.510...32D). '
                                                    'Set to -9999 if source is not covered by HST filter. See also the '
                                                    '`no_coverage_flag\' column.',
                                     'tab_comment': None
                                     },
            'PHANGS_F814W_mJy_TOT_ERR': {'col_name': 'PHANGS_F814W_mJy_ERR',
                                         # 'col_name': 'f814w_mJy_err',
                                         'unit': u.mJy,
                                         'unit_str': 'mJy',
                                         'doc_comment': 'Uncertainty of `F814W_mJy\'',
                                         'tab_comment': None
                                         },
            # detection and coverage flags
            'PHANGS_NON_DETECTION_FLAG': {'col_name': 'PHANGS_NON_DETECTION_FLAG',
                                          # 'col_name': 'non_detection_flag',
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
            'PHANGS_NO_COVERAGE_FLAG': {'col_name': 'PHANGS_NO_COVERAGE_FLAG',
                                        # 'col_name': 'no_coverage_flag',
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
            'PHANGS_CI': {'col_name': 'PHANGS_CI',
                          # 'col_name': 'ci',
                          'unit': u.mag,
                          'unit_str': 'mag',
                          'doc_comment': 'Concentration index: difference in magnitudes measured in 1 pix and 3 pix '
                                         'radii apertures (r=0.03962", r=3.*0.03962"=0.11886").',
                          'tab_comment': 'Concentration index: difference in magnitudes measured in 1 pix and 3 pix '
                                         'radii apertures.',
                          },
            # Concentration index
            'cc_class': {'col_name': 'PHANGS_COLORCOLOR_CLASS_UBVI',
                          'unit': None,
                          'unit_str': 'str',
                          'doc_comment': 'Flag to identify in which region on the color-color diagram the object was '
                                         'associated with. Values are \"ycl\" (young cluster locus), \"map\" '
                                         '(middle aged plume) \"ogcc\" (old globular cluster clump) or \"outside\" '
                                         '(outside the main regions and therefore not classified). '
                                         'A detailed description is found in Section 4.4 of Maschmann et al. (in prep)',
                          'tab_comment': 'Flag to identify in which region on the color-color diagram the object was '
                                         'associated with. Values are `ycl\' (young cluster locus), `map\' '
                                         '(middle aged plume) `ogcc\' (old globular cluster clump) or `outside\' '
                                         '(outside the main regions and therefore not classified). '
                                         'A detailed description is found in Section\\,\\ref{ssect:cc_regions}.',
                          },

            # Fit results
            'PHANGS_AGE_MINCHISQ': {'col_name': 'PHANGS_AGE_MINCHISQ',
                                    # 'col_name': 'age',
                                    'unit': u.Myr,
                                    'unit_str': 'Myr',
                                    'doc_comment': 'Cluster age corresponding to SED fit with minimum reduced chisq',
                                    'tab_comment': 'Cluster age corresponding to SED fit with minimum reduced chisq',
                                    },
            'PHANGS_AGE_MINCHISQ_ERR': {'col_name': 'PHANGS_AGE_MINCHISQ_ERR',
                                        # 'col_name': 'age_err',
                                        'unit': u.Myr,
                                        'unit_str': 'Myr',
                                        'doc_comment': 'Uncertainty of `age\'',
                                        'tab_comment': 'Uncertainty of `age\'',
                                        },
            'PHANGS_MASS_MINCHISQ': {'col_name': 'PHANGS_MASS_MINCHISQ',
                                     # 'col_name': 'mass',
                                     'unit': u.M_sun,
                                     'unit_str': 'M_sun',
                                     'doc_comment': 'Cluster mass corresponding to SED fit with minimum reduced chisq',
                                     'tab_comment': 'Cluster mass corresponding to SED fit with minimum reduced chisq',
                                     },
            'PHANGS_MASS_MINCHISQ_ERR': {'col_name': 'PHANGS_MASS_MINCHISQ_ERR',
                                         # 'col_name': 'mass_err',
                                         'unit': u.M_sun,
                                         'unit_str': 'M_sun',
                                         'doc_comment': 'Uncertainty of `mass\'',
                                         'tab_comment': 'Uncertainty of `mass\'',
                                         },
            'PHANGS_EBV_MINCHISQ': {'col_name': 'PHANGS_EBV_MINCHISQ',
                                    # 'col_name': 'ebv',
                                    'unit': u.mag,
                                    'unit_str': r'mag',
                                    'doc_comment': 'Cluster reddening E(B-V) corresponding to SED fit with minimum '
                                                   'reduced chisq',
                                    'tab_comment': 'Cluster reddening E(B-V) corresponding to SED fit with minimum '
                                                   'reduced chisq',
                                    },
            'PHANGS_EBV_MINCHISQ_ERR': {'col_name': 'PHANGS_EBV_MINCHISQ_ERR',
                                        # 'col_name': 'ebv_err',
                                        'unit': u.mag,
                                        'unit_str': 'mag',
                                        'doc_comment': 'Uncertainty of `ebv\'',
                                        'tab_comment': 'Uncertainty of `ebv\'',
                                        },
            'PHANGS_REDUCED_MINCHISQ': {'col_name': 'PHANGS_REDUCEDCHISQ_MINCHISQ',
                                        # 'col_name': 'chisq',
                                        'unit': None,
                                        'unit_str': 'float',
                                        'doc_comment': 'The reduced chisq value from the SED fit that computed the '
                                                       'age, mass, and reddening of the cluster',
                                        'tab_comment': 'The reduced chisq value from the SED fit',
                                        },


            'NBHa_intensity': {'col_name': 'NBHa_intensity', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'NBHa_intensity_medsub': {'col_name': 'NBHa_intensity_medsub', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'NBHa_mask_medsub_lev1': {'col_name': 'NBHa_mask_medsub_lev1', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'NBHa_mask_medsub_lev2': {'col_name': 'NBHa_mask_medsub_lev2', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'NBHa_mask_medsub_lev3': {'col_name': 'NBHa_mask_medsub_lev3', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'NBHa_mask_medsub_lev5': {'col_name': 'NBHa_mask_medsub_lev5', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'NBHa_mask_medsub_lev10': {'col_name': 'NBHa_mask_medsub_lev10', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'NBHa_intensity_medsub_inseg': {'col_name': 'NBHa_intensity_medsub_inseg', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'NBHa_HIIreg': {'col_name': 'NBHa_HIIreg', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'OGCpoly_UBVI': {'col_name': 'OGCpoly_UBVI', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'FlagOGC_UBVI': {'col_name': 'FlagOGC_UBVI', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'OGC_Cull1': {'col_name': 'OGC_Cull1', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'OGCpoly_BVI': {'col_name': 'OGCpoly_BVI', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'FlagOGC_BVI': {'col_name': 'FlagOGC_BVI', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'OGC_Cull2': {'col_name': 'OGC_Cull2', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'OGCrpoly_UBVI': {'col_name': 'OGCrpoly_UBVI', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'FlagOGCr_UBVI': {'col_name': 'FlagOGCr_UBVI', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'OGCr_Cull1': {'col_name': 'OGCr_Cull1', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'OGCrpoly_BVI': {'col_name': 'OGCrpoly_BVI', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'FlagOGCr_BVI': {'col_name': 'FlagOGCr_BVI', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'OGCr_Cull2': {'col_name': 'OGCr_Cull2', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'YROpoly_UBVI': {'col_name': 'YROpoly_UBVI', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'FlagYRO_UBVI': {'col_name': 'FlagYRO_UBVI', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'YRO_Cull1': {'col_name': 'YRO_Cull1', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'YROpoly_BVI': {'col_name': 'YROpoly_BVI', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'FlagYRO_BVI': {'col_name': 'FlagYRO_BVI', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'YRO_Cull2': {'col_name': 'YRO_Cull2', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_youngestmode_age': {'col_name': 'mm_Zyoung_youngestmode_age', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_youngestmode_ebv': {'col_name': 'mm_Zyoung_youngestmode_ebv', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_youngestmode_mass': {'col_name': 'mm_Zyoung_youngestmode_mass', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_youngestmode_likelihood': {'col_name': 'mm_Zyoung_youngestmode_likelihood', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_bestmode_age': {'col_name': 'mm_Zyoung_bestmode_age', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_bestmode_ebv': {'col_name': 'mm_Zyoung_bestmode_ebv', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_bestmode_mass': {'col_name': 'mm_Zyoung_bestmode_mass', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_bestmode_likelihood': {'col_name': 'mm_Zyoung_bestmode_likelihood', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_youngestmode_age_limlo': {'col_name': 'mm_Zyoung_youngestmode_age_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_youngestmode_ebv_limlo': {'col_name': 'mm_Zyoung_youngestmode_ebv_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_youngestmode_mass_limlo': {'col_name': 'mm_Zyoung_youngestmode_mass_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_bestmode_age_limlo': {'col_name': 'mm_Zyoung_bestmode_age_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_bestmode_ebv_limlo': {'col_name': 'mm_Zyoung_bestmode_ebv_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_bestmode_mass_limlo': {'col_name': 'mm_Zyoung_bestmode_mass_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_youngestmode_age_limhi': {'col_name': 'mm_Zyoung_youngestmode_age_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_youngestmode_ebv_limhi': {'col_name': 'mm_Zyoung_youngestmode_ebv_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_youngestmode_mass_limhi': {'col_name': 'mm_Zyoung_youngestmode_mass_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_bestmode_age_limhi': {'col_name': 'mm_Zyoung_bestmode_age_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_bestmode_ebv_limhi': {'col_name': 'mm_Zyoung_bestmode_ebv_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Zyoung_bestmode_mass_limhi': {'col_name': 'mm_Zyoung_bestmode_mass_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_youngestmode_age': {'col_name': 'mm_Z_youngestmode_age', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_youngestmode_ebv': {'col_name': 'mm_Z_youngestmode_ebv', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_youngestmode_mass': {'col_name': 'mm_Z_youngestmode_mass', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_youngestmode_likelihood': {'col_name': 'mm_Z_youngestmode_likelihood', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_bestmode_age': {'col_name': 'mm_Z_bestmode_age', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_bestmode_ebv': {'col_name': 'mm_Z_bestmode_ebv', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_bestmode_mass': {'col_name': 'mm_Z_bestmode_mass', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_bestmode_likelihood': {'col_name': 'mm_Z_bestmode_likelihood', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_youngestmode_age_limlo': {'col_name': 'mm_Z_youngestmode_age_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_youngestmode_ebv_limlo': {'col_name': 'mm_Z_youngestmode_ebv_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_youngestmode_mass_limlo': {'col_name': 'mm_Z_youngestmode_mass_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_bestmode_age_limlo': {'col_name': 'mm_Z_bestmode_age_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_bestmode_ebv_limlo': {'col_name': 'mm_Z_bestmode_ebv_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_bestmode_mass_limlo': {'col_name': 'mm_Z_bestmode_mass_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_youngestmode_age_limhi': {'col_name': 'mm_Z_youngestmode_age_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_youngestmode_ebv_limhi': {'col_name': 'mm_Z_youngestmode_ebv_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_youngestmode_mass_limhi': {'col_name': 'mm_Z_youngestmode_mass_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_bestmode_age_limhi': {'col_name': 'mm_Z_bestmode_age_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_bestmode_ebv_limhi': {'col_name': 'mm_Z_bestmode_ebv_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_Z_bestmode_mass_limhi': {'col_name': 'mm_Z_bestmode_mass_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_oldestmode_age': {'col_name': 'mm_lowZ_oldestmode_age', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_oldestmode_ebv': {'col_name': 'mm_lowZ_oldestmode_ebv', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_oldestmode_mass': {'col_name': 'mm_lowZ_oldestmode_mass', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_oldestmode_likelihood': {'col_name': 'mm_lowZ_oldestmode_likelihood', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_bestmode_age': {'col_name': 'mm_lowZ_bestmode_age', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_bestmode_ebv': {'col_name': 'mm_lowZ_bestmode_ebv', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_bestmode_mass': {'col_name': 'mm_lowZ_bestmode_mass', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_bestmode_likelihood': {'col_name': 'mm_lowZ_bestmode_likelihood', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_oldestmode_age_limlo': {'col_name': 'mm_lowZ_oldestmode_age_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_oldestmode_ebv_limlo': {'col_name': 'mm_lowZ_oldestmode_ebv_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_oldestmode_mass_limlo': {'col_name': 'mm_lowZ_oldestmode_mass_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_bestmode_age_limlo': {'col_name': 'mm_lowZ_bestmode_age_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_bestmode_ebv_limlo': {'col_name': 'mm_lowZ_bestmode_ebv_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_bestmode_mass_limlo': {'col_name': 'mm_lowZ_bestmode_mass_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_oldestmode_age_limhi': {'col_name': 'mm_lowZ_oldestmode_age_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_oldestmode_ebv_limhi': {'col_name': 'mm_lowZ_oldestmode_ebv_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_oldestmode_mass_limhi': {'col_name': 'mm_lowZ_oldestmode_mass_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_bestmode_age_limhi': {'col_name': 'mm_lowZ_bestmode_age_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_bestmode_ebv_limhi': {'col_name': 'mm_lowZ_bestmode_ebv_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_lowZ_bestmode_mass_limhi': {'col_name': 'mm_lowZ_bestmode_mass_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_youngestmode_age': {'col_name': 'mm_ZyoungallowhighEBV_youngestmode_age', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_youngestmode_ebv': {'col_name': 'mm_ZyoungallowhighEBV_youngestmode_ebv', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_youngestmode_mass': {'col_name': 'mm_ZyoungallowhighEBV_youngestmode_mass', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_youngestmode_likelihood': {'col_name': 'mm_ZyoungallowhighEBV_youngestmode_likelihood', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_bestmode_age': {'col_name': 'mm_ZyoungallowhighEBV_bestmode_age', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_bestmode_ebv': {'col_name': 'mm_ZyoungallowhighEBV_bestmode_ebv', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_bestmode_mass': {'col_name': 'mm_ZyoungallowhighEBV_bestmode_mass', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_bestmode_likelihood': {'col_name': 'mm_ZyoungallowhighEBV_bestmode_likelihood', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_youngestmode_age_limlo': {'col_name': 'mm_ZyoungallowhighEBV_youngestmode_age_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_youngestmode_ebv_limlo': {'col_name': 'mm_ZyoungallowhighEBV_youngestmode_ebv_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_youngestmode_mass_limlo': {'col_name': 'mm_ZyoungallowhighEBV_youngestmode_mass_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_bestmode_age_limlo': {'col_name': 'mm_ZyoungallowhighEBV_bestmode_age_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_bestmode_ebv_limlo': {'col_name': 'mm_ZyoungallowhighEBV_bestmode_ebv_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_bestmode_mass_limlo': {'col_name': 'mm_ZyoungallowhighEBV_bestmode_mass_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_youngestmode_age_limhi': {'col_name': 'mm_ZyoungallowhighEBV_youngestmode_age_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_youngestmode_ebv_limhi': {'col_name': 'mm_ZyoungallowhighEBV_youngestmode_ebv_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_youngestmode_mass_limhi': {'col_name': 'mm_ZyoungallowhighEBV_youngestmode_mass_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_bestmode_age_limhi': {'col_name': 'mm_ZyoungallowhighEBV_bestmode_age_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_bestmode_ebv_limhi': {'col_name': 'mm_ZyoungallowhighEBV_bestmode_ebv_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'mm_ZyoungallowhighEBV_bestmode_mass_limhi': {'col_name': 'mm_ZyoungallowhighEBV_bestmode_mass_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'OGC_fix': {'col_name': 'OGC_fix', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'OGC_spiraldigit': {'col_name': 'OGC_spiraldigit', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'OGCclass_okay': {'col_name': 'OGCclass_okay', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'OGCcc_okay': {'col_name': 'OGCcc_okay', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'envmask': {'col_name': 'envmask', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'YRO_fix': {'col_name': 'YRO_fix', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'YROclass_okay': {'col_name': 'YROclass_okay', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'YROcc_okay': {'col_name': 'YROcc_okay', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'highEBV_fix': {'col_name': 'highEBV_fix', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'highEBVclass_okay': {'col_name': 'highEBVclass_okay', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'highEBVcc_okay': {'col_name': 'highEBVcc_okay', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'SEDfix_age': {'col_name': 'SEDfix_age', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'SEDfix_ebv': {'col_name': 'SEDfix_ebv', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'SEDfix_mass': {'col_name': 'SEDfix_mass', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'SEDfix_age_limlo': {'col_name': 'SEDfix_age_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'SEDfix_ebv_limlo': {'col_name': 'SEDfix_ebv_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'SEDfix_mass_limlo': {'col_name': 'SEDfix_mass_limlo', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'SEDfix_age_limhi': {'col_name': 'SEDfix_age_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'SEDfix_ebv_limhi': {'col_name': 'SEDfix_ebv_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'SEDfix_mass_limhi': {'col_name': 'SEDfix_mass_limhi', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'SEDfix_likelihood_norm': {'col_name': 'SEDfix_likelihood_norm', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'PHANGS_DM': {'col_name': 'PHANGS_DM', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'PHANGS_DMPC': {'col_name': 'PHANGS_DMPC', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            'PHANGS_DERR': {'col_name': 'PHANGS_DERR', 'unit': None, 'unit_str': '', 'doc_comment': '', 'tab_comment': ''},
            
            
            #
            #
            # 'PHANGS_SEDFIX_CATEGORY_DR4': {# 'col_name': 'PHANGS_SEDFIX_CATEGORY_DR4',
            #                    'col_name': 'PHANGS_SEDFIX_CATEGORY_DR4',
            #                    'unit': None,
            #                    'unit_str': 'str',
            #                    'doc_comment': 'Category used to assign parameter grid for SED fit based on '
            #                                   'color-color topology and H-alpha surface brightness. '
            #                                   'Possible values can be `YRO\', `UNCHANGED\', or `OGC\'',
            #                    'tab_comment': 'Category used to assign parameter grid for SED fit. '
            #                                   'Possible values can be `YRO\', `UNCHANGED\', or `OGC\''
            #                    },
            #
            #
            #
            # 'sed_fix_category': {
            #                'col_name': 'sed_fix_category',
            #                'unit': None,
            #                'unit_str': 'str',
            #     'doc_comment': 'Category used to assign parameter grid for SED fit. Possible values can be '
            #                    '`YRO\', `UNCHANGED\', or `OGC\' ',
            #     'tab_comment': 'Category used to assign parameter grid for SED fit. Possible values can be '
            #                    '`YRO\', `UNCHANGED\', or `OGC\' ',
            #                },
            # 'SEDfix_age': {# 'col_name': 'SEDfix_age',
            #                'col_name': 'sed_fix_age',
            #                'unit': u.Myr,
            #                'unit_str': 'Myr',
            #     'doc_comment': 'age corresponding to fixed SED fit.',
            #     'tab_comment': 'age corresponding to fixed SED fit.',
            #                },
            # 'SEDfix_age_err': {# 'col_name': 'SEDfix_age_err',
            #                'col_name': 'sed_fix_age_err',
            #                'unit': u.Myr,
            #                'unit_str': 'Myr',
            #     'doc_comment':  'Uncertainty of `sed_fix_age\'',
            #     'tab_comment':  'Uncertainty of `sed_fix_age\'',
            #                },
            # 'SEDfix_ebv': {# 'col_name': 'SEDfix_ebv',
            #                'col_name': 'sed_fix_ebv',
            #                'unit': u.mag,
            #                'unit_str': 'mag',
            #     'doc_comment': 'Reddening E(B-V) corresponding to fixed SED fit.',
            #     'tab_comment': 'Reddening E(B-V) corresponding to fixed SED fit.',
            #                },
            # 'SEDfix_ebv_err': {'col_name': 'sed_fix_ebv_err',
            #                    'unit': u.mag,
            #                'unit_str': 'mag',
            #     'doc_comment':  'Uncertainty of `sed_fix_ebv\'',
            #     'tab_comment':  'Uncertainty of `sed_fix_ebv\'',
            #                },
            # 'SEDfix_mass': {
            #     'col_name': 'sed_fix_mass',
            #     'unit': u.M_sun,
            #     'unit_str': 'M_sun',
            #     'doc_comment': 'Stellar mass corresponding to fixed SED fit.',
            #     'tab_comment': 'Stellar mass corresponding to fixed SED fit.',
            #                 },
            # 'SEDfix_mass_err': {
            #                'col_name': 'sed_fix_mass_err',
            #                'unit': u.M_sun,
            #                'unit_str': 'M_sun',
            #     'doc_comment':  'Uncertainty of `sed_fix_mass\'',
            #     'tab_comment':  'Uncertainty of `sed_fix_mass\'',
            #                },
            # 'PHANGS_GALAXY':{
            #     'col_name': 'PHANGS_GALAXY',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'PHANGS_DMPC':{
            #     'col_name': 'PHANGS_DMPC',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'PHANGS_DERR':{
            #     'col_name': 'PHANGS_DERR',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            #
            # 'mm_Z_youngestmode_age': {
            #     'col_name': 'mm_Z_youngestmode_age',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_Z_youngestmode_ebv': {
            #     'col_name': 'mm_Z_youngestmode_ebv',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_Z_youngestmode_mass': {
            #     'col_name': 'mm_Z_youngestmode_mass',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_Z_bestmode_age': {
            #     'col_name': 'mm_Z_bestmode_age',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_Z_bestmode_ebv': {
            #     'col_name': 'mm_Z_bestmode_ebv',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_Z_bestmode_mass': {
            #     'col_name': 'mm_Z_bestmode_mass',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_lowZ_oldestmode_age': {
            #     'col_name': 'mm_lowZ_oldestmode_age',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_lowZ_oldestmode_ebv': {
            #     'col_name': 'mm_lowZ_oldestmode_ebv',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_lowZ_oldestmode_mass': {
            #     'col_name': 'mm_lowZ_oldestmode_mass',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_lowZ_bestmode_age': {
            #     'col_name': 'mm_lowZ_bestmode_age',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_lowZ_bestmode_ebv': {
            #     'col_name': 'mm_lowZ_bestmode_ebv',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_lowZ_bestmode_mass': {
            #     'col_name': 'mm_lowZ_bestmode_mass',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_ZhighEBV_youngestmode_age': {
            #     'col_name': 'mm_ZhighEBV_youngestmode_age',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_ZhighEBV_youngestmode_ebv': {
            #     'col_name': 'mm_ZhighEBV_youngestmode_ebv',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_ZhighEBV_youngestmode_mass': {
            #     'col_name': 'mm_ZhighEBV_youngestmode_mass',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_ZhighEBV_bestmode_age': {
            #     'col_name': 'mm_ZhighEBV_bestmode_age',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_ZhighEBV_bestmode_ebv': {
            #     'col_name': 'mm_ZhighEBV_bestmode_ebv',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # },
            # 'mm_ZhighEBV_bestmode_mass': {
            #     'col_name': 'mm_ZhighEBV_bestmode_mass',
            #     'unit': None,
            #     'unit_str': ' ',
            #     'doc_comment': ' ',
            #     'tab_comment': ' '
            # }
        }

        self.cand_identifier_columns = ['ID_PHANGS_CLUSTERS_v1p2',
                                        'ID_PHANGS_ALLSOURCES_v1p2',
                                        'PHANGS_X', 'PHANGS_Y', 'PHANGS_RA', 'PHANGS_DEC']

        self.identifier_columns = ['INDEX', 'ID_PHANGS_CLUSTERS_v1p2',
                                   # 'ID_PHANGS_ALLSOURCES',
                                   'PHANGS_X', 'PHANGS_Y', 'PHANGS_RA', 'PHANGS_DEC']

        self.classification_columns_doc = ['PHANGS_CLUSTER_CLASS_HUMAN', 'PHANGS_CLUSTER_CLASS_ML_VGG',
                                       'PHANGS_CLUSTER_CLASS_ML_VGG_QUAL', 'class_ml_vgg_corr']

        self.classification_columns = ['PHANGS_CLUSTER_CLASS_HUMAN', 'PHANGS_CLUSTER_CLASS_ML_VGG',
                                       'PHANGS_CLUSTER_CLASS_ML_VGG_QUAL']

        self.example_photometry_columns = ['PHANGS_[BAND]_VEGA_TOT', 'PHANGS_[BAND]_VEGA_TOT_ERR',
                                           'PHANGS_[BAND]_mJy_TOT', 'PHANGS_[BAND]_mJy_TOT_ERR']

        self.detect_shape_columns = ['PHANGS_NON_DETECTION_FLAG', 'PHANGS_NO_COVERAGE_FLAG', 'PHANGS_CI']

        # list of columns entering table 1 (Observational properties)
        self.tab1_columns = (['INDEX', 'id_phangs_cluster', 'ID_PHANGS_CLUSTERS_v1p2', 'ID_PHANGS_ALLSOURCES_v1p2',
                              'PHANGS_X', 'PHANGS_Y', 'PHANGS_RA', 'PHANGS_DEC'] +
                             self.classification_columns + self.example_photometry_columns +
                             self.detect_shape_columns + ['cc_class'])

        self.sed_columns_doc = ['PHANGS_AGE_MINCHISQ', 'PHANGS_AGE_MINCHISQ_ERR',
                            'PHANGS_MASS_MINCHISQ', 'PHANGS_MASS_MINCHISQ_ERR',
                            'PHANGS_EBV_MINCHISQ', 'PHANGS_EBV_MINCHISQ_ERR',
                            'PHANGS_REDUCED_MINCHISQ',
                            'NBHa_intensity', 'NBHa_intensity_medsub', 'NBHa_mask_medsub_lev1', 'NBHa_mask_medsub_lev2',
                                'NBHa_mask_medsub_lev3', 'NBHa_mask_medsub_lev5', 'NBHa_mask_medsub_lev10',
                                'NBHa_intensity_medsub_inseg', 'NBHa_HIIreg', 'OGCpoly_UBVI', 'FlagOGC_UBVI',
                                'OGC_Cull1', 'OGCpoly_BVI', 'FlagOGC_BVI', 'OGC_Cull2', 'OGCrpoly_UBVI',
                                'FlagOGCr_UBVI', 'OGCr_Cull1', 'OGCrpoly_BVI', 'FlagOGCr_BVI', 'OGCr_Cull2',
                                'YROpoly_UBVI', 'FlagYRO_UBVI', 'YRO_Cull1', 'YROpoly_BVI', 'FlagYRO_BVI', 'YRO_Cull2',
                                'mm_Zyoung_youngestmode_age', 'mm_Zyoung_youngestmode_ebv',
                                'mm_Zyoung_youngestmode_mass', 'mm_Zyoung_youngestmode_likelihood',
                                'mm_Zyoung_bestmode_age', 'mm_Zyoung_bestmode_ebv', 'mm_Zyoung_bestmode_mass',
                                'mm_Zyoung_bestmode_likelihood', 'mm_Zyoung_youngestmode_age_limlo',
                                'mm_Zyoung_youngestmode_ebv_limlo', 'mm_Zyoung_youngestmode_mass_limlo',
                                'mm_Zyoung_bestmode_age_limlo', 'mm_Zyoung_bestmode_ebv_limlo',
                                'mm_Zyoung_bestmode_mass_limlo', 'mm_Zyoung_youngestmode_age_limhi',
                                'mm_Zyoung_youngestmode_ebv_limhi', 'mm_Zyoung_youngestmode_mass_limhi',
                                'mm_Zyoung_bestmode_age_limhi', 'mm_Zyoung_bestmode_ebv_limhi',
                                'mm_Zyoung_bestmode_mass_limhi', 'mm_Z_youngestmode_age', 'mm_Z_youngestmode_ebv',
                                'mm_Z_youngestmode_mass', 'mm_Z_youngestmode_likelihood', 'mm_Z_bestmode_age',
                                'mm_Z_bestmode_ebv', 'mm_Z_bestmode_mass', 'mm_Z_bestmode_likelihood',
                                'mm_Z_youngestmode_age_limlo', 'mm_Z_youngestmode_ebv_limlo',
                                'mm_Z_youngestmode_mass_limlo', 'mm_Z_bestmode_age_limlo', 'mm_Z_bestmode_ebv_limlo',
                                'mm_Z_bestmode_mass_limlo', 'mm_Z_youngestmode_age_limhi',
                                'mm_Z_youngestmode_ebv_limhi', 'mm_Z_youngestmode_mass_limhi',
                                'mm_Z_bestmode_age_limhi', 'mm_Z_bestmode_ebv_limhi', 'mm_Z_bestmode_mass_limhi',
                                'mm_lowZ_oldestmode_age', 'mm_lowZ_oldestmode_ebv', 'mm_lowZ_oldestmode_mass',
                                'mm_lowZ_oldestmode_likelihood', 'mm_lowZ_bestmode_age', 'mm_lowZ_bestmode_ebv',
                                'mm_lowZ_bestmode_mass', 'mm_lowZ_bestmode_likelihood', 'mm_lowZ_oldestmode_age_limlo',
                                'mm_lowZ_oldestmode_ebv_limlo', 'mm_lowZ_oldestmode_mass_limlo',
                                'mm_lowZ_bestmode_age_limlo', 'mm_lowZ_bestmode_ebv_limlo',
                                'mm_lowZ_bestmode_mass_limlo', 'mm_lowZ_oldestmode_age_limhi',
                                'mm_lowZ_oldestmode_ebv_limhi', 'mm_lowZ_oldestmode_mass_limhi',
                                'mm_lowZ_bestmode_age_limhi', 'mm_lowZ_bestmode_ebv_limhi',
                                'mm_lowZ_bestmode_mass_limhi', 'mm_ZyoungallowhighEBV_youngestmode_age',
                                'mm_ZyoungallowhighEBV_youngestmode_ebv', 'mm_ZyoungallowhighEBV_youngestmode_mass',
                                'mm_ZyoungallowhighEBV_youngestmode_likelihood', 'mm_ZyoungallowhighEBV_bestmode_age',
                                'mm_ZyoungallowhighEBV_bestmode_ebv', 'mm_ZyoungallowhighEBV_bestmode_mass',
                                'mm_ZyoungallowhighEBV_bestmode_likelihood',
                                'mm_ZyoungallowhighEBV_youngestmode_age_limlo',
                                'mm_ZyoungallowhighEBV_youngestmode_ebv_limlo',
                                'mm_ZyoungallowhighEBV_youngestmode_mass_limlo',
                                'mm_ZyoungallowhighEBV_bestmode_age_limlo',
                                'mm_ZyoungallowhighEBV_bestmode_ebv_limlo',
                                'mm_ZyoungallowhighEBV_bestmode_mass_limlo',
                                'mm_ZyoungallowhighEBV_youngestmode_age_limhi',
                                'mm_ZyoungallowhighEBV_youngestmode_ebv_limhi',
                                'mm_ZyoungallowhighEBV_youngestmode_mass_limhi',
                                'mm_ZyoungallowhighEBV_bestmode_age_limhi',
                                'mm_ZyoungallowhighEBV_bestmode_ebv_limhi',
                                'mm_ZyoungallowhighEBV_bestmode_mass_limhi', 'OGC_fix', 'OGC_spiraldigit',
                                'OGCclass_okay', 'OGCcc_okay', 'envmask', 'YRO_fix', 'YROclass_okay', 'YROcc_okay',
                                'highEBV_fix', 'highEBVclass_okay', 'highEBVcc_okay', 'SEDfix_age', 'SEDfix_ebv',
                                'SEDfix_mass', 'SEDfix_age_limlo', 'SEDfix_ebv_limlo', 'SEDfix_mass_limlo',
                                'SEDfix_age_limhi', 'SEDfix_ebv_limhi', 'SEDfix_mass_limhi', 'SEDfix_likelihood_norm',
                                'PHANGS_DM', 'PHANGS_DMPC', 'PHANGS_DERR']

        self.sed_columns = ['PHANGS_AGE_MINCHISQ', 'PHANGS_AGE_MINCHISQ_ERR',
                            'PHANGS_MASS_MINCHISQ', 'PHANGS_MASS_MINCHISQ_ERR',
                            'PHANGS_EBV_MINCHISQ', 'PHANGS_EBV_MINCHISQ_ERR',
                            'PHANGS_REDUCED_MINCHISQ',
                            'NBHa_intensity', 'NBHa_intensity_medsub', 'NBHa_mask_medsub_lev1', 'NBHa_mask_medsub_lev2',
                                'NBHa_mask_medsub_lev3', 'NBHa_mask_medsub_lev5', 'NBHa_mask_medsub_lev10',
                                'NBHa_intensity_medsub_inseg', 'NBHa_HIIreg', 'OGCpoly_UBVI', 'FlagOGC_UBVI',
                                'OGC_Cull1', 'OGCpoly_BVI', 'FlagOGC_BVI', 'OGC_Cull2', 'OGCrpoly_UBVI',
                                'FlagOGCr_UBVI', 'OGCr_Cull1', 'OGCrpoly_BVI', 'FlagOGCr_BVI', 'OGCr_Cull2',
                                'YROpoly_UBVI', 'FlagYRO_UBVI', 'YRO_Cull1', 'YROpoly_BVI', 'FlagYRO_BVI', 'YRO_Cull2',
                                'mm_Zyoung_youngestmode_age', 'mm_Zyoung_youngestmode_ebv',
                                'mm_Zyoung_youngestmode_mass', 'mm_Zyoung_youngestmode_likelihood',
                                'mm_Zyoung_bestmode_age', 'mm_Zyoung_bestmode_ebv', 'mm_Zyoung_bestmode_mass',
                                'mm_Zyoung_bestmode_likelihood', 'mm_Zyoung_youngestmode_age_limlo',
                                'mm_Zyoung_youngestmode_ebv_limlo', 'mm_Zyoung_youngestmode_mass_limlo',
                                'mm_Zyoung_bestmode_age_limlo', 'mm_Zyoung_bestmode_ebv_limlo',
                                'mm_Zyoung_bestmode_mass_limlo', 'mm_Zyoung_youngestmode_age_limhi',
                                'mm_Zyoung_youngestmode_ebv_limhi', 'mm_Zyoung_youngestmode_mass_limhi',
                                'mm_Zyoung_bestmode_age_limhi', 'mm_Zyoung_bestmode_ebv_limhi',
                                'mm_Zyoung_bestmode_mass_limhi', 'mm_Z_youngestmode_age', 'mm_Z_youngestmode_ebv',
                                'mm_Z_youngestmode_mass', 'mm_Z_youngestmode_likelihood', 'mm_Z_bestmode_age',
                                'mm_Z_bestmode_ebv', 'mm_Z_bestmode_mass', 'mm_Z_bestmode_likelihood',
                                'mm_Z_youngestmode_age_limlo', 'mm_Z_youngestmode_ebv_limlo',
                                'mm_Z_youngestmode_mass_limlo', 'mm_Z_bestmode_age_limlo', 'mm_Z_bestmode_ebv_limlo',
                                'mm_Z_bestmode_mass_limlo', 'mm_Z_youngestmode_age_limhi',
                                'mm_Z_youngestmode_ebv_limhi', 'mm_Z_youngestmode_mass_limhi',
                                'mm_Z_bestmode_age_limhi', 'mm_Z_bestmode_ebv_limhi', 'mm_Z_bestmode_mass_limhi',
                                'mm_lowZ_oldestmode_age', 'mm_lowZ_oldestmode_ebv', 'mm_lowZ_oldestmode_mass',
                                'mm_lowZ_oldestmode_likelihood', 'mm_lowZ_bestmode_age', 'mm_lowZ_bestmode_ebv',
                                'mm_lowZ_bestmode_mass', 'mm_lowZ_bestmode_likelihood', 'mm_lowZ_oldestmode_age_limlo',
                                'mm_lowZ_oldestmode_ebv_limlo', 'mm_lowZ_oldestmode_mass_limlo',
                                'mm_lowZ_bestmode_age_limlo', 'mm_lowZ_bestmode_ebv_limlo',
                                'mm_lowZ_bestmode_mass_limlo', 'mm_lowZ_oldestmode_age_limhi',
                                'mm_lowZ_oldestmode_ebv_limhi', 'mm_lowZ_oldestmode_mass_limhi',
                                'mm_lowZ_bestmode_age_limhi', 'mm_lowZ_bestmode_ebv_limhi',
                                'mm_lowZ_bestmode_mass_limhi', 'mm_ZyoungallowhighEBV_youngestmode_age',
                                'mm_ZyoungallowhighEBV_youngestmode_ebv', 'mm_ZyoungallowhighEBV_youngestmode_mass',
                                'mm_ZyoungallowhighEBV_youngestmode_likelihood', 'mm_ZyoungallowhighEBV_bestmode_age',
                                'mm_ZyoungallowhighEBV_bestmode_ebv', 'mm_ZyoungallowhighEBV_bestmode_mass',
                                'mm_ZyoungallowhighEBV_bestmode_likelihood',
                                'mm_ZyoungallowhighEBV_youngestmode_age_limlo',
                                'mm_ZyoungallowhighEBV_youngestmode_ebv_limlo',
                                'mm_ZyoungallowhighEBV_youngestmode_mass_limlo',
                                'mm_ZyoungallowhighEBV_bestmode_age_limlo',
                                'mm_ZyoungallowhighEBV_bestmode_ebv_limlo',
                                'mm_ZyoungallowhighEBV_bestmode_mass_limlo',
                                'mm_ZyoungallowhighEBV_youngestmode_age_limhi',
                                'mm_ZyoungallowhighEBV_youngestmode_ebv_limhi',
                                'mm_ZyoungallowhighEBV_youngestmode_mass_limhi',
                                'mm_ZyoungallowhighEBV_bestmode_age_limhi',
                                'mm_ZyoungallowhighEBV_bestmode_ebv_limhi',
                                'mm_ZyoungallowhighEBV_bestmode_mass_limhi', 'OGC_fix', 'OGC_spiraldigit',
                                'OGCclass_okay', 'OGCcc_okay', 'envmask', 'YRO_fix', 'YROclass_okay', 'YROcc_okay',
                                'highEBV_fix', 'highEBVclass_okay', 'highEBVcc_okay', 'SEDfix_age', 'SEDfix_ebv',
                                'SEDfix_mass', 'SEDfix_age_limlo', 'SEDfix_ebv_limlo', 'SEDfix_mass_limlo',
                                'SEDfix_age_limhi', 'SEDfix_ebv_limhi', 'SEDfix_mass_limhi', 'SEDfix_likelihood_norm',
                                'PHANGS_DM', 'PHANGS_DMPC', 'PHANGS_DERR']

        self.tab2_columns = self.identifier_columns + self.sed_columns


        # colnames_for documentation
        self.doc_identifier_names = ['INDEX', 'id_phangs_cluster', 'ID_PHANGS_CLUSTERS_v1p2',
                                     'ID_PHANGS_ALLSOURCES_v1p2']
        self.doc_coord_names = ['PHANGS_X', 'PHANGS_Y', 'PHANGS_RA', 'PHANGS_DEC']

        self.doc_photometry_columns = ['PHANGS_F275W_VEGA_TOT', 'PHANGS_F275W_VEGA_TOT_ERR',
                                       'PHANGS_F336W_VEGA_TOT', 'PHANGS_F336W_VEGA_TOT_ERR',
                                       'PHANGS_F435W_VEGA_TOT', 'PHANGS_F435W_VEGA_TOT_ERR',
                                       'PHANGS_F438W_VEGA_TOT', 'PHANGS_F438W_VEGA_TOT_ERR',
                                       'PHANGS_F555W_VEGA_TOT', 'PHANGS_F555W_VEGA_TOT_ERR',
                                       'PHANGS_F814W_VEGA_TOT', 'PHANGS_F814W_VEGA_TOT_ERR',

                                       'PHANGS_F275W_mJy_TOT', 'PHANGS_F275W_mJy_TOT_ERR',
                                       'PHANGS_F336W_mJy_TOT', 'PHANGS_F336W_mJy_TOT_ERR',
                                       'PHANGS_F435W_mJy_TOT', 'PHANGS_F435W_mJy_TOT_ERR',
                                       'PHANGS_F438W_mJy_TOT', 'PHANGS_F438W_mJy_TOT_ERR',
                                       'PHANGS_F555W_mJy_TOT', 'PHANGS_F555W_mJy_TOT_ERR',
                                       'PHANGS_F814W_mJy_TOT', 'PHANGS_F814W_mJy_TOT_ERR']

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
        #columns_list += ['SEDfix_age', 'SEDfix_ebv', 'SEDfix_mass']

        return columns_list

    def get_cand_table_column_list(self, target):
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
        columns_list = self.cand_identifier_columns + self.classification_columns

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
        columns_list += self.sed_columns

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
            unit_str = self.cat_info[key]['unit_str'].replace('_', '\\_')
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
