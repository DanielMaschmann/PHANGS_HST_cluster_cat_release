"""
All production routines for the final PHANGS-HST cluster catalog data release will be gathered here
"""

import os
from datetime import datetime
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits, ascii
from astropy.table import Column, hstack
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import ConvexHull

from cat_info import CatalogInfo
import helper_func
from visualization_tool import PhotVisualize


class DataReleaseRoutines(CatalogInfo):
    def __init__(self, hst_cc_ver, path2ir, catalog_output_path,
                 path2artifact=None, artifact_removal_flag=False, artifact_rais_file_not_found_flag=True,
                 plot_removed_artifacts_flag=False, data_release=4, catalog_release=2,
                 existing_artifact_removal_flag=False,
                 path2questionable_artifacts=None,
                 path2diffraction_spike_masks=None,
                 v_i_color_lim=None, ci_lim=None,
                 hst_data_path=None, nircam_data_path=None, miri_data_path=None,
                 hst_data_ver='v1', nircam_data_ver='v0p4p2', miri_data_ver='v0p5'):

        # add input keys to attributes
        self.hst_cc_ver = hst_cc_ver
        self.path2ir = path2ir
        self.catalog_output_path = catalog_output_path
        self.path2artifact = path2artifact
        self.artifact_removal_flag = artifact_removal_flag
        self.artifact_rais_file_not_found_flag = artifact_rais_file_not_found_flag
        self.plot_removed_artifacts_flag = plot_removed_artifacts_flag
        self.data_release = data_release
        self.catalog_release = catalog_release
        self.existing_artifact_removal_flag = existing_artifact_removal_flag
        self.path2questionable_artifacts = path2questionable_artifacts
        self.path2diffraction_spike_masks = path2diffraction_spike_masks
        self.v_i_color_lim = v_i_color_lim
        self.ci_lim = ci_lim

        self.questionable_artefact_table = None

        # for data visualization
        self.hst_data_path = hst_data_path
        self.nircam_data_path = nircam_data_path
        self.miri_data_path = miri_data_path
        self.hst_data_ver = hst_data_ver
        self.nircam_data_ver = nircam_data_ver
        self.miri_data_ver = miri_data_ver


        # correction factor for aperture correction bug in ngc 1512 and ngc 1510
        self.ngc1512_app_corr_offset_mag = 1.724
        self.ngc1512_app_corr_offset_flux = 10**(-0.4*self.ngc1512_app_corr_offset_mag)


        # load constructor of parent class
        super().__init__()

    def create_final_data_release(self):
        """
        Function to run all steps to create the PHANGs-HST cluster catalog data release
        """
        # create catalogs
        # loop over all targets for which a cluster catalog exists

        removal_statistics_dict = {}

        for target in self.phangs_hst_cluster_cat_target_list:

            # get artifact removal table
            table_artifact = self.get_artifact_cat(target=target)
            # get table with second classification
            table_re_classified = self.get_second_classify_cat(target=target)
            # get table for diffraction spikes
            if target == 'ngc1512':
                diffraction_spike_mask, diffraction_spike_wcs = self.load_diffraction_spike_masks(target=target, target_str='ngc1512c')
            else:
                diffraction_spike_mask, diffraction_spike_wcs = self.load_diffraction_spike_masks(target=target, target_str=target)

            if (target == 'ngc0628e') | (target == 'ngc0628c'):
                diffraction_spike_mask_2, diffraction_spike_wcs_2 = (
                    self.load_diffraction_spike_masks(target='ngc0628', target_str='ngc0628'))
            else:
                diffraction_spike_mask_2 = None
                diffraction_spike_wcs_2 = None

            # create candidate table
            candidate_table_ir = self.get_ir_cat(target, classify=None, cl_class='candidates')
            # find cross-matching artifact
            for artifact_index in range(len(table_artifact)):
                artifact_in_table_ir = ((candidate_table_ir['ID_PHANGS_CLUSTERS_v1p2'] ==
                                         table_artifact['ID_PHANGS_CLUSTERS_v1p2'][artifact_index]) &
                                        (candidate_table_ir['PHANGS_X'] ==
                                         table_artifact['PHANGS_X'][artifact_index]) &
                                        (candidate_table_ir['PHANGS_Y'] ==
                                         table_artifact['PHANGS_Y'][artifact_index]) &
                                        (candidate_table_ir['PHANGS_RA'] ==
                                         table_artifact['PHANGS_RA'][artifact_index]) &
                                        (candidate_table_ir['PHANGS_DEC'] ==
                                         table_artifact['PHANGS_DEC'][artifact_index]))
                artifact_in_table_re_classified = \
                    ((table_re_classified['ID_PHANGS_CLUSTERS_v1p2'] ==
                      table_artifact['ID_PHANGS_CLUSTERS_v1p2'][artifact_index]))
                # update human classification:
                # print('BCW classification: ',
                # candidate_table_ir['PHANGS_CLUSTER_CLASS_HUMAN'][artifact_in_table_ir][0])
                if np.invert(np.isnan(candidate_table_ir['PHANGS_CLUSTER_CLASS_HUMAN'][artifact_in_table_ir][0])):
                    # print('already classified')
                    continue
                # print('Chris classification: ', table_artifact['NEW_CLASS'][artifact_index])
                # change the human classification
                if table_artifact['NEW_CLASS'][artifact_index] != -999:
                    candidate_table_ir['PHANGS_CLUSTER_CLASS_HUMAN'][artifact_in_table_ir] = (
                        table_artifact)['NEW_CLASS'][artifact_index]
                else:
                    # # print('BCW second classification',
                    # #       int(table_re_classified['BCW_estimate'][artifact_in_table_re_classified][0]))
                    # # if there is no data in the column, the type is masked
                    if not np.ma.isMaskedArray(table_re_classified['BCW_estimate'][artifact_in_table_re_classified][0]):
                        if int(table_re_classified['BCW_estimate'][artifact_in_table_re_classified][0]) != -999:
                            # print('reclassified with ',
                            #       int(table_re_classified['BCW_estimate'][artifact_in_table_re_classified][0]))
                            candidate_table_ir['PHANGS_CLUSTER_CLASS_HUMAN'][artifact_in_table_ir] = \
                                int(table_re_classified['BCW_estimate'][artifact_in_table_re_classified][0])

            # update diffraction spike classes
            if diffraction_spike_mask is not None:
                # x_pixel_coords = np.array(np.rint(candidate_table_ir['PHANGS_X']), dtype=int)
                # y_pixel_coords = np.array(np.rint(candidate_table_ir['PHANGS_Y']), dtype=int)
                ra_pixel_coords = candidate_table_ir['PHANGS_RA']
                dec_pixel_coords = candidate_table_ir['PHANGS_DEC']
                pos = SkyCoord(ra=ra_pixel_coords, dec=dec_pixel_coords, unit=(u.degree, u.degree), frame='fk5')
                pos_pix = diffraction_spike_wcs.world_to_pixel(pos)
                x_pixel_coords = np.array(np.rint(pos_pix[0]), dtype=int)
                y_pixel_coords = np.array(np.rint(pos_pix[1]), dtype=int)
                mask_covered_coordinates = ((x_pixel_coords > 0) & (y_pixel_coords > 0) &
                                            (x_pixel_coords < diffraction_spike_mask.shape[0]) &
                                            (y_pixel_coords < diffraction_spike_mask.shape[1]))
                artifact_in_diffraction_spike = (
                        diffraction_spike_mask[y_pixel_coords[mask_covered_coordinates],
                                               x_pixel_coords[mask_covered_coordinates]] > 0)
                candidate_table_ir['PHANGS_CLUSTER_CLASS_HUMAN'][mask_covered_coordinates][artifact_in_diffraction_spike] = 8
            if diffraction_spike_mask_2 is not None:
                ra_pixel_coords = candidate_table_ir['PHANGS_RA']
                dec_pixel_coords = candidate_table_ir['PHANGS_DEC']
                pos = SkyCoord(ra=ra_pixel_coords, dec=dec_pixel_coords, unit=(u.degree, u.degree), frame='fk5')
                pos_pix = diffraction_spike_wcs_2.world_to_pixel(pos)
                x_pixel_coords = np.array(np.rint(pos_pix[0]), dtype=int)
                y_pixel_coords = np.array(np.rint(pos_pix[1]), dtype=int)
                mask_covered_coordinates = ((x_pixel_coords > 0) & (y_pixel_coords > 0) &
                                            (x_pixel_coords < diffraction_spike_mask_2.shape[0]) &
                                            (y_pixel_coords < diffraction_spike_mask_2.shape[1]))
                artifact_in_diffraction_spike = (
                        diffraction_spike_mask_2[y_pixel_coords[mask_covered_coordinates],
                                                 x_pixel_coords[mask_covered_coordinates]] > 0)
                candidate_table_ir['PHANGS_CLUSTER_CLASS_HUMAN'][mask_covered_coordinates][artifact_in_diffraction_spike] = 8
            # update very red stars
            vi_color = candidate_table_ir['PHANGS_F555W_vega_tot'] - candidate_table_ir['PHANGS_F814W_vega_tot']
            ci = candidate_table_ir['PHANGS_CI']
            very_red_star_mask_candidate_table_ir = (vi_color > self.v_i_color_lim) & (ci < self.ci_lim)
            candidate_table_ir['PHANGS_CLUSTER_CLASS_HUMAN'][very_red_star_mask_candidate_table_ir] = 19

            # save candidate table
            # create candidate table
            candidate_table = self.create_cand_table(target=target, table=candidate_table_ir)
            # sort them by increasing Y pixel
            increase_y_sort = np.argsort(candidate_table['PHANGS_Y'])
            candidate_table = candidate_table[increase_y_sort]
            # get ids from all objects which have been classified as a cluster
            indexes_with_clusters = np.where(
                # can be human cluster
                (candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] == 1) |
                (candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] == 2) |
                (candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] == 3) |
                # or ML cluster
                (((candidate_table['PHANGS_CLUSTER_CLASS_ML_VGG'] == 1) |
                  (candidate_table['PHANGS_CLUSTER_CLASS_ML_VGG'] == 2) |
                  (candidate_table['PHANGS_CLUSTER_CLASS_ML_VGG'] == 3)) &
                 # but must be no human classified artefact otherwise this object is no cluster and will be sorted out
                 ((candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] < 4) |
                  np.isnan(candidate_table['PHANGS_CLUSTER_CLASS_HUMAN']))))
            # create data for column
            cluster_id_column_data = np.ones(len(candidate_table), dtype=int) * -999
            # only give increasing values for clusters
            cluster_id_column_data[indexes_with_clusters] = np.arange(1, len(indexes_with_clusters[0])+1)
            # create cluster id column
            cluster_id_column = Column(data=cluster_id_column_data, name=self.cat_info['id_phangs_cluster']['col_name'],
                                       dtype=cluster_id_column_data.dtype,
                                       description=self.cat_info['id_phangs_cluster']['doc_comment'])
            # get index column
            index_column = Column(data=np.arange(1, len(candidate_table)),
                                  name=self.cat_info['INDEX']['col_name'],
                                  dtype=int,
                                  description=self.cat_info['INDEX']['doc_comment'])

            candidate_table = hstack([index_column, cluster_id_column, candidate_table])

            # add corrected ML class column
            ml_class_corr_data = candidate_table['PHANGS_CLUSTER_CLASS_ML_VGG']
            mask_identified_artefacts_cand = candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] > 3
            ml_class_corr_data[mask_identified_artefacts_cand] = candidate_table[mask_identified_artefacts_cand]['PHANGS_CLUSTER_CLASS_HUMAN']
            ml_class_corr_column = Column(data=ml_class_corr_data,
                                          name=self.cat_info['class_ml_vgg_corr']['col_name'],
                                          dtype=int,
                                          description=self.cat_info['class_ml_vgg_corr']['doc_comment'])
            index_last_classification = np.where(np.array(candidate_table.colnames) == 'PHANGS_CLUSTER_CLASS_ML_VGG_QUAL')
            candidate_table.add_column(ml_class_corr_column, index=index_last_classification[0][0]+1)

            # set SED fits of not covered objects to -999
            mask_bad_coverage_candidates = ((candidate_table['PHANGS_NON_DETECTION_FLAG'] >= 2) |
                                            (candidate_table['PHANGS_NO_COVERAGE_FLAG'] >= 2))
            list_names_sed_fixes = ['SEDfix_age', 'SEDfix_ebv', 'SEDfix_mass', 'SEDfix_age_limlo', 'SEDfix_ebv_limlo',
                              'SEDfix_mass_limlo', 'SEDfix_age_limhi', 'SEDfix_ebv_limhi', 'SEDfix_mass_limhi']
            for names_sed_fixes in list_names_sed_fixes:
                candidate_table[names_sed_fixes][mask_bad_coverage_candidates] = -999

            # save table
            # check if folder already exists
            if not os.path.isdir(self.catalog_output_path):
                os.mkdir(self.catalog_output_path)
            if not os.path.isdir(self.catalog_output_path + '/catalogs'):
                os.mkdir(self.catalog_output_path + '/catalogs')
            # get table name
            table_name_candidate = self.get_data_release_table_name(target=target, classify=None, cl_class=None,
                                                                    table_type='cand')
            # save table
            candidate_table.write(Path(self.catalog_output_path + '/catalogs') / Path(table_name_candidate), overwrite=True)

            # if plotting is enabled load photo visualization access
            if self.plot_removed_artifacts_flag:
                visualization_access = PhotVisualize(
                    target_name=target,
                    hst_data_path=self.hst_data_path,
                    nircam_data_path=self.nircam_data_path,
                    miri_data_path=self.miri_data_path,
                    hst_data_ver=self.hst_data_ver,
                    nircam_data_ver=self.nircam_data_ver,
                    miri_data_ver=self.miri_data_ver
                )
                visualization_access.load_hst_nircam_miri_bands(flux_unit='MJy/sr', load_err=False)
            else:
                visualization_access = None

            # loop over human and ML catalogs
            for classify in ['human', 'ml']:
                # loop over class 1+2 and class 3 tables
                for cl_class in ['class12', 'class3']:
                    print('creating cluster catalog for ', target, classify, cl_class)
                    # load the internal release
                    table_ir = self.get_ir_cat(target, classify, cl_class)
                    # create an empty artifact mask
                    artifact_mask_table_ir = np.zeros(len(table_ir), dtype=bool)
                    if self.artifact_removal_flag:
                        used_artifacts_mask_table_artifact = np.zeros(len(table_artifact), dtype=bool)
                        # if the artifact removal flag is set we will loop over all manually IDed artifacts.
                        for artifact_index in range(len(table_artifact)):
                            # now we only use the identified object, when the human classification is not 1,2 or 3
                            if table_artifact['PHANGS_CLUSTER_CLASS_HUMAN'][artifact_index] in [1, 2, 3]:
                                continue
                            # find cross-matching artifacts
                            artifact_in_table_ir = ((table_ir['ID_PHANGS_CLUSTERS_v1p2'] ==
                                                     table_artifact['ID_PHANGS_CLUSTERS_v1p2'][artifact_index]) &
                                                    (table_ir['PHANGS_X'] ==
                                                     table_artifact['PHANGS_X'][artifact_index]) &
                                                    (table_ir['PHANGS_Y'] ==
                                                     table_artifact['PHANGS_Y'][artifact_index]) &
                                                    (table_ir['PHANGS_RA'] ==
                                                     table_artifact['PHANGS_RA'][artifact_index]) &
                                                    (table_ir['PHANGS_DEC'] ==
                                                     table_artifact['PHANGS_DEC'][artifact_index]))
                            artifact_in_table_re_classified = \
                                ((table_re_classified['ID_PHANGS_CLUSTERS_v1p2'] ==
                                  table_artifact['ID_PHANGS_CLUSTERS_v1p2'][artifact_index]))
                            # now check if the artefact index was re-evaluated. otherwise directly add them to the mask
                            if np.sum(np.array(artifact_in_table_re_classified, dtype=bool)) > 0:
                                bcw_classify = table_re_classified[artifact_in_table_re_classified]['BCW_estimate']
                                # check if BCW classified it as a cluster
                                if bcw_classify in [1, 2, 3, -999]:
                                    print('BCW: No artefact')
                                    continue
                                # # check if MF classified them as a cluster if BCW did not classified them
                                # elif bcw_classify == -999:
                                #     mf_classify = table_re_classified[artifact_in_table_re_classified]['NEW_CLASS']
                                #     if mf_classify in ['1', '2', '3', '1/2', '2/3']:
                                #         print('MF: No artefact')
                                #         continue
                                #     else:
                                #         print('MF: artefact')
                                #         artifact_mask_table_ir += artifact_in_table_ir
                                #         # add identified artifacts into a mask for plotting
                                #         if np.sum(np.array(artifact_in_table_ir, dtype=int)) > 0:
                                #             used_artifacts_mask_table_artifact[artifact_index] = True
                                else:
                                    print('BCW: artefact')
                                    artifact_mask_table_ir += artifact_in_table_ir
                                    # add identified artifacts into a mask for plotting
                                    if np.sum(np.array(artifact_in_table_ir, dtype=int)) > 0:
                                        used_artifacts_mask_table_artifact[artifact_index] = True
                            else:
                                artifact_mask_table_ir += artifact_in_table_ir
                                # add identified artifacts into a mask for plotting
                                if np.sum(np.array(artifact_in_table_ir, dtype=int)) > 0:
                                    used_artifacts_mask_table_artifact[artifact_index] = True

                        # flagg objects inside the diffraction spikes
                        if diffraction_spike_mask is not None:
                            ra_pixel_coords = table_ir['PHANGS_RA']
                            dec_pixel_coords = table_ir['PHANGS_DEC']
                            pos = SkyCoord(ra=ra_pixel_coords, dec=dec_pixel_coords, unit=(u.degree, u.degree), frame='fk5')
                            pos_pix = diffraction_spike_wcs.world_to_pixel(pos)
                            x_pixel_coords = np.array(np.rint(pos_pix[0]), dtype=int)
                            y_pixel_coords = np.array(np.rint(pos_pix[1]), dtype=int)
                            mask_covered_coordinates = ((x_pixel_coords > 0) & (y_pixel_coords > 0) &
                                                        (x_pixel_coords < diffraction_spike_mask.shape[0]) &
                                                        (y_pixel_coords < diffraction_spike_mask.shape[1]))
                            artifact_in_diffraction_spike = np.zeros(len(table_ir['PHANGS_X']), dtype=bool)

                            artifact_in_diffraction_spike[mask_covered_coordinates] = \
                                (diffraction_spike_mask[y_pixel_coords[mask_covered_coordinates],
                                                        x_pixel_coords[mask_covered_coordinates]] > 0)

                            # x_pixel_coords = np.array(np.rint(table_ir['PHANGS_X']), dtype=int)
                            # y_pixel_coords = np.array(np.rint(table_ir['PHANGS_Y']), dtype=int)
                            # artifact_in_diffraction_spike = diffraction_spike_mask[y_pixel_coords, x_pixel_coords] > 0
                        else:
                            artifact_in_diffraction_spike = np.zeros(len(table_ir['PHANGS_X']), dtype=bool)

                        if diffraction_spike_mask_2 is not None:
                            ra_pixel_coords = table_ir['PHANGS_RA']
                            dec_pixel_coords = table_ir['PHANGS_DEC']
                            pos = SkyCoord(ra=ra_pixel_coords, dec=dec_pixel_coords, unit=(u.degree, u.degree), frame='fk5')
                            pos_pix = diffraction_spike_wcs_2.world_to_pixel(pos)
                            x_pixel_coords = np.array(np.rint(pos_pix[0]), dtype=int)
                            y_pixel_coords = np.array(np.rint(pos_pix[1]), dtype=int)
                            mask_covered_coordinates = ((x_pixel_coords > 0) & (y_pixel_coords > 0) &
                                                        (x_pixel_coords < diffraction_spike_mask_2.shape[0]) &
                                                        (y_pixel_coords < diffraction_spike_mask_2.shape[1]))

                            artifact_in_diffraction_spike[mask_covered_coordinates] += \
                                (diffraction_spike_mask_2[y_pixel_coords[mask_covered_coordinates],
                                                          x_pixel_coords[mask_covered_coordinates]] > 0)

                        # plot artifacts
                        if self.plot_removed_artifacts_flag:
                            # check if photometry path were given
                            if ((self.hst_data_path is None) & (self.nircam_data_path is None) &
                                    (self.miri_data_path is None)):
                                print('No Path was given for cluster visualization')
                            else:
                                for artifact_index in range(len(table_artifact)):
                                    if not used_artifacts_mask_table_artifact[artifact_index]:
                                        continue
                                    # if an object was already classified as an artefact we do not need to plot it again
                                    if table_artifact['PHANGS_CLUSTER_CLASS_HUMAN'][artifact_index] != -999:
                                        continue
                                    ra = table_artifact['PHANGS_RA'][artifact_index]
                                    dec = table_artifact['PHANGS_DEC'][artifact_index]
                                    cluster_id = table_artifact['ID_PHANGS_CLUSTERS_v1p2'][artifact_index]
                                    new_class = table_artifact['NEW_CLASS'][artifact_index]
                                    ml_class = table_artifact['PHANGS_CLUSTER_CLASS_ML_VGG'][artifact_index]
                                    ml_class_qual = (
                                        table_artifact)['PHANGS_CLUSTER_CLASS_ML_VGG_QUAL'][artifact_index]
                                    hum_class = table_artifact['PHANGS_CLUSTER_CLASS_HUMAN'][artifact_index]

                                    color_vi = (table_artifact['PHANGS_F555W_VEGA_TOT'][artifact_index] -
                                                table_artifact['PHANGS_F814W_VEGA_TOT'][artifact_index])

                                    str_line_1 = ('%s, '
                                                  'ID_PHANGS_CLUSTERS_v1p2 = %i, '
                                                  'HUM_CLASS = %s, '
                                                  'ML_CLASS = %s, '
                                                  'ML_CLASS_QUAL = %.1f, '
                                                  'NEW_CLASS = %s '
                                                  % (target, cluster_id, hum_class, ml_class, ml_class_qual,
                                                     new_class))
                                    # add V-I color to the
                                    str_line_2 = ('V-I = %.2f' % color_vi)

                                    fig = visualization_access.plot_multi_band_artifacts(ra=ra, dec=dec,
                                                                                         str_line_1=str_line_1,
                                                                                         str_line_2=str_line_2)
                                    # plot for each galaxy
                                    if not os.path.isdir('plot_output/' + target):
                                        os.makedirs('plot_output/' + target)
                                    fig.savefig('plot_output/' + target + '/artifact_%i.png' % cluster_id)

                                    # plot if this object is not yet decided
                                    if (new_class == -999) & (hum_class == -999):
                                        if not os.path.isdir('plot_output/unclear'):
                                            os.makedirs('plot_output/unclear')
                                        fig.savefig('plot_output/unclear' + '/artifact_%s_%i.png' %
                                                    (target, cluster_id))

                                    plt.clf()
                                    plt.close("all")

                    print('number cross matched artefacts ', sum(artifact_mask_table_ir))
                    print('number objects in diffraction spikes', sum(artifact_in_diffraction_spike))

                    if self.existing_artifact_removal_flag:
                        # select all artifacts which have been already classified but still can occur in the ML sample
                        existing_artifact_mask_table_ir = table_ir['PHANGS_CLUSTER_CLASS_HUMAN'] > 3
                        print('existing_artifact_mask_table_ir', sum(existing_artifact_mask_table_ir))
                    else:
                        existing_artifact_mask_table_ir = np.zeros(len(table_ir), dtype=bool)

                    # select stars which are too red
                    vi_color = table_ir['PHANGS_F555W_vega_tot'] - table_ir['PHANGS_F814W_vega_tot']
                    ci = table_ir['PHANGS_CI']
                    very_red_star_mask_table_ir = (vi_color > self.v_i_color_lim) & (ci < self.ci_lim)
                    print('number of red stars (V-I > %.1f & CI < %.1f)' % (self.v_i_color_lim, self.ci_lim),
                          sum(very_red_star_mask_table_ir))

                    print('total artefact removal: ', sum(artifact_mask_table_ir + existing_artifact_mask_table_ir +
                                                          artifact_in_diffraction_spike + very_red_star_mask_table_ir))
                    if target in removal_statistics_dict.keys():
                        removal_statistics_dict[target].update({
                            'first_insp_%s_%s' % (classify, cl_class): sum(existing_artifact_mask_table_ir),
                            'second_insp_%s_%s' % (classify, cl_class): sum(artifact_mask_table_ir * np.invert(existing_artifact_mask_table_ir)),
                            'red_star_%s_%s' % (classify, cl_class): sum(very_red_star_mask_table_ir * np.invert(existing_artifact_mask_table_ir + artifact_mask_table_ir)),
                            'diff_%s_%s' % (classify, cl_class): sum(artifact_in_diffraction_spike * np.invert(very_red_star_mask_table_ir + existing_artifact_mask_table_ir + artifact_mask_table_ir)),
                        })
                    else:
                        removal_statistics_dict.update({target: {
                            'first_insp_%s_%s' % (classify, cl_class): sum(existing_artifact_mask_table_ir),
                            'second_insp_%s_%s' % (classify, cl_class): sum(artifact_mask_table_ir * np.invert(existing_artifact_mask_table_ir)),
                            'red_star_%s_%s' % (classify, cl_class): sum(very_red_star_mask_table_ir * np.invert(existing_artifact_mask_table_ir + artifact_mask_table_ir)),
                            'diff_%s_%s' % (classify, cl_class): sum(artifact_in_diffraction_spike * np.invert(very_red_star_mask_table_ir + existing_artifact_mask_table_ir + artifact_mask_table_ir)),
                        }})

                    # apply artifact mask
                    table_ir = table_ir[np.invert(artifact_mask_table_ir + existing_artifact_mask_table_ir +
                                                  artifact_in_diffraction_spike + very_red_star_mask_table_ir)]

                    # make exception for NGC 1512
                    if target == 'ngc1512':
                        x_pixel = table_ir['PHANGS_X']
                        y_pixel = table_ir['PHANGS_Y']

                        # center_ngc1510_x = 9340
                        # center_ngc1510_y = 3320
                        # selection_rad = 1000
                        # mask_ngc1510_select = (np.sqrt((x_pixel-center_ngc1510_x)**2 +
                        #                                (y_pixel-center_ngc1510_y)**2) < selection_rad)

                        # get straight line to divide the points
                        x1 = 6107
                        y1 = 2580
                        x2 = 9408
                        y2 = 7010
                        slope = (y1-y2)/(x1-x2)
                        intersect = y1 - slope * x1
                        mask_ngc1510_select = y_pixel < slope * x_pixel + intersect

                        # plt.plot([x1, x2], [y1, y2], color='k', linestyle='--')
                        # plt.scatter(x_pixel, y_pixel)
                        # plt.scatter(x_pixel[mask_ngc1510_select], y_pixel[mask_ngc1510_select])
                        # plt.show()

                        table_ir_ngc1512 = table_ir[np.invert(mask_ngc1510_select)]
                        table_ir_ngc1510 = table_ir[mask_ngc1510_select]
                        obs_table_1 = self.create_obs_table(target='ngc1512', table=table_ir_ngc1512,
                                                            cand_table=candidate_table, classify=classify)
                        obs_table_2 = self.create_obs_table(target='ngc1510', table=table_ir_ngc1510,
                                                            cand_table=candidate_table, classify=classify)
                        sed_table_1 = self.create_sed_table(target=target,
                                                            table=table_ir_ngc1512,
                                                            cand_table=candidate_table)
                        sed_table_2 = self.create_sed_table(target=target,
                                                            table=table_ir_ngc1510,
                                                            cand_table=candidate_table)
                        obs_table_1['INDEX'] = range(1, len(obs_table_1['INDEX'])+1)
                        obs_table_2['INDEX'] = range(1, len(obs_table_2['INDEX'])+1)

                        sed_table_1['INDEX'] = range(1, len(sed_table_1['INDEX'])+1)
                        sed_table_2['INDEX'] = range(1, len(sed_table_2['INDEX'])+1)
                                            # get table name
                        table_name_obs_1 = self.get_data_release_table_name(target='ngc1512', classify=classify,
                                                                            cl_class=cl_class)
                        table_name_obs_2 = self.get_data_release_table_name(target='ngc1510', classify=classify,
                                                                            cl_class=cl_class)

                        table_name_sed_1 = self.get_data_release_table_name(target='ngc1512', classify=classify,
                                                                            cl_class=cl_class, table_type='sed')
                        table_name_sed_2 = self.get_data_release_table_name(target='ngc1510', classify=classify,
                                                                            cl_class=cl_class, table_type='sed')
                    else:
                        obs_table_1 = self.create_obs_table(target=target, table=table_ir,
                                                            cand_table=candidate_table, classify=classify)
                        sed_table_1 = self.create_sed_table(target=target,
                                                            table=table_ir,
                                                            cand_table=candidate_table)
                        obs_table_2 = None
                        sed_table_2 = None

                        table_name_obs_1 = self.get_data_release_table_name(target=target, classify=classify,
                                                                            cl_class=cl_class)
                        table_name_sed_1 = self.get_data_release_table_name(target=target, classify=classify,
                                                                            cl_class=cl_class, table_type='sed')
                        table_name_obs_2 = None
                        table_name_sed_2 = None

                    # set SED fits of not covered objects to -999
                    list_names_sed_fixes = ['SEDfix_age', 'SEDfix_ebv', 'SEDfix_mass', 'SEDfix_age_limlo',
                                            'SEDfix_ebv_limlo', 'SEDfix_mass_limlo', 'SEDfix_age_limhi',
                                            'SEDfix_ebv_limhi', 'SEDfix_mass_limhi']
                    mask_bad_coverage_sed = ((obs_table_1['PHANGS_NON_DETECTION_FLAG'] >= 2) |
                                             (obs_table_1['PHANGS_NO_COVERAGE_FLAG'] >= 2))
                    for names_sed_fixes in list_names_sed_fixes:
                        sed_table_1[names_sed_fixes][mask_bad_coverage_sed] = -999

                    if sed_table_2 is not None:
                        mask_bad_coverage_sed = ((obs_table_2['PHANGS_NON_DETECTION_FLAG'] >= 2) |
                                                 (obs_table_2['PHANGS_NO_COVERAGE_FLAG'] >= 2))
                        for names_sed_fixes in list_names_sed_fixes:
                            sed_table_2[names_sed_fixes][mask_bad_coverage_sed] = -999

                    # save table
                    # check if folder already exists
                    if not os.path.isdir(self.catalog_output_path + '/catalogs'):
                        os.mkdir(self.catalog_output_path + '/catalogs')

                    # save table
                    obs_table_1.write(Path(self.catalog_output_path + '/catalogs') / Path(table_name_obs_1), overwrite=True)
                    if obs_table_2 is not None:
                        obs_table_2.write(Path(self.catalog_output_path + '/catalogs') / Path(table_name_obs_2), overwrite=True)

                    sed_table_1.write(Path(self.catalog_output_path + '/catalogs') / Path(table_name_sed_1), overwrite=True)
                    if sed_table_2 is not None:
                        sed_table_2.write(Path(self.catalog_output_path + '/catalogs') / Path(table_name_sed_2), overwrite=True)


                    # # create observation table
                    # obs_table = self.create_obs_table(target=target, table=table_ir, cand_table=candidate_table,
                    #                                   classify=classify)
                    #
                    # sed_table = self.create_sed_table(table=table_ir, cand_table=candidate_table)
                    # sed_table.write(Path(self.catalog_output_path + '/catalogs') / Path(table_name_sed), overwrite=True)


        print(removal_statistics_dict)
        for target in self.phangs_hst_cluster_cat_target_list:
            print(target, 'hum_class12',
                  'first-insp:', removal_statistics_dict[target]['first_insp_human_class12'],
                  'second-insp: ', removal_statistics_dict[target]['second_insp_human_class12'],
                  'red-star: ', removal_statistics_dict[target]['red_star_human_class12'],
                  'diffraction: ', removal_statistics_dict[target]['diff_human_class12'])
            print(target, 'hum_class3',
                  'first-insp:', removal_statistics_dict[target]['first_insp_human_class3'],
                  'second-insp: ', removal_statistics_dict[target]['second_insp_human_class3'],
                  'red-star: ', removal_statistics_dict[target]['red_star_human_class3'],
                  'diffraction: ', removal_statistics_dict[target]['diff_human_class3'])
            print(target, 'ml_class12',
                  'first-insp:', removal_statistics_dict[target]['first_insp_ml_class12'],
                  'second-insp: ', removal_statistics_dict[target]['second_insp_ml_class12'],
                  'red-star: ', removal_statistics_dict[target]['red_star_ml_class12'],
                  'diffraction: ', removal_statistics_dict[target]['diff_ml_class12'])
            print(target, 'ml_class3',
                  'first-insp:', removal_statistics_dict[target]['first_insp_ml_class3'],
                  'second-insp: ', removal_statistics_dict[target]['second_insp_ml_class3'],
                  'red-star: ', removal_statistics_dict[target]['red_star_ml_class3'],
                  'diffraction: ', removal_statistics_dict[target]['diff_ml_class3'])

    def get_data_release_table_name(self, target, classify, cl_class, table_type='obs'):
        """
        Function to create final catalog name
        ----------
        target : str
            name of PHANGS-HST target. Must be in self.phangs_hst_target_list
        classify : str
            classification `human` or `ml`
        cl_class : str
            class group specification either `class12` or `class3`
        table_type : str
            to distinguish between obs and sed table
        Returns
        -------
        table_name : str
            table name
        """

        cat_str = 'hlsp_phangs-cat_hst_'

        instruments = ''
        if self.phangs_hst_obs_band_dict[target]['acs_wfc1_observed_bands']:
            instruments += 'acs'
            if self.phangs_hst_obs_band_dict[target]['wfc3_uvis_observed_bands']:
                instruments += '-uvis'
        else:
            instruments += 'uvis'

        cat_str += instruments + '_'
        cat_str += target.lower() + '_'
        cat_str += 'multi' + '_'
        cat_str += 'v1' + '_'
        if table_type == 'obs':
            cat_str += 'obs' + '-'
        if table_type == 'sed':
            cat_str += 'sed' + '-'
        if table_type == 'cand':
            cat_str += 'obs-sed' + '-'
        if classify == 'human':
            cat_str += 'human' + '-'
        if classify == 'ml':
            cat_str += 'machine' + '-'
        if cl_class == 'class12':
            cat_str += 'cluster-class12'
        if cl_class == 'class3':
            cat_str += 'compact-association-class3'
        if cl_class is None:
            cat_str += 'candidates'
        cat_str += '.fits'

        return cat_str

        # # old version
        # if table_type == 'obs':
        #     return ('PHANGS_HST_cluster_catalog_dr_%s_cat_release_%s_%s_obs_%s_%s.fits' %
        #             (self.data_release, self.catalog_release, target, classify, cl_class))
        # elif table_type == 'sed':
        #     return ('PHANGS_HST_cluster_catalog_dr_%s_cat_release_%s_%s_sed_%s_%s.fits' %
        #             (self.data_release, self.catalog_release, target, classify, cl_class))
        # elif table_type == 'candidates':
        #     return ('PHANGS_HST_candidate_catalog_dr_%s_cat_release_%s_%s_obs_sed.fits' %
        #             (self.data_release, self.catalog_release, target))

    def create_obs_table(self, target, table, cand_table, classify):
        """
        Function to convert an internal data release table into a final data release table
        ----------
        target : str
            target for which the IR catalog should be identified. Must be in self.phangs_hst_target_list
        table : type ``astropy.io.fits.fitsrec.FITS_rec``
            input fits table
        cand_table : type ``astropy.io.fits.fitsrec.FITS_rec``
            candiate table for index crossmatch fits table
        classify : str
            classification human or ml
        Returns
        -------
        table_ir : ``astropy.table.Table``
            Final data release table for one object
        """
        # get column names for the catalog
        column_name_list = self.get_obs_table_column_list(target)
        # create table
        obs_table = None
        for col_name in column_name_list:
            column_content = table[col_name]
            # flux correction
            if target in ['ngc1512', 'ngc1510']:
                if (self.cat_info[col_name]['col_name'][-3:] == 'mJy') | (self.cat_info[col_name]['col_name'][-7:] == 'mJy_ERR'):
                    # now check which entrances have a detection
                    mask_content = column_content != -9999.0
                    column_content[mask_content] *= self.ngc1512_app_corr_offset_flux
                if self.cat_info[col_name]['col_name'][-4:] == 'VEGA':
                    mask_content = column_content != -9999.0
                    column_content[mask_content] += self.ngc1512_app_corr_offset_mag

            if self.cat_info[col_name]['unit'] is not None:
                column_content *= self.cat_info[col_name]['unit']
            column = Column(data=column_content,
                            name=self.cat_info[col_name]['col_name'],
                            dtype=column_content.dtype,
                            description=self.cat_info[col_name]['doc_comment'])

            if obs_table is None:
                obs_table = column
            else:
                obs_table = hstack([obs_table, column])
        # now add cluster and all source id
        data_cluster_id = np.ones(len(obs_table), dtype=int) * -999
        data_allsource_id = np.ones(len(obs_table), dtype=int) * -999
        for obj_index in range(len(obs_table)):
            obj_mask = ((cand_table['PHANGS_X'] == obs_table['PHANGS_X'][obj_index]) &
                        (cand_table['PHANGS_Y'] == obs_table['PHANGS_Y'][obj_index]))
            data_cluster_id[obj_index] = cand_table['ID_PHANGS_CLUSTER'][obj_mask]
            data_allsource_id[obj_index] = cand_table['ID_PHANGS_ALLSOURCES'][obj_mask]

        cluster_id_column = Column(data=data_cluster_id,
                                   name=self.cat_info['id_phangs_cluster']['col_name'],
                                   dtype=data_cluster_id.dtype,
                                   description=self.cat_info['id_phangs_cluster']['doc_comment'])
        allsource_id_column = Column(data=data_allsource_id,
                                     name=self.cat_info['ID_PHANGS_ALLSOURCES_v1p2']['col_name'],
                                     dtype=data_allsource_id.dtype,
                                     description=self.cat_info['ID_PHANGS_ALLSOURCES_v1p2']['doc_comment'])

        obs_table.add_column(cluster_id_column, index=1)
        obs_table.add_column(allsource_id_column, index=3)

        # add region column:
        young_mask = table['SEDfix_age'] < 10
        vi_color = table['PHANGS_F555W_VEGA_TOT'] - table['PHANGS_F814W_VEGA_TOT']
        color_u = table['PHANGS_F336W_VEGA_TOT']
        if 'F438W' in self.phangs_hst_obs_band_dict[target]['wfc3_uvis_observed_bands']:
            color_b = table['PHANGS_F438W_VEGA_TOT']
        else:
            color_b = table['PHANGS_F435W_VEGA_TOT']
        ub_color = color_u - color_b

        vi_hull_ycl, ub_hull_ycl = self.load_hull(region_type='ycl', classify=classify, y_color='ub')
        vi_hull_map, ub_hull_map = self.load_hull(region_type='map', classify=classify, y_color='ub')
        vi_hull_ogcc, ub_hull_ogcc = self.load_hull(region_type='ogcc', classify=classify, y_color='ub')

        hull_ycl = ConvexHull(np.array([vi_hull_ycl, ub_hull_ycl]).T)
        hull_map = ConvexHull(np.array([vi_hull_map, ub_hull_map]).T)
        hull_ogcc = ConvexHull(np.array([vi_hull_ogcc, ub_hull_ogcc]).T)

        in_hull_ycl = helper_func.points_in_hull(np.array([vi_color, ub_color]).T, hull_ycl)
        in_hull_map = helper_func.points_in_hull(np.array([vi_color, ub_color]).T, hull_map)
        in_hull_ogcc = helper_func.points_in_hull(np.array([vi_color, ub_color]).T, hull_ogcc)

        mask_ycl = in_hull_ycl * np.invert(in_hull_map) + young_mask * (in_hull_map * in_hull_ogcc)
        mask_map = in_hull_map * np.invert(young_mask)
        mask_ogcc = in_hull_ogcc * np.invert(young_mask)

        column_data_region = np.array(['outside'] * len(table['SEDfix_age']), dtype=str)

        column_data_region[mask_ycl] = 'YCL'
        column_data_region[mask_map] = 'MAP'
        column_data_region[mask_ogcc] = 'OGCC'

        region_column = Column(data=column_data_region,
                               name=self.cat_info['cc_class']['col_name'],
                               dtype=str,
                               description=self.cat_info['cc_class']['doc_comment'])
        obs_table.add_column(region_column, index=-1)

        return obs_table

    def create_sed_table(self, target, table, cand_table):
        """
        Function to convert an internal data release table into a final data release table
        ----------
        target : str
            target for which the IR catalog should be identified. Must be in self.phangs_hst_target_list
        table : type ``astropy.io.fits.fitsrec.FITS_rec``
            input fits table
        cand_table : type ``astropy.io.fits.fitsrec.FITS_rec``
            candiate table for index crossmatch fits table
        Returns
        -------
        table_ir : ``astropy.table.Table``
            Final data release table for one object
        """
        # get column names for the catalog
        column_name_list = self.tab2_columns
        # create table
        sed_table = None
        for col_name in column_name_list:
            column_content = table[col_name]
            # mass correction
            if target in ['ngc1512', 'ngc1510']:
                if (self.cat_info[col_name]['col_name'] in ['PHANGS_MASS_MINCHISQ', 'PHANGS_MASS_MINCHISQ_ERR']) | ('mass' in self.cat_info[col_name]['col_name']):
                    # now check which entrances have a detection
                    print(self.cat_info[col_name]['col_name'])
                    mask_content = (column_content != -9999.0) & (column_content != -999.0)
                    column_content[mask_content] *= self.ngc1512_app_corr_offset_flux

            if self.cat_info[col_name]['unit'] is not None:
                column_content *= self.cat_info[col_name]['unit']
            column = Column(data=column_content,
                            name=self.cat_info[col_name]['col_name'],
                            dtype=column_content.dtype,
                            description=self.cat_info[col_name]['doc_comment']
                            )

            if sed_table is None:
                sed_table = column
            else:
                sed_table = hstack([sed_table, column])

        # now add cluster and all source id
        data_cluster_id = np.ones(len(sed_table), dtype=int) * -999
        data_allsource_id = np.ones(len(sed_table), dtype=int) * -999
        for obj_index in range(len(sed_table)):
            obj_mask = ((cand_table['PHANGS_X'] == sed_table['PHANGS_X'][obj_index]) &
                        (cand_table['PHANGS_Y'] == sed_table['PHANGS_Y'][obj_index]))
            data_cluster_id[obj_index] = cand_table['ID_PHANGS_CLUSTER'][obj_mask]
            data_allsource_id[obj_index] = cand_table['ID_PHANGS_ALLSOURCES'][obj_mask]

        cluster_id_column = Column(data=data_cluster_id,
                                   name=self.cat_info['id_phangs_cluster']['col_name'],
                                   dtype=data_cluster_id.dtype,
                                   description=self.cat_info['id_phangs_cluster']['doc_comment'])
        allsource_id_column = Column(data=data_allsource_id,
                                     name=self.cat_info['ID_PHANGS_ALLSOURCES_v1p2']['col_name'],
                                     dtype=data_allsource_id.dtype,
                                     description=self.cat_info['ID_PHANGS_ALLSOURCES_v1p2']['doc_comment'])

        sed_table.add_column(cluster_id_column, index=1)
        sed_table.add_column(allsource_id_column, index=3)

        return sed_table

    def create_cand_table(self, target, table):
        """
        Function to convert an internal data release table into a final data release table
        ----------
        target : str
            target for which the IR catalog should be identified. Must be in self.phangs_hst_target_list
        table : type ``astropy.io.fits.fitsrec.FITS_rec``
            input fits table

        Returns
        -------
        table_ir : ``astropy.table.Table``
            Final data release table for one object
        """
        # get column names for the catalog
        column_name_list = self.get_cand_table_column_list(target)

        # create table
        cand_table = None
        for col_name in column_name_list:
            column_content = table[col_name]
            # flux correction
            if target in ['ngc1512', 'ngc1510']:
                if (self.cat_info[col_name]['col_name'][-3:] == 'mJy') | (self.cat_info[col_name]['col_name'][-7:] == 'mJy_ERR'):
                    # now check which entrances have a detection
                    mask_content = column_content != -9999.0
                    column_content[mask_content] *= self.ngc1512_app_corr_offset_flux
                if self.cat_info[col_name]['col_name'][-4:] == 'VEGA':
                    mask_content = column_content != -9999.0
                    column_content[mask_content] += self.ngc1512_app_corr_offset_mag

                # mass correction
                if (self.cat_info[col_name]['col_name'] in ['PHANGS_MASS_MINCHISQ', 'PHANGS_MASS_MINCHISQ_ERR']) | ('mass' in self.cat_info[col_name]['col_name']):
                    # now check which entrances have a detection
                    print(self.cat_info[col_name]['col_name'])
                    mask_content = (column_content != -9999.0) & (column_content != -999.0)
                    column_content[mask_content] *= self.ngc1512_app_corr_offset_flux

            if self.cat_info[col_name]['unit'] is not None:
                column_content *= self.cat_info[col_name]['unit']
            column = Column(data=column_content,
                            name=self.cat_info[col_name]['col_name'],
                            dtype=column_content.dtype,
                            description=self.cat_info[col_name]['doc_comment']
                            )

            if cand_table is None:
                cand_table = column
            else:
                cand_table = hstack([cand_table, column])

        return cand_table

    def get_ir_cat(self, target, classify, cl_class):
        """
        Function to get internal release catalog identifier
        ----------
        target : str
            target for which the IR catalog should be identified. Must be in self.phangs_hst_target_list
        classify : str
            cluster classification either `human` or `ml`
        cl_class : str
            cluster classes must be either `class12` or `class3`

        Returns
        -------
        table_ir : ``astropy.io.fits.fitsrec.FITS_rec``
            internal release table
        """

        # check if version is supported and get file name
        if self.hst_cc_ver == 'SEDfix_final_test_catalogs':
            # check if 0 id after ngc. this naming convention was not applied to the IR catalog naming files
            if (target[:3] == 'ngc') & (target[3] == '0'):
                target_str = target[0:3] + target[4:]
            else:
                target_str = target
            if cl_class == 'candidates':
                cat_file_name_ir = (
                    Path('SEDfix_%s_Ha1_inclusiveGCcc_inclusiveGCclass_phangshst_candidates_bcw_v1p2_IR4.fits' %
                         target_str))
            else:
                cat_file_name_ir = (
                    Path('SEDfix_PHANGS_IR4_%s_Ha1_inclusiveGCcc_inclusiveGCclass_phangs_hst_v1p2_%s_%s.fits'
                         % (target_str, classify, cl_class)))
        else:
            raise AttributeError(' the specified attribute hst_cc_ver is ', self.hst_cc_ver,
                                 ' Which is not supported by this version.')

        # get file path
        file_path_ir = Path(self.path2ir) / cat_file_name_ir
        # check if file exists
        if not os.path.isfile(file_path_ir):
            print(file_path_ir, ' not found ! ')
            raise FileNotFoundError('there is no HST cluster catalog for the target ', target,
                                    ' make sure that the file ', file_path_ir, ' exists.')
        # open table and get column names
        table_hdu = fits.open(file_path_ir)
        table_ir = table_hdu[1].data
        table_hdu.close()

        return table_ir

    def get_artifact_cat(self, target):
        """
        Function to get artifact catalog name
        ----------
        target : str
            target for which the IR catalog should be identified. Must be in ``self.phangs_hst_target_list``

        Returns
        -------
        table_artifact : ``astropy.io.fits.fitsrec.FITS_rec`` or None
            artifact table
        """
        # make sure the zero is removed after NGC for some objects
        if (target[:3] == 'ngc') & (target[3] == '0'):
            target_str = target[:3] + target[4:]
        else:
            target_str = target

        # get artifact table file path
        file_path_artifact = self.path2artifact / Path('%s_artifacts.fits' % target_str)
        # check if file exists
        if (not os.path.isfile(file_path_artifact)) & self.artifact_rais_file_not_found_flag:
            print(file_path_artifact, ' not found ! ')
            raise FileNotFoundError('there is no artigfact table for the target ', target,
                                    ' make sure that the file ', file_path_artifact, ' exists.')
        elif not os.path.isfile(file_path_artifact):
            print('No artifact table found for ', target)
            return None
        else:

            # open table and get column names
            table_hdu = fits.open(file_path_artifact)
            table_artifact = table_hdu[1].data
            table_hdu.close()

            return table_artifact

    def load_second_classify_table(self):
        """
        Function to load the secondary classification table into the attributes
        ----------
        ...
        Returns
        -------
        None
        """
        self.questionable_artefact_table = ascii.read(self.path2questionable_artifacts)

    def get_second_classify_cat(self, target):
        """
        Function to get artifact catalog that was reclassified because the first artefact classification was
        questionable.
        ----------
        target : str
            target for which the IR catalog should be identified. Must be in ``self.phangs_hst_target_list``

        Returns
        -------
        table_artifact : ``astropy.io.fits.fitsrec.FITS_rec`` or None
            artifact table
        """
        # make sure the zero is removed after NGC for some objects

        if self.questionable_artefact_table is None:
            self.load_second_classify_table()

        if (target[:3] == 'ngc') & (target[3] == '0'):
            target_str = target[:3] + target[4:]
        else:
            target_str = target

        mask_target = self.questionable_artefact_table['GALAXY'] == target_str

        return self.questionable_artefact_table[mask_target]

    def load_diffraction_spike_masks(self, target, target_str):
        """
        Function to get masks for diffraction spikes
        ----------
        target : str
            target for which the IR catalog should be identified. Must be in ``self.phangs_hst_target_list``

        Returns
        -------
        diffraction_spike_mask : ndarray
        """

        if (target_str[:3] == 'ngc') & (target_str[3] == '0'):
            target_str = target[:3] + target[4:]

        diffraction_mask = None
        diffraction_wcs = None
        for band in self.phangs_hst_obs_band_dict[target]['wfc3_uvis_observed_bands']:
            diff_spike_mask_file_name = (self.path2diffraction_spike_masks + '/' '%s_%s_uvis_mask.fits' %
                                         (target_str, band.lower()))
            if os.path.isfile(diff_spike_mask_file_name):
                hdu = fits.open(diff_spike_mask_file_name)
                if diffraction_mask is None:
                    diffraction_mask = hdu[0].data
                    diffraction_wcs = WCS(hdu[0].header)
                else:
                    diffraction_mask += hdu[0].data

        for band in self.phangs_hst_obs_band_dict[target]['acs_wfc1_observed_bands']:
            diff_spike_mask_file_name = (self.path2diffraction_spike_masks + '/' '%s_%s_acs_mask.fits' %
                                         (target_str, band.lower()))
            if os.path.isfile(diff_spike_mask_file_name):
                hdu = fits.open(diff_spike_mask_file_name)
                if diffraction_mask is None:
                    diffraction_mask = hdu[0].data
                    diffraction_wcs = WCS(hdu[0].header)

                else:
                    diffraction_mask += hdu[0].data

        return diffraction_mask, diffraction_wcs

    def write_string_to_file(self, file, string, max_length):
        if len(string) < max_length:
            file.writelines(string + ' \n')
        else:
            # in order to shorten the description line we need to split the description string
            shortened_str, current_str = self.split_line(line_string=string, max_length=max_length)
            # write first line of the string
            file.writelines(shortened_str + ' \n')

            # use while loop to keep on dividing the rest of the description string till it is short enough
            while_flag = True
            while while_flag:
                if len(current_str) > max_length:
                    shortened_str, current_str = self.split_line(line_string=current_str,
                                                                 max_length=max_length)
                else:
                    shortened_str = current_str
                    # after this condition is fulfilled the while loop will stop after the current iteration
                    while_flag = False
                file.writelines(shortened_str + ' \n')

    def create_documentation(self, txt_col_width=80):

        # get the documentation_file_name
        # doc_name = 'PHANGS_DR%s_CATR%s_compact_clusters_README.txt' % (self.data_release, self.catalog_release)
        doc_name = 'hlsp_phangs-cat_hst_multi_all_multi_v1_readme.txt'
        doc_file = open(self.catalog_output_path + '/' + doc_name, "w")

        # divider for different sections
        section_divider = '-' * txt_col_width
        header_line = '-' * int(txt_col_width/2)

        # start filling documentation
        doc_file.writelines(section_divider + ' \n')

        # write header
        header_str_intro = ('This is the documentation to describe the data release of PHANGS-HST star clusters and '
                            'compact associations for 39 galaxies.')
        split_header_line = self.split_line(line_string=header_str_intro, max_length=txt_col_width)
        header_str_ver = 'Data Release %s, Catalog Release %s' % (self.data_release, self.catalog_release)
        header_str_date = datetime.now().strftime("%d %B, %Y")
        for line in split_header_line:
            doc_file.writelines(line + ' \n')
        doc_file.writelines(header_str_ver + ' \n')
        doc_file.writelines(header_str_date + ' \n')
        doc_file.writelines(section_divider + ' \n')
        doc_file.writelines(' \n')

        # adding descriptive text
        descriptive_text_intro = \
            ('This README provides details pertaining to the content of all data files included in this data release. '
             'Further information is provided in the article Maschmann et al.(submitted). ')

        descriptive_text_table_types_0 = \
            ('Here, we provide a set of 9 catalogs for each of the 38 galaxies. '
             'We note that the galaxy NGC628 is devided into '
             'two pointings NGC628C and NGC628E which are treated here as individual targets. ')

        table_naming_list = \
            ['For each target the following 9 catalogs are provided:',
             # 'Candidate table'
             ]
        test_target = 'ngc7496'
        # candidate_name_str = self.get_data_release_table_name(target=test_target, classify=None,
        #                                                               cl_class=None, table_type='cand')
        # candidate_name_str = (
        #     candidate_name_str.replace('uvis', '<ins>').replace(test_target, '<target>'))
        # table_naming_list.append(candidate_name_str)
        table_naming_list.append('')
        table_naming_list.append('Human classified tables')
        table_naming_list.append(
            self.get_data_release_table_name(
                target=test_target, classify='human', cl_class='class12', table_type='obs')
            .replace('uvis', '<ins>').replace(test_target, '<target>'))
        table_naming_list.append(
            self.get_data_release_table_name(
                target=test_target, classify='human', cl_class='class3', table_type='obs')
            .replace('uvis', '<ins>').replace(test_target, '<target>'))
        # table_naming_list.append(
        #     self.get_data_release_table_name(
        #         target=test_target, classify='human', cl_class='class12', table_type='sed')
        #     .replace('uvis', '<ins>').replace(test_target, '<target>'))
        # table_naming_list.append(
        #     self.get_data_release_table_name(
        #         target=test_target, classify='human', cl_class='class3', table_type='sed')
        #     .replace('uvis', '<ins>').replace(test_target, '<target>'))
        table_naming_list.append('')
        table_naming_list.append('Machine learning classified tables')
        table_naming_list.append(
            self.get_data_release_table_name(
                target=test_target, classify='ml', cl_class='class12', table_type='obs')
            .replace('uvis', '<ins>').replace(test_target, '<target>'))
        table_naming_list.append(
            self.get_data_release_table_name(
                target=test_target, classify='ml', cl_class='class3', table_type='obs')
            .replace('uvis', '<ins>').replace(test_target, '<target>'))
        # table_naming_list.append(
        #     self.get_data_release_table_name(
        #         target=test_target, classify='ml', cl_class='class12', table_type='sed')
        #     .replace('uvis', '<ins>').replace(test_target, '<target>'))
        # table_naming_list.append(
        #     self.get_data_release_table_name(
        #         target=test_target, classify='ml', cl_class='class3', table_type='sed')
        #     .replace('uvis', '<ins>').replace(test_target, '<target>'))
        table_naming_list.append('')

        table_names_explanation_list = \
            ('<ins> denotes the instrument used for the observation. It can be can be  WFC and/or UVIS. '
             'This is specified in Table 1 in Maschmann et al.(submitted)).',
             ' ',
             'Catalogs with the keyword \"obs\" contain observational properties of the objects. For more details see '
             'Maschmann et al.(submitted)).',
             # ' ',
             # 'Catalogs with the keyword \"sed\" contain the results of the SED fitting for more details see '
             # 'Thilker et al.(submitted)).',
             ' ',
             'keywords \"human\" and \"machine\" refer to the method used to classify the objects.',
             ' ',
             'The last keyword is either \"cluster_class12\" or \"compact_association_class3\" '
             'and denotes the object’s class',
             ' ',
             ' These keywords for the classification method and object class are not used in the name of the '
             'candidate table since this table contains all candidates before the classification step',
             )
        # write the intro lines:
        self.write_string_to_file(file=doc_file, string=descriptive_text_intro, max_length=txt_col_width)
        doc_file.writelines(' \n')
        self.write_string_to_file(file=doc_file, string=descriptive_text_table_types_0, max_length=txt_col_width)
        doc_file.writelines(' \n')
        for line in table_naming_list:
            doc_file.writelines(line + ' \n')
        for line in table_names_explanation_list:
            if len(line) > txt_col_width:
                self.write_string_to_file(file=doc_file, string=line, max_length=txt_col_width)
            else:
                doc_file.writelines(line + ' \n')

        doc_file.writelines(section_divider + ' \n')
        doc_file.writelines(' \n')

        # print column names description
        running_column_index = 1
        doc_file.writelines(header_line + ' \n')
        doc_file.writelines('Object Identifiers \n')
        doc_file.writelines(header_line + ' \n')

        for col_name in self.doc_identifier_names:
            print(col_name)
            print(self.cat_info[col_name]['col_name'])
            doc_file.writelines('[%i] ' % running_column_index +
                                self.cat_info[col_name]['col_name'] + ' [' +
                                self.cat_info[col_name]['unit_str'] + ']' + ' \n')
            self.write_string_to_file(file=doc_file, string=self.cat_info[col_name]['doc_comment'], max_length=txt_col_width)
            doc_file.writelines(' \n')
            running_column_index += 1

        doc_file.writelines(header_line + ' \n')
        doc_file.writelines('Coordinates \n')
        doc_file.writelines(header_line + ' \n')
        for col_name in self.doc_coord_names:
            print(col_name)
            print(self.cat_info[col_name]['col_name'])
            doc_file.writelines('[%i] ' % running_column_index +
                                self.cat_info[col_name]['col_name'] + ' [' +
                                self.cat_info[col_name]['unit_str'] + ']' + ' \n')
            self.write_string_to_file(file=doc_file, string=self.cat_info[col_name]['doc_comment'], max_length=txt_col_width)
            doc_file.writelines(' \n')
            running_column_index += 1

        doc_file.writelines(header_line + ' \n')
        doc_file.writelines('Morphological Classification \n')
        doc_file.writelines(header_line + ' \n')

        # explanation for morphological classification
        morph_explanation_str = \
            ('We use the following numbers to distinguish different object types: ',
             ' 1	Star cluster - single peak, circularly symmetric',
             ' 2	Star cluster - single peak, but elongated or asymmetric',
             ' 3	Compact stellar association – asymmetric, multiple peaks',
             ' 4	Artifact: star',
             ' 5	Artifact: pair of stars',
             ' 6	Artifact: saturated star',
             ' 7	Artifact: redundant - other object within 5 pixels',
             ' 8	Artifact: diffraction spike and bleeding',
             ' 9	Artifact: background galaxy',
             ' 10	Artifact: too faint to tell',
             ' 11	Artifact: clean stars',
             ' 12	Artifact: edge',
             ' 13	Triple',
             ' 14	Bad pixel',
             ' 15	Nucleus of galaxy',
             ' 16	Fluff (i.e. no clear peaks)',
             ' 17	Uncategorized artifact',
             ' 19	Very red star with V-I color > 2.0 and CI < 1.45')

        for line in morph_explanation_str:
            doc_file.writelines(line + ' \n')
        doc_file.writelines(header_line + ' \n')

        for col_name in self.classification_columns_doc:
            print(col_name)
            print(self.cat_info[col_name]['col_name'])
            doc_file.writelines('[%i] ' % running_column_index +
                                self.cat_info[col_name]['col_name'] + ' [' +
                                self.cat_info[col_name]['unit_str'] + ']' + ' \n')
            self.write_string_to_file(file=doc_file, string=self.cat_info[col_name]['doc_comment'],
                                      max_length=txt_col_width)
            doc_file.writelines(' \n')
            running_column_index += 1

        doc_file.writelines(header_line + ' \n')
        doc_file.writelines('Photometry \n')
        doc_file.writelines(header_line + ' \n')
        for col_name in self.doc_photometry_columns:
            print(col_name)
            print(self.cat_info[col_name]['col_name'])
            doc_file.writelines('[%i] ' % running_column_index +
                                self.cat_info[col_name]['col_name'] + ' [' +
                                self.cat_info[col_name]['unit_str'] + ']' + ' \n')
            self.write_string_to_file(file=doc_file, string=self.cat_info[col_name]['doc_comment'],
                                      max_length=txt_col_width)
            doc_file.writelines(' \n')
            running_column_index += 1
        for col_name in self.detect_shape_columns:
            print(col_name)
            print(self.cat_info[col_name]['col_name'])
            doc_file.writelines('[%i] ' % running_column_index +
                                self.cat_info[col_name]['col_name'] + ' [' +
                                self.cat_info[col_name]['unit_str'] + ']' + ' \n')
            self.write_string_to_file(file=doc_file, string=self.cat_info[col_name]['doc_comment'], max_length=txt_col_width)
            doc_file.writelines(' \n')
            running_column_index += 1


        # doc_file.writelines(header_line + ' \n')
        # doc_file.writelines('SED fitting results \n')
        # doc_file.writelines(header_line + ' \n')
        # for col_name in self.sed_columns_doc:
        #     print(col_name)
        #     print(self.cat_info[col_name]['col_name'])
        #     doc_file.writelines('[%i] ' % running_column_index +
        #                         self.cat_info[col_name]['col_name'] + ' [' +
        #                         self.cat_info[col_name]['unit_str'] + ']' + ' \n')
        #     self.write_string_to_file(file=doc_file, string=self.cat_info[col_name]['doc_comment'], max_length=txt_col_width)
        #     doc_file.writelines(' \n')
        #     running_column_index += 1


        # now print galaxy lists
        doc_file.writelines(section_divider + ' \n')
        doc_file.writelines(' \n')
        doc_file.writelines(header_line + ' \n')
        doc_file.writelines('Color-color Hulls \n')
        doc_file.writelines(header_line + ' \n')

        # write explanation for the distance table
        description_str = ('We further provide the hulls describing the regions found in the color-color diagrams '
                           'in Figures 10 and 11 of Maschmann et al. (submitted). '
                           'These regions are computed for U-B vs V-I and NUV-B vs V-I diagrams. '
                           'We identified three regions: Young Cluster Locus (YCL), Middle Aged Plume (MAP) and '
                           'Old Globular Cluster Clump (OGCC). We computed these regions for the human and ML sample '
                           'separately. '
                           'In the names of the region files we specify the region (ycl, map or ogcc), which sample '
                           'was used to compute them (human or ml) and which color-color diagram type was used.')
        self.write_string_to_file(file=doc_file, string=description_str, max_length=txt_col_width)
        doc_file.writelines(' \n')

        table_naming_list = \
            ('The following region files are available: ',
             'hlsp_phangs-cat_hst_multi_hull_multi_v1_ycl-human-nuvbvi',
             'hlsp_phangs-cat_hst_multi_hull_multi_v1_map-human-nuvbvi',
             'hlsp_phangs-cat_hst_multi_hull_multi_v1_ogcc-human-nuvbvi',
             'hlsp_phangs-cat_hst_multi_hull_multi_v1_ycl-ml-nuvbvi',
             'hlsp_phangs-cat_hst_multi_hull_multi_v1_map-ml-nuvbvi',
             'hlsp_phangs-cat_hst_multi_hull_multi_v1_ogcc-ml-nuvbvi',
             'hlsp_phangs-cat_hst_multi_hull_multi_v1_ycl-human-ubvi',
             'hlsp_phangs-cat_hst_multi_hull_multi_v1_map-human-ubvi',
             'hlsp_phangs-cat_hst_multi_hull_multi_v1_ogcc-human-ubvi',
             'hlsp_phangs-cat_hst_multi_hull_multi_v1_ycl-ml-ubvi',
             'hlsp_phangs-cat_hst_multi_hull_multi_v1_map-ml-ubvi',
             'hlsp_phangs-cat_hst_multi_hull_multi_v1_ogcc-ml-ubvi',
             '',)
        for line in table_naming_list:
            doc_file.writelines(line + ' \n')



        # now print galaxy lists
        doc_file.writelines(section_divider + ' \n')
        doc_file.writelines(' \n')
        doc_file.writelines(header_line + ' \n')
        doc_file.writelines('Target Distances \n')
        doc_file.writelines(header_line + ' \n')

        # write explanation for the distance table
        description_str = ('The distances are listed for each PHANGS-HST galaxy. The columns are: (1) target name, '
                           '(2) distance in Mpc, (3) distance uncertainty in Mpc, (4) method used for the distance '
                           'measurement and (5) reference. Column 4 and 5 are described in more detail in '
                           'Lee et al. (2022) (2022ApJS..258...10L) see Table 1. '
                           'For a detailed description of the distance measurements '
                           'see Anand et al. (2021) (2021MNRAS.501.3621A)')
        self.write_string_to_file(file=doc_file, string=description_str, max_length=txt_col_width)
        doc_file.writelines(' \n')

        path_sample_table = ('/home/benutzer/data/PHANGS_products/sample_table/v1p6' + '/phangs_sample_table_v1p6.fits')
        hdu_sample_table = fits.open(path_sample_table)
        data_sample_table = hdu_sample_table[1].data
        for target in self.phangs_galaxy_list:
            mask_target = data_sample_table['name'] == target

            dist = data_sample_table['dist'][mask_target]
            dist_unc_dex = data_sample_table['dist_unc'][mask_target]

            dist_unc = 10 ** (np.log10(dist) + dist_unc_dex) - dist

            dist_label = data_sample_table['dist_label'][mask_target]
            dist_ref = data_sample_table['dist_ref'][mask_target]
            if dist_ref == 'PHANGSTRGB':
                dist_ref = '2021MNRAS.501.3621A'
            elif dist_ref == 'PHANGSPNLF':
                dist_ref = '2022MNRAS.511.6087S'
            else:
                dist_ref = dist_ref[0]

            if target == 'ngc1512':
                doc_file.writelines('NGC1510 %.2f %.2f %s %s \n' % (dist, dist_unc, dist_label[0], dist_ref))

            doc_file.writelines('%s %.2f %.2f %s %s \n' % (target.upper(), dist, dist_unc, dist_label[0], dist_ref))

        doc_file.writelines(' \n')


        doc_file.writelines(section_divider + ' \n')
        doc_file.writelines(' \n')
        doc_file.writelines(header_line + ' \n')
        doc_file.writelines('Photometric Corrections \n')
        doc_file.writelines(header_line + ' \n')

        # write explanation for the distance table
        description_str = ('Photometry has been corrected for foreground Milky Way extinction. '
                           'Values for A(V) toward each target were first taken from '
                           'Schlafly and Finkbeiner (2011) (2011ApJ...737..103S) using the IRSA-hosted web service, '
                           'and then specifically scaled to our HST bandpasses using a Fitzpatrick (1999) '
                           '(1999PASP..111...63F) extinction curve with R_V = 3.1. We also provide this information as '
                           'a standalone file in the main directory.')
        self.write_string_to_file(file=doc_file, string=description_str, max_length=txt_col_width)
        doc_file.writelines(' \n')

        mw_corr_list = ('galaxy A_NUV A_U A_B A_V A_I',
                        'ic1954 0.0886 0.0718 0.0586 0.0453 0.0241',
                        'ic5332 0.0906 0.0734 0.0599 0.0463 0.0247',
                        'ngc0628c 0.379 0.307 0.25 0.191 0.103',
                        'ngc0628e 0.379 0.307 0.25 0.1940 0.103',
                        'ngc0685 0.125 0.1010 0.0825 0.0638 0.034',
                        'ngc1087 0.187 0.152 0.124 0.0957 0.0509',
                        'ngc1097 0.145 0.1170 0.0958 0.0741 0.0394',
                        'ngc1300 0.163 0.132 0.108 0.0821 0.0444',
                        'ngc1317 0.113 0.0914 0.0745 0.0576 0.0307',
                        'ngc1365 0.111 0.0897 0.0732 0.0566 0.0301',
                        'ngc1385 0.109 0.0881 0.0719 0.0556 0.0296',
                        'ngc1433 0.0483 0.0392 0.0319 0.0247 0.0131',
                        'ngc1510 0.0564 0.0457 0.0373 0.0288 0.0153',
                        'ngc1512 0.0564 0.0457 0.0373 0.0288 0.0153',
                        'ngc1559 0.161 0.131 0.106 0.0824 0.0438',
                        'ngc1566 0.0483 0.0392 0.0319 0.0247 0.0131',
                        'ngc1672 0.125 0.1010 0.0824 0.0638 0.034',
                        'ngc1792 0.121 0.0979 0.0798 0.0618 0.0329',
                        'ngc2775 0.2340 0.1890 0.154 0.119 0.0635',
                        'ngc2835 0.54 0.437 0.357 0.276 0.147',
                        'ngc2903 0.1670 0.135 0.11 0.0854 0.0455',
                        'ngc3351 0.149 0.121 0.0985 0.0762 0.0405',
                        'ngc3621 0.437 0.354 0.289 0.22 0.119',
                        'ngc3627 0.181 0.147 0.12 0.0926 0.0493',
                        'ngc4254 0.209 0.17 0.138 0.107 0.057',
                        'ngc4298 0.191 0.155 0.126 0.0978 0.052',
                        'ngc4303 0.121 0.0979 0.0798 0.0618 0.0329',
                        'ngc4321 0.141 0.114 0.0932 0.0721 0.0383',
                        'ngc4535 0.105 0.0848 0.0692 0.0535 0.0285',
                        'ngc4536 0.0987 0.0799 0.0652 0.0504 0.0268',
                        'ngc4548 0.205 0.166 0.136 0.105 0.0559',
                        'ngc4569 0.254 0.206 0.168 0.13 0.069',
                        'ngc4571 0.25 0.2020 0.165 0.128 0.0679',
                        'ngc4654 0.141 0.114 0.0932 0.0721 0.0383',
                        'ngc4689 0.123 0.0995 0.0812 0.0628 0.0334',
                        'ngc4826 0.2240 0.181 0.148 0.114 0.0608',
                        'ngc5068 0.5540 0.449 0.366 0.283 0.151',
                        'ngc5248 0.133 0.108 0.0878 0.0679 0.0362',
                        'ngc6744 0.2320 0.188 0.153 0.118 0.063',
                        'ngc7496 0.0524 0.0424 0.0346 0.0268 0.0142')
        for mw_corr in mw_corr_list:
            doc_file.writelines(mw_corr + ' \n')
        doc_file.writelines(' \n')

        doc_file.writelines(section_divider + ' \n')
        doc_file.writelines(' \n')
        doc_file.writelines(header_line + ' \n')
        doc_file.writelines('Aperture Corrections \n')
        doc_file.writelines(header_line + ' \n')

        # write explanation for the distance table
        description_str_1 = ('Aperture corrections have been applied to the photometry provided in the '
                             'compact cluster catalogs. The details of our aperture correction computation methodology '
                             'have been presented in Deger et al. (2022) (2022MNRAS.510...32D) and '
                             'Lee et al. (2022) (2022ApJS..258...10L). The rationale for '
                             'choosing to use a constant offset to derive aperture corrections in bands other than '
                             'the V-band, rather than deriving corrections directly in each filter is given in'
                             ' Deger et al. (2022) (2022MNRAS.510...32D). '
                             'The offsets we apply to measured V-band aperture correction '
                             'for NUV, U, B, and I bands are -0.19, -0.12, -0.03, and -0.12 mags, respectively. '
                             'These corrections are determined on the basis of the FWHM of the WFC3 point spread'
                             ' function at their corresponding wavelengths.')
        description_str_2 = ('The aperture corrections applied to the five band photometry of the galaxies '
                             'published in this release are as follows. The values provided in the table below '
                             'are in Vega magnitudes.')

        self.write_string_to_file(file=doc_file, string=description_str_1, max_length=txt_col_width)
        doc_file.writelines(' \n')
        self.write_string_to_file(file=doc_file, string=description_str_2, max_length=txt_col_width)
        doc_file.writelines(' \n')

        ap_corr_list = ('galaxy NUV_corr U_corr B_corr V_corr I_corr',
                        'ic1954 -0.81 -0.74 -0.65 -0.62 -0.74',
                        'ic5332 -0.78 -0.71 -0.62 -0.59 -0.71',
                        'ngc628c -0.94 -0.87 -0.78 -0.75 -0.87',
                        'ngc628e -0.91 -0.84 -0.75 -0.72 -0.84',
                        'ngc0685 -0.85 -0.78 -0.69 -0.66 -0.78',
                        'ngc1087 -0.86 -0.79 -0.70 -0.67 -0.79',
                        'ngc1097 -0.82 -0.75 -0.66 -0.63 -0.75',
                        'ngc1300 -0.8 -0.73 -0.64 -0.61 -0.73',
                        'ngc1317 -0.8 -0.73 -0.64 -0.61 -0.73',
                        'ngc1365 -0.8 -0.73 -0.64 -0.61 -0.73',
                        'ngc1385 -0.88 -0.81 -0.72 -0.69 -0.81',
                        'ngc1433 -0.81 -0.74 -0.65 -0.62 -0.74',
                        'ngc1510 -2.534 -2.464 -2.374 -2.344 -2.464',
                        'ngc1512 -2.534 -2.464 -2.374 -2.344 -2.464',
                        'ngc1559 -0.86 -0.79 -0.70 -0.67 -0.79',
                        'ngc1566 -0.81 -0.74 -0.65 -0.62 -0.74',
                        'ngc1672 -0.79 -0.72 -0.63 -0.6 -0.72',
                        'ngc1792 -0.99 -0.92 -0.83 -0.8 -0.92',
                        'ngc2775 -0.64 -0.57 -0.48 -0.45 -0.57',
                        'ngc2835 -1.0 -0.93 -0.84 -0.81 -0.93',
                        'ngc2903 -0.89 -0.82 -0.73 -0.7 -0.82',
                        'ngc3351 -0.87 -0.8 -0.71 -0.68 -0.8',
                        'ngc3621 -0.96 -0.89 -0.8 -0.77 -0.89',
                        'ngc3627 -1.0 -0.93 -0.84 -0.81 -0.93',
                        'ngc4254 -0.99 -0.92 -0.83 -0.8 -0.92',
                        'ngc4298 -0.79 -0.72 -0.63 -0.6 -0.72',
                        'ngc4303 -0.93 -0.86 -0.77 -0.74 -0.86',
                        'ngc4321 -0.88 -0.81 -0.72 -0.69 -0.81',
                        'ngc4535 -0.90 -0.83 -0.74 -0.71 -0.83',
                        'ngc4536 -0.82 -0.75 -0.66 -0.63 -0.75',
                        'ngc4548 -0.83 -0.76 -0.67 -0.64 -0.76',
                        'ngc4569 -0.95 -0.88 -0.79 -0.76 -0.88',
                        'ngc4571 -0.8 -0.73 -0.64 -0.61 -0.73',
                        'ngc4654 -1.02 -0.95 -0.86 -0.83 -0.95',
                        'ngc4689 -0.82 -0.75 -0.66 -0.63 -0.75',
                        'ngc4826 -0.91 -0.84 -0.75 -0.72 -0.84',
                        'ngc5068 -0.96 -0.89 -0.8 -0.77 -0.89',
                        'ngc5248 -0.84 -0.77 -0.68 -0.65 -0.77',
                        'ngc6744 -0.73 -0.66 -0.57 -0.54 -0.66',
                        'ngc7496 -0.73 -0.66 -0.57 -0.54 -0.66')
        for ap_corr in ap_corr_list:
            doc_file.writelines(ap_corr + ' \n')
        doc_file.writelines(' \n')

    def load_hull(self, region_type, classify, y_color='ub'):
        if region_type == 'ycl':
            region_str = 'young'
            class_number = 3
        elif region_type == 'map':
            region_str = 'mid'
            class_number = 1
        elif region_type == 'ogcc':
            region_str = 'ogc'
            class_number = 1
        else:
            raise KeyError('region string not understood')

        if classify == 'human':
            classify_str = 'hum'
        elif classify == 'ml':
            classify_str = 'ml'
        else:
            raise KeyError('classify string not understood')

        x_hull = np.load('/home/benutzer/Documents/projects/hst_cluster_catalog/analysis/segmentation/'
                         'data_output/vi_hull_%s_%svi_%s_%i.npy' % (region_str, y_color, classify_str, class_number))
        y_hull = np.load('/home/benutzer/Documents/projects/hst_cluster_catalog/analysis/segmentation/'
                         'data_output/%s_hull_%s_%svi_%s_%i.npy' %
                         (y_color, region_str, y_color, classify_str, class_number))

        return x_hull, y_hull

    def create_hull_files(self):

        for region_type in ['ycl', 'map', 'ogcc']:
            for classify in ['human', 'ml']:
                for y_color in ['nuvb', 'ub']:
                    x_hull, y_hull = self.load_hull(region_type=region_type, classify=classify, y_color=y_color)

                    # get the documentation_file_name
                    # doc_name_human = 'PHANGS_DR%s_CATR%s_hull_ycl_human_ubvi.txt' % (self.data_release, self.catalog_release)
                    doc_name_human = ('hlsp_phangs-cat_hst_multi_hull_multi_v1_%s-%s-%svi.txt' %
                                      (region_type, classify, y_color))

                    if not os.path.isdir(self.catalog_output_path + '/hull'):
                        os.makedirs(self.catalog_output_path + '/hull')
                    doc_file = open(self.catalog_output_path + '/hull/' + doc_name_human, "w")

                    doc_file.writelines('# This file contains the color values for the hulls identified in Maschmann et al.(submitted)) Section 4.4 \n')
                    doc_file.writelines('# The colors are in Vega magnitude. The underlying photometry was MW extintion corrected. \n')
                    doc_file.writelines('# Since these data points are describing a hull, the first and the last datapoints are the same. \n')
                    doc_file.writelines('# Column 1: V-I, colum 2: %s-B  \n' % (y_color[:-1].upper()))
                    for vi, ub in zip(x_hull, y_hull):
                        doc_file.writelines('%.4f %.4f \n' % (vi, ub))

                    doc_file.close()

