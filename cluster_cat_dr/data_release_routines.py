"""
All production routines for the final PHANGS-HST cluster catalog data release will be gathered here
"""

import os
from pathlib import Path
import numpy as np
from astropy.io import fits
from cat_info import CatalogInfo
from astropy.table import Column, hstack


class DataReleaseRoutines(CatalogInfo):
    def __init__(self, hst_cc_ver, path2ir, catalog_output_path,
                 path2artifact=None, artifact_removal_flag=False, artifact_rais_file_not_found_flag=True,
                 plot_removed_artifacts_flag=False, data_release=4, catalog_release=2):
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

        # load constructor of parent class
        super().__init__()

    def create_final_data_release(self):
        """
        Function to run all steps to create the PHANGs-HST cluster catalog data release
        """
        # create catalogs
        for target in self.phangs_hst_target_list:
            # get artifact removal table
            table_artifact = self.get_artifact_cat(target=target)

            for classify in ['human', 'ml']:
                for cl_class in ['class12', 'class3']:
                    print('creating cluster catalog for ', target, classify, cl_class)
                    table_ir = self.get_ir_cat(target, classify, cl_class)
                    # create artifact mask
                    artifact_mask_table_ir = np.zeros(len(table_ir), dtype=bool)
                    if table_artifact is not None:
                        for artifact_index in range(len(table_artifact)):
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
                            artifact_mask_table_ir += artifact_in_table_ir
                    print('number cross matched artefacts ', sum(artifact_mask_table_ir))

                    # plot artifacts !!!
                    # TBD !!!!!

                    # apply artifact mask
                    table_ir = table_ir[np.invert(artifact_mask_table_ir)]

                    # create observation table
                    obs_table = self.create_obs_table(target=target, table=table_ir)

                    # save table
                    # check if table already exists
                    if not os.path.isdir(self.catalog_output_path):
                        os.mkdir(self.catalog_output_path)
                    # get table name
                    table_name = self.get_data_release_table_name(target=target, classify=classify, cl_class=cl_class)
                    # save table
                    obs_table.write(Path(self.catalog_output_path) / Path(table_name), overwrite=True)

    def get_data_release_table_name(self, target, classify, cl_class):
        """
        Function to create final catalog name
        ----------
        target : str
            name of PHANGS-HST target. Must be in self.phangs_hst_target_list
        classify : str
            classification `human` or `ml`
        cl_class : str
            class group specification either `class12` or `class3`
        Returns
        -------
        table_name : str
            table name
        """
        return ('PHANGS_HST_cluster_catalog_dr_%s_cat_release_%s_%s_%s_%s.fits' %
                (self.data_release, self.catalog_release, target, classify, cl_class))

    def create_obs_table(self, target, table):
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
        column_name_list = self.get_obs_table_column_list(target)
        # create table
        obs_table = None
        for col_name in column_name_list:
            column_content = table[col_name]
            if self.cat_info[col_name]['unit'] is not None:
                column_content *= self.cat_info[col_name]['unit']
            column = Column(data=column_content,
                            name=self.cat_info[col_name]['col_name'],
                            dtype=column_content.dtype,
                            description=self.cat_info[col_name]['doc_comment']
                            )

            if obs_table is None:
                obs_table = column
            else:
                obs_table = hstack([obs_table, column])

        return obs_table


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
            cat_file_name_ir = Path('SEDfix_PHANGS_IR4_%s_Ha1_inclusiveGCcc_inclusiveGCclass_phangs_hst_v1p2_%s_%s.fits'
                                    % (target, classify, cl_class))
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
        table_ir = fits.open(file_path_ir)[1].data

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
        # get artifact table file path
        file_path_artifact = self.path2artifact / Path('%s_artifacts.fits' % target)
        # check if file exists
        if (not os.path.isfile(file_path_artifact)) & self.artifact_rais_file_not_found_flag:
            print(file_path_artifact, ' not found ! ')
            raise FileNotFoundError('there is no artigfact table for the target ', target,
                                    ' make sure that the file ', file_path_artifact, ' exists.')
        elif not os.path.isfile(file_path_artifact):
            print('No artifact table found for ', target)
            return None
        else:
            return fits.open(file_path_artifact)[1].data
