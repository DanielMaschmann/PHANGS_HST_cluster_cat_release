"""
Script to create full data release for PHANGS-HST cluster catalogs

"""
import numpy as np

_author = 'Daniel Maschmann'
_contact = 'danielmaschmann(at)arizona(dot)edu'
_script_version = 1.0
_data_release = 2


import os
from pathlib import Path
from sys import argv
import catalog_content
from astropy.io import fits

# Usage:
# 1 file path to internal release
# 2


# specify the used Version
hst_cc_ver = 'SEDfix_final_test_catalogs'

# set data path to internal release data
if len(argv) == 1:
    path2ir = '/home/benutzer/data/PHANGS_products/HST_catalogs'
    path2artifact = '/home/benutzer/data/PHANGS_products/HST_catalogs/Artifact_Removal/AR0'
else:
    path2ir = argv[1]
    path2artifact = argv[2]


# loop over all galaxies
# target = 'ic1954'
# target = 'ngc1559'
target = 'ngc1566'
# loop over human and ML
classify = 'ml'
# loop over class12 and class 3
cl_class = 'class12'


# get the catalog
cat_file_name_ir = Path('SEDfix_PHANGS_IR4_%s_Ha1_inclusiveGCcc_inclusiveGCclass_phangs_hst_v1p2_%s_%s.fits'
                       % (target, classify, cl_class))
cluster_dict_path = Path(path2ir) / Path(hst_cc_ver)
file_path_ir = cluster_dict_path / cat_file_name_ir
# check if file exists
if not os.path.isfile(file_path_ir):
    print(file_path_ir, ' not found ! ')
    raise FileNotFoundError('there is no HST cluster catalog for the target ', target,
                            ' make sure that the file ', file_path_ir, ' exists')
# load table
table_ir = fits.open(file_path_ir)[1].data
col_names_ir = table_ir.names

# get artifacts
file_name_artifact = Path('%s_artifacts.fits' % target)
file_path_artifact = path2artifact / file_name_artifact
# check if file exists
if not os.path.isfile(file_path_artifact):
    print(file_path_artifact, ' not found ! ')
    raise FileNotFoundError('there is no HST cluster catalog for the target ', target,
                            ' make sure that the file ', file_path_artifact, ' exists')
# get artifact tables
table_artifact = fits.open(file_path_artifact)[1].data
col_names_artifact = table_artifact.names

print(col_names_ir)
print(col_names_artifact)

print(table_artifact['PHANGS_CLUSTER_CLASS_ML_VGG'])
print(table_artifact['PHANGS_CLUSTER_CLASS_HUMAN'])

# create artifact mask
artifact_mask_table_ir = np.zeros(len(table_ir), dtype=bool)
for artifact_index in range(len(table_artifact)):

    # artifact_in_table_ir = (table_ir['ID_PHANGS_CLUSTERS_v1p2'] ==
    #                         table_artifact['ID_PHANGS_CLUSTERS_v1p2'][artifact_index])
    artifact_in_table_ir = ((table_ir['ID_PHANGS_CLUSTERS_v1p2'] ==
                             table_artifact['ID_PHANGS_CLUSTERS_v1p2'][artifact_index]) &
                            (table_ir['PHANGS_X'] == table_artifact['PHANGS_X'][artifact_index]) &
                            (table_ir['PHANGS_Y'] == table_artifact['PHANGS_Y'][artifact_index]) &
                            (table_ir['PHANGS_RA'] == table_artifact['PHANGS_RA'][artifact_index]) &
                            (table_ir['PHANGS_DEC'] == table_artifact['PHANGS_DEC'][artifact_index]))
    if sum(artifact_in_table_ir) == 0:
        print('No cross match for Object ', table_artifact['ID_PHANGS_CLUSTERS_v1p2'][artifact_index])
        print('PHANGS_CLUSTER_CLASS_ML_VGG ', table_artifact['PHANGS_CLUSTER_CLASS_ML_VGG'][artifact_index],
              'PHANGS_CLUSTER_CLASS_HUMAN ', table_artifact['PHANGS_CLUSTER_CLASS_HUMAN'][artifact_index],
              'NEW_CLASS ', table_artifact['NEW_CLASS'][artifact_index])
    artifact_mask_table_ir += artifact_in_table_ir

print('size of artifact table ', len(table_artifact))
print('Number cross matched artifacts in the IR table ', sum(artifact_mask_table_ir))

