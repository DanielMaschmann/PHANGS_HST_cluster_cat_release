"""
Script to create full data release for PHANGS-HST cluster catalogs
"""

import os
from data_release_routines import DataReleaseRoutines

try:
    import data_access_config
except ImportError:
    raise ImportError('No data_access_config.py file found. This file needs to be created to specify data paths to '
                      'internal release etc. An example is shown in the file data_access_config_example.py ')

_author = 'Daniel Maschmann'
_contact = 'danielmaschmann(at)arizona(dot)edu'
_data_release = 4
_catalog_release = 3

# get local data path structure
config_dict = data_access_config.config_dict
# check if data_path_exists
if not os.path.isdir(config_dict['path2ir']):
    raise FileNotFoundError('The folder ', config_dict['path2ir'], ' does not exists on this machine.')
if config_dict['artifact_removal_flag'] and not os.path.isdir(config_dict['path2artifact']):
    raise FileNotFoundError('The folder ', config_dict['path2artifact'], ' does not exists on this machine ',
                            'You need either specify path2artifact or switch artifact_removal_flag to False.')

# Create object of DataReleaseRoutine class and give it all variables from the config_dict
data_release_ojt = DataReleaseRoutines(**config_dict)

# create the final data release (all the catalogues)
data_release_ojt.create_final_data_release()

# create documentaion
data_release_ojt.create_documentation()

data_release_ojt.create_hull_files()

# data_release_ojt.print_content_table(table_name='tab2', max_length_description=70)
data_release_ojt.print_content_table(table_name='tab1', max_length_description=70)
