"""
specification configurations and data structures for the local machine
"""

""" all information and configurations which is specific for the local machine must be gathered in this dictionary """
config_dict = {
    # Internal release
    # specify the used Version. This is important since the naming of catalog files can change between IR versions.
    'hst_cc_ver': 'SEDfix_final_test_catalogs',
    # path to internal release
    # 'path2ir': '/home/benutzer/data/PHANGS_products/HST_catalogs/SEDfix_final_test_catalogs',
    'path2ir': '/home/benutzer/data/PHANGS_products/HST_catalogs/SEDfix_Ha1_inclusiveGCcc_inclusiveGCclass',

    # Data output
    # data release
    'data_release': '4',
    # catalog release
    'catalog_release': '2',
    # path to Data Release 4, Catalog Release 2 at MAST
    'catalog_output_path': '/home/benutzer/data/PHANGS_products/HST_catalogs/phangs_hst_cc_dr4_cr2',

    # Identified artifact in the internal release
    # path to tables/catalogs of artifacts
    'path2artifact': '/home/benutzer/data/PHANGS_products/HST_catalogs/Artifact_Removal/AR1',
    # a re-evaluation of questionable artefacts.
    'path2questionable_artifacts': '/home/benutzer/data/PHANGS_products/HST_catalogs/Artifact_Removal/'
                                   'Questionable Artifacts - questionable_artifacts.csv',
    # masks for diffraction spikes
    'path2diffraction_spike_masks': '/home/benutzer/data/PHANGS_products/HST_catalogs/Artifact_Removal/'
                                    'diffraction_spikes/fits',
    # limit for V-I color and concentration index
    'v_i_color_lim': 2.0,
    'ci_lim': 1.45,

    # flag if artifact removal step should be done
    'artifact_removal_flag': True,
    # flag to raise error if a file is not found. This can be handy for early stages in the development as not every
    # target has yet an artifact catalog. For the final version, however, this should be set to True
    'artifact_rais_file_not_found_flag': True,
    # flag to plot the removed artefacts
    'plot_removed_artifacts_flag': False,

    # additional
    'existing_artifact_removal_flag': True,

    # for artifact visualization
    'hst_data_path': '/media/benutzer/Sicherung/data/phangs_hst',
    'nircam_data_path': '/media/benutzer/Sicherung/data/phangs_jwst',
    'miri_data_path': '/media/benutzer/Sicherung/data/phangs_jwst',
    'hst_data_ver': 'v1.0',
    'nircam_data_ver': 'v0p9',
    'miri_data_ver': 'v0p9',

}
