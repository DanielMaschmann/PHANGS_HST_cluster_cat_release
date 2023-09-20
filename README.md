# PHANGS_HST_cluster_cat_release
The main purpose of this repository is to create the PHANGS-HST cluster catalogs for the public release

In general the data release should be produced by running the pythonscript cluster_cat_dr/create_data_release.py

Before doing this all data paths need to be configured in the file cluster_cat_dr/data_access_config.py. 
Since this file is depending on your local data management, you have to produce this on your own and this file is not 
tracked by the version control. You can orientate yourself with the file cluster_cat_dr/data_access_config_example.py 

Once you run the creation of the data release the catalog files will be dropped in the specified folder. 
Furthermore, metadata such as documentaries and tables for the articles will be written to metadata_output folder.
For the artifact removal step, you can specify the artifact_plotting flag, which will produce cutout images to visually 
inspect the objects excluded from the final data release.  
