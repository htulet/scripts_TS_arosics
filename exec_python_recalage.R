library(reticulate)

path_to_env <- "C:/Users/tulet/miniconda3/envs/test_d"       #Put the name of your python environment here
use_python(path_to_env)

init_path <- "D:/Phenologie/Pipeline/scripts_TS_arosics/__init__.py"

source_python(init_path)



Time_SIFT_process(pathDIR="Y:/RGB/EPOCHS",
                     data_type='RGB',
                     resol_ref=0.05, 
                     crs="EPSG::32622", 
                     site_name = "",
                     out_dir_ortho = "Z:/shared/PhenOBS/Paracou/Metashape/RGB_Broad_Mosaics/Test_sensibilite/matin/ORTHO",
                     out_dir_DEM = "Z:/shared/PhenOBS/Paracou/Metashape/RGB_Broad_Mosaics/Test_sensibilite/matin/DEM",
                     calibrate_col = 'True')
                     #out_dir_project = 'None',
              

complete_arosics_process(path_in = "Z:/shared/PhenOBS/Paracou/Metashape/RGB_Broad_Mosaics/Test_sensibilite/AM/ORTHO",
                         ref_filepath = "Z:/shared/PhenOBS/Paracou/Metashape/RGB_Broad_Mosaics/MNS50cm_Par2019.tif",
                         out_dir_path = "Z:/shared/PhenOBS/Paracou/Metashape/RGB_Broad_Mosaics/Test_sensibilite_rect/Local/AM",
                         corr_type = 'local',
                         dynamic_corr = 'False',
                         apply_matrix = 'False',
                         save_vector_plot='True',
                         save_csv='True')
