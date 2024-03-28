library(reticulate)

path_to_env <- "C:/Users/U051-U219/miniconda3/envs/env_R"       #Put the name of your python environment here
use_python(path_to_env)

arosics_script_path <- "D:/pipeline_TS/scripts_TS_arosics/arosics_chain.py"
TS_script_path <- '//amap-data.cirad.fr/work/users/HadrienTulet/scripts_TS_arosics/2023-07-13_TimeSIFT_scripts_auto.py'

source_python(arosics_script_path)
source_python(TS_script_path)



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