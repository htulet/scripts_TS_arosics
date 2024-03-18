library(reticulate)
path_to_env <- "C:/Users/U051-U219/miniconda3/envs/env_R"
use_python(path_to_env)

arosics_script_path <- "arosics_chain.py"
TS_script_path <- '2023-07-13_TimeSIFT_scripts_auto.py'

source_python(arosics_script_path)
source_python(TS_script_path)



Time_SIFT_process(pathDIR="Y:/RGB/EPOCHS",
                     data_type='RGB',
                     resol_ref=0.05, 
                     crs="EPSG::32622", 
                     site_name = "",
                     out_dir_ortho = "Z:/users/HadrienTulet/keypoints_test/ORTHO",
                     out_dir_DEM = "Z:/users/HadrienTulet/keypoints_test/DEM",
                     calibrate_col = 'True')
                     #out_dir_project = 'None',
              

complete_arosics_process(path_in = "Z:/users/HadrienTulet/tests_arosics/data",
                         ref_filepath = "Z:/shared/PhenOBS/Paracou/Metashape/RGB_Broad_Mosaics/PARACOU2023_DSM_50cm_filt.tif",
                         out_dir_path = "Z:/users/HadrienTulet/tests_arosics/apply_matrix_no_mask_v2",
                         corr_type = 'global',
                         dynamic_corr = 'False',
                         apply_matrix = 'True')

