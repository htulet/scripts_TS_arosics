library(reticulate)

path_to_env <- "D:/tulet/envs/pipeline"       #Put the name of your python environment here
use_python(path_to_env)

init_path <- "D:/Phenologie/Pipeline/scripts_TS_arosics/__init__.py"

source_python(init_path)



Time_SIFT_process(pathDIR='Z:/users/HadrienTulet/data_pipe',
                     data_type='RGB',
                     out_dir_ortho = 'Z:/users/HadrienTulet/data_pipeline/ORTHO',
                     calibrate_col = True)
                     #out_dir_project = 'None',
              

complete_arosics_process(path_in = "Z:/shared/PhenOBS/Paracou/Metashape/RGB_Broad_Mosaics/Test_sensibilite/AM/ORTHO",
                         ref_filepath = "Z:/shared/PhenOBS/Paracou/Metashape/RGB_Broad_Mosaics/MNS50cm_Par2019.tif",
                         out_dir_path = "Z:/shared/PhenOBS/Paracou/Metashape/RGB_Broad_Mosaics/Test_sensibilite_rect/Local/AM",
                         corr_type = 'local',
                         dynamic_corr = 'False',
                         apply_matrix = 'False',
                         save_vector_plot='True',
                         save_csv='True')



script_path <- "D:/Phenologie/Pipeline/scripts_TS_arosics/subprocess_test.py"
source_python(script_path)
