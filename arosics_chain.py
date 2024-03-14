import rasterio 
import numpy as np
from pyproj import CRS
import sys
import subprocess
import argparse
import os
from arosics import COREG, COREG_LOCAL, DESHIFTER
import rasterio.features as features
import time
from shapely import Polygon

parser = argparse.ArgumentParser()
parser.add_argument('--path_in')
parser.add_argument('--ref_filepath')
parser.add_argument('--out_dir_path')
parser.add_argument('--corr_type', default='global')
parser.add_argument('--dynamic_corr', default=False)
parser.add_argument('--mp', default=1)
parser.add_argument('--max_shift', default=250)
parser.add_argument('--max_iter', default=100)
parser.add_argument('--ws', default=None)
parser.add_argument('--grid_res', default=1000)
parser.add_argument('--apply_matrix', default=True)
args = parser.parse_args()

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def harmonize_crs(input_path, ref_path, check_ref=True):
    with rasterio.open(ref_path) as ds_ref:
        metadata_ref = ds_ref.meta.copy()
        crs_ref = metadata_ref['crs']
        with rasterio.open(input_path) as ds_in:
            img_in = ds_in.read()
            metadata_in = ds_in.meta.copy()
            metadata_in['crs'] = crs_ref
            ds_in.close()  
        if check_ref:
            img_ref = ds_ref.read()
            ds_ref.close()
            with rasterio.open(ref_path, "w", **metadata_ref) as ds_ref_out:
                ds_ref_out.write(img_ref)
                ds_ref_out.close()
        metadata_ref['crs'] = crs_ref
        ds_ref.close()
    with rasterio.open(input_path, "w", **metadata_in) as ds_out:
        ds_out.write(img_in)
        ds_out.close()

def call_arosics(path_in, path_ref, path_out, corr_type = 'global', max_shift=250, max_iter=100, window_size=None, mp=0, grid_res=1000):
    assert corr_type in ['global', 'local']
    if corr_type == 'global':
        if window_size is None :
            window_size = 1500
        grid_res = ""
    elif corr_type == 'local':
        if window_size is None :
            window_size = 4000 
    command = f"arosics {corr_type} -max_shift {max_shift} -max_iter {max_iter} -ws {window_size} {window_size} -mp {mp} -fmt_out GTIFF -o {path_out} {path_ref} {path_in} {grid_res}"
    print(command)
    subprocess.run(command, shell=True)

def complete_arosics_process_v1(path_in, ref_filepath, out_dir_path, corr_type = 'global', dynamic_corr = False):
    dynamic_corr = str2bool(dynamic_corr)
    extensions = ('.tif', '.tiff', '.TIF', '.TIFF')
    if os.path.isfile(path_in):
        harmonize_crs(path_in, ref_filepath)
        path_out = os.path.join(out_dir_path, path_in.split('/')[-1].split('\\')[-1].split('.')[0] + '_aligned.tif')
        call_arosics(path_in, ref_filepath, path_out, corr_type)
    elif os.path.isdir(path_in):
        for file in sorted(os.listdir(path_in)):
            if file.endswith(extensions):
                current_file_path = os.path.join(path_in, file)
                harmonize_crs(current_file_path, ref_filepath)
                path_out = os.path.join(out_dir_path, file.split('.')[0] + '_aligned.tif')
                call_arosics(current_file_path, ref_filepath, path_out, corr_type)
                if dynamic_corr :
                    ref_filepath = path_out
                
    else:
        raise ValueError(f"The specified path '{path_in}' is not a file nor a directory.")
        
def complete_arosics_process(path_in, ref_filepath, out_dir_path, corr_type = 'global', dynamic_corr = False, apply_matrix=True, max_shift=250, max_iter=100, window_size=None, mp=1, grid_res=1000):
    assert corr_type in ['global', 'local']
    rm_temp_files = False
    dynamic_corr = str2bool(dynamic_corr)
    apply_matrix = str2bool(apply_matrix)
    print(bool(apply_matrix))
    #Set default values for window_size
    if corr_type == 'global':
        if window_size is None :
            window_size = 1500
        grid_res = ""
    elif corr_type == 'local':
        if window_size is None :
            window_size = 4000 

    dynamic_corr = str2bool(dynamic_corr)

    extensions = ('.tif', '.tiff', '.TIF', '.TIFF')

    if os.path.isfile(path_in):
        if not file.endswith(extensions):
            raise ValueError(f"The specified file '{path_in}' must be of GeoTiff format")
        else:
            harmonize_crs(path_in, ref_filepath)
            path_out = os.path.join(out_dir_path, path_in.split('/')[-1].split('\\')[-1].split('.')[0] + f'_aligned_{corr_type}.tif')
            call_arosics(path_in, ref_filepath, path_out, corr_type, mp=mp, window_size=window_size, max_shift=max_shift, max_iter=max_iter, grid_res=grid_res)

    elif os.path.isdir(path_in):
        files = [file for file in sorted(os.listdir(path_in)) if file.endswith(extensions)]
        print("files : ", files)

        if len(files)==0:
            raise ValueError(f"The specified directory '{path_in}' does not contain any GeoTiff files.")
        
        elif not apply_matrix:
            for file in files:
                current_file_path = os.path.join(path_in, file)
                harmonize_crs(current_file_path, ref_filepath)
                path_out = os.path.join(out_dir_path, file.split('.')[0].replace("_temp", "") + f'_aligned_{corr_type}.tif')
                call_arosics(current_file_path, ref_filepath, path_out, corr_type, mp=mp, window_size=window_size, max_shift=max_shift, max_iter=max_iter, grid_res=grid_res)
                
        else:
            hlist, blist, glist, dlist = [], [], [], []
            for file in files:
                meta = rasterio.open(os.path.join(path_in, file)).meta.copy()
                tf = list(meta['transform'])
                hlist.append(tf[5])
                blist.append(tf[5] + tf[4]*meta['height'])
                glist.append(tf[2])
                dlist.append(tf[2] + tf[0]*meta['width'])
            """
            if not (all(x == hlist[0] for x in hlist) and all(x == glist[0] for x in glist)):     #and (np.max(dlist)-np.max(glist)) < 0.9*(np.max(dlist)-np.min(glist)) and (np.min(hlist)-np.max(blist)) < 0.9*(np.max(hlist)-np.max(blist))
                rm_temp_files=True
                mask_coords = [np.max(glist), np.min(blist), np.max(dlist), np.min(hlist)]
                geom = Polygon([(mask_coords[0], mask_coords[1]), (mask_coords[0], mask_coords[3]), (mask_coords[2], mask_coords[3]), (mask_coords[2], mask_coords[1])])
                print("Mask : ", geom.bounds)
                #All images are cropped to fit the geometry of the mask
                for i in range(len(files)) :
                    file = files[i]
                    rast = rasterio.open(os.path.join(path_in, file))
                    img = rast.read()
                    mask = features.rasterize([geom],
                                                out_shape = rast.shape,
                                                transform = rast.transform,
                                                default_value=1)
                    new_img = np.tile(mask, (3,1,1)) * img
                    nonzero_rows, nonzero_cols = np.any(new_img != 0, axis=(0,2)), np.any(new_img != 0, axis=(0,1))
                    red_new_img = new_img[:, nonzero_rows][:,:, nonzero_cols]
                    print(red_new_img.shape)
                    file_meta = rast.meta.copy()
                    target_tf = rasterio.Affine(list(rast.meta.copy()['transform'])[0], 0.0, mask_coords[0], 0.0, list(rast.meta.copy()['transform'])[4], mask_coords[3])
                    file_meta['transform'] = target_tf       
                    #file_meta['crs'] = rasterio.open(ref_filepath).meta['crs']
                    file_meta['height'], file_meta['width'] = red_new_img.shape[1:]
                    out_path = os.path.join(path_in, f"{file.split('.')[0]}_temp.tif")
                    files[i] = f"{file.split('.')[0]}_temp.tif"
                    rast.close()
                    with rasterio.open(out_path, "w", **file_meta) as ds_ref:
                            ds_ref.write(red_new_img)
                            ds_ref.close()
            """
            if dynamic_corr :
                for file in files:
                    current_file_path = os.path.join(path_in, file)
                    harmonize_crs(current_file_path, ref_filepath)
                    path_out = os.path.join(out_dir_path, file.split('.')[0].replace("_temp", "") + f'_aligned_{corr_type}.tif')
                    call_arosics(current_file_path, ref_filepath, path_out, corr_type, mp=mp, window_size=window_size, max_shift=max_shift, max_iter=max_iter, grid_res=grid_res)
                    ref_filepath = path_out

            else:
                first_file = files[0]
                harmonize_crs(os.path.join(path_in, first_file), ref_filepath)
                path_out = os.path.join(out_dir_path, first_file.split('.')[0].replace("_temp", "") + f'_aligned_{corr_type}.tif')              
                CPUs = None if mp else 1
                if corr_type=='global':
                    CR = COREG(ref_filepath, os.path.join(path_in, first_file), path_out=path_out, fmt_out="GTIFF", ws=(window_size, window_size), max_shift=max_shift, max_iter=max_iter, CPUs=CPUs)
                elif corr_type=='local':
                    CR = COREG_LOCAL(ref_filepath, os.path.join(path_in, first_file), path_out=path_out, fmt_out="GTIFF", window_size=(window_size, window_size), max_shift=max_shift, max_iter=max_iter, CPUs=CPUs, grid_res=grid_res)
                CR.correct_shifts()
                for file in files[1:]:
                    current_file_path = os.path.join(path_in, file)
                    harmonize_crs(current_file_path, ref_filepath, check_ref=False)
                    path_out = os.path.join(out_dir_path, file.split('.')[0].replace("_temp", "") + f'_aligned_{corr_type}.tif')
                    DESHIFTER(current_file_path, CR.coreg_info, path_out=path_out, fmt_out="GTIFF").correct_shifts() 
        
        if rm_temp_files:
            for file in files:
                os.remove(os.path.join(path_in, file))
        
    else:
        raise ValueError(f"The specified path '{path_in}' is not a file nor a directory.")
"""   
if __name__ == '__main__':
    complete_arosics_process(path_in = args.path_in, ref_filepath = args.ref_filepath, out_dir_path = args.out_dir_path, corr_type = args.corr_type, dynamic_corr=str2bool(args.dynamic_corr), apply_matrix=str2bool(args.apply_matrix))
"""    
 

