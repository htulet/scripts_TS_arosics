import rasterio 
import numpy as np
from pyproj import CRS
import sys
import subprocess
import argparse
import os
from arosics import COREG, COREG_LOCAL, DESHIFTER
import rasterio.features as features
from rasterio.windows import Window
import time
import pandas as pd
from shapely.geometry import Polygon


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
parser.add_argument('--wp', default=(None, None))
parser.add_argument('--grid_res', default=1000)
parser.add_argument('--apply_matrix', default=False)
parser.add_argument('--save_plot', default=False)
parser.add_argument('--save_csv', default=True)
args = parser.parse_args()


def str2bool(v):
    """
    Converts string to bool. Ex : str2bool('True') = True
    """
    if v is None or isinstance(v, bool):
        return v
    if v.lower()=='none':
        return None
    elif v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    elif v in [0, 1]:
        return bool(v)
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def harmonize_crs(input_path, ref_path, check_ref=True):
    """
    Forces two raster files to have the same coordinate system : takes the crs of the reference file and writes it into the input file

    Parameters:
    input file (str): 
        Path to the first image
    ref_path (str): 
        Path to the second image (reference)
    check_ref (bool, optional):  
        If True (default), perform an additional safety measure by rewriting the crs of the reference image aswell. May prevent errors if both files have their crs defined from different libraries (rasterio.CRS and pyproj.CRS)
    """
    with rasterio.open(ref_path) as ds_ref:
        metadata_ref = ds_ref.meta.copy()
        crs_ref = metadata_ref['crs']
        with rasterio.open(input_path) as ds_in:
            img_in = ds_in.read()
            metadata_in = ds_in.meta.copy()
            metadata_in['crs'] = crs_ref
            correction_needed = metadata_in['crs']!=metadata_ref['crs']
            print("harmonization needed : ", correction_needed)
            if correction_needed:
                metadata_in['crs'] = crs_ref
            ds_in.close()  
        if check_ref and not correction_needed:
            img_ref = ds_ref.read()
            ds_ref.close()
            with rasterio.open(ref_path, "w", **metadata_ref) as ds_ref_out:
                ds_ref_out.write(img_ref)
                ds_ref_out.close()
        metadata_ref['crs'] = crs_ref
        ds_ref.close()
    if correction_needed:
        with rasterio.open(input_path, "w", **metadata_in) as ds_out:
            ds_out.write(img_in)
            ds_out.close()


def call_arosics(path_in, path_ref, path_out=None, corr_type = 'global', max_shift=250, max_iter=100, window_size=1500, window_pos = (None, None), mp=None, grid_res=1000, save_csv = True, save_vector_plot = False):
    """
    Calls arosics functions to perform a global or local co-registration between two images. Option to save the coregistrated image, and in the case of a local CoReg, the tie points data and the vector shift map.

    Parameters:
        path_in (str): 
            source path of the target image, i.e. the image to be shifted
        path_ref (str): 
            Path to the refernce image 
        path_out (str, optional): 
            target path of the coregistered image. Defaults to None, in which case nothing will be written to disk
        corr_type (str): 
            Type of co-registration. Either 'global' (default) or 'local'
        max_shift (int): 
            maximum shift distance in reference image pixel units
        max_iter (int): 
            maximum number of iterations for matching (default: 5)
        window_size (int): 
            custom matching window size [pixels] as (X, Y) tuple (default: (256,256))
        window_pos (tuple(int)): 
            custom matching window position as (X, Y) map coordinate in the same projection as the reference image (default: central position of image overlap). Only used when performing global co-registration
        grid_res (int): 
            tie point grid resolution in pixels of the target image (x-direction). Only applies to local co-registration
        mp (int): 
            Number of CPUs to use. If None (default), all available CPUs are used. If mp=1, no multiprocessing is done. 
        save_csv (bool): 
            If True (default), saves the tie points data in a csv file. Has an effect only when performing local co-registration
        save_vector_plot (bool): 
            If True (default), saves the a map of the calculated tie point grid in a JPEG file. Has an effect only when performing local co-registration

    Returns:
        A COREG object containing all info on the calculated shifts
    """
    #CPUs = None if mp else 1
    CPUs = mp if mp is None else int(mp)
    print("CPUs : ", CPUs)
    if corr_type=='global':
        CR = COREG(path_ref, path_in, path_out=path_out, fmt_out="GTIFF", ws=(window_size, window_size), wp=window_pos, max_shift=max_shift, max_iter=max_iter, CPUs=CPUs)
        CR.correct_shifts()
        if save_csv :
            shifts = CR.coreg_info['corrected_shifts_map']
            shift_x, shift_y = shifts['x'], shifts['y']
            df = pd.DataFrame({'Shift_X':[shift_x], 'Shift_Y':[shift_y]})
            df.to_csv(os.path.join(os.path.dirname(path_out), os.path.basename(path_out).split('.')[0] + '_shift.csv'), index=False)

    elif corr_type=='local':
        CR = COREG_LOCAL(path_ref, path_in, path_out=path_out, fmt_out="GTIFF", window_size=(window_size, window_size), max_shift=max_shift, max_iter=max_iter, CPUs=CPUs, grid_res=grid_res)
        CR.correct_shifts()
        if save_csv:
            df = CR.CoRegPoints_table
            df.to_csv(os.path.join(os.path.dirname(path_out), os.path.basename(path_out).split('.')[0] + '_CoRegPoints_table.csv'), index=False)
        if save_vector_plot:
            DPI=300
            vector_scale=15
            CR.view_CoRegPoints(shapes2plot = 'vectors', savefigPath = path_out.split('.')[0] + f"_vector_map_{DPI}DPI.JPEG", savefigDPI=DPI, vector_scale=vector_scale, backgroundIm='tgt')
    return CR

        
def complete_arosics_process(path_in, ref_filepath, out_dir_path, corr_type = 'global', max_shift=250, max_iter=100, grid_res=1000, window_size=None, window_pos = (None, None), mp=None, save_csv = True, save_vector_plot = False, dynamic_corr = False, apply_matrix=False):
    """
    Complete pipeline that uses arosics to perform a global or local co-registration on a file or a group of files located inside a folder. In the case of a local CoReg, option to save the tie points data and the vector shift map.

    :param str path_in: Path to the target image, or to a folder containing multiple target images. Images must be of Geotiff format.
    :param str ref_filepath: Path to the reference image.
    :param str out_dir_path: Directory where the outputs will be saved.
    :param str corr_type: Type of co-registration. Either 'global' (default) or 'local'.
    :param int max_shift: Maximum shift distance in reference image pixel units.
    :param int max_iter: Maximum number of iterations for matching (default: 5).
    :param int grid_res: Tie point grid resolution in pixels of the target image (x-direction). Only applies to local co-registration.
    :param int window_size: Custom matching window size [pixels] as (X, Y) tuple (default: (256,256)).
    :param tuple window_pos: Custom matching window position as (X, Y) map coordinate in the same projection as the reference image (default: central position of image overlap). Only used when performing global co-registration.
    :param int mp: Number of CPUs to use. If None (default), all available CPUs are used. If mp=1, no multiprocessing is done.
    :param bool save_csv: If True (default), saves the tie points data in a CSV file. Has an effect only when performing local co-registration.
    :param bool save_vector_plot: If True (default), saves the a map of the calculated tie point grid in a JPEG file. Has an effect only when performing local co-registration.
    :param bool dynamic_corr: When correcting multiple images, whether or not to use the last corrected image as reference for the next co-registration.
        If False (default), all images are corrected using 'ref_filepath' as the reference image.
        If True, image 1 will use 'ref_filepath' as a reference, then image N (N>=2) will use the corrected version of image N-1 as reference.
    :param bool apply_matrix: When correcting multiple images, whether or not to directly apply the shifts computed for the first image to all the remaining ones, instead of computing the shifts for each one independently. Defaults to False.
        !! Warning !! : Using this option allows faster computing time and better alignment between input images, but will create problems if those images have different bounds.
        In this version, only the intersection of the input images is kept after coregistration.
                        
    :returns: A COREG object (if path_in is a file) or a list of COREG objects (if path_in is a folder) containing all info on the calculated shifts.
    """

    assert corr_type in ['global', 'local']
    rm_temp_files = False
    dynamic_corr = str2bool(dynamic_corr)
    apply_matrix = str2bool(apply_matrix)
    save_vector_plot = str2bool(save_vector_plot)
    save_csv = str2bool(save_csv)
    mp = mp if mp is None else int(mp)
    #Set default values for window_size
    if corr_type == 'global':
        if window_size is None :
            window_size = 1500
        grid_res = ""
    elif corr_type == 'local':
        if window_size is None :
            window_size = 4000 

    if not os.path.exists(out_dir_path):
        os.mkdir(out_dir_path)

    extensions = ('.tif', '.tiff', '.TIF', '.TIFF')

    if os.path.isfile(path_in):
        if not path_in.endswith(extensions):
            raise ValueError(f"The specified file '{path_in}' must be of GeoTiff format")
        else:
            harmonize_crs(path_in, ref_filepath)
            path_out = os.path.join(out_dir_path, path_in.split('/')[-1].split('\\')[-1].split('.')[0] + f'_aligned_{corr_type}.tif')
            CR = call_arosics(path_in, ref_filepath, path_out=path_out, corr_type=corr_type, mp=mp, window_size=window_size, window_pos=window_pos, max_shift=max_shift, max_iter=max_iter, grid_res=grid_res, save_vector_plot=save_vector_plot, save_csv=save_csv)
            return CR
            
    elif os.path.isdir(path_in):
        list_CR = []
        files = [file for file in sorted(os.listdir(path_in)) if file.endswith(extensions)]
        print("files : ", files)

        if len(files)==0:
            raise ValueError(f"The specified directory '{path_in}' does not contain any GeoTiff files.")

        elif dynamic_corr or not apply_matrix :
            for i in range(len(files)):
                file = files[i]
                current_file_path = os.path.join(path_in, file)
                harmonize_crs(current_file_path, ref_filepath, check_ref = True if i==0 else False)
                path_out = os.path.join(out_dir_path, file.split('.')[0].replace("_temp", "") + f'_aligned_{corr_type}.tif')
                CR = call_arosics(current_file_path, ref_filepath, path_out=path_out, corr_type=corr_type, mp=mp, window_size=window_size, window_pos=window_pos, max_shift=max_shift, max_iter=max_iter, grid_res=grid_res, save_vector_plot=save_vector_plot, save_csv=save_csv)
                list_CR.append(CR)
                if dynamic_corr:
                    ref_filepath = path_out
                
        else:
            hlist, blist, glist, dlist, band_count = [], [], [], [], []
            for file in files:
                meta = rasterio.open(os.path.join(path_in, file)).meta.copy()
                tf = list(meta['transform'])
                hlist.append(tf[5])
                blist.append(tf[5] + tf[4]*meta['height'])
                glist.append(tf[2])
                dlist.append(tf[2] + tf[0]*meta['width'])
                band_count.append(meta['count'])
            
            
            if not (all(x == hlist[0] for x in hlist) and all(x == glist[0] for x in glist)):     #and (np.max(dlist)-np.max(glist)) < 0.9*(np.max(dlist)-np.min(glist)) and (np.min(hlist)-np.max(blist)) < 0.9*(np.max(hlist)-np.max(blist))  #if apply_mask
                rm_temp_files=True
                mask_coords = [np.min(glist), np.min(blist), np.max(dlist), np.max(hlist)]
                geom = Polygon([(mask_coords[0], mask_coords[1]), (mask_coords[0], mask_coords[3]), (mask_coords[2], mask_coords[3]), (mask_coords[2], mask_coords[1])])
                num_cols = int(np.ceil((mask_coords[2]-mask_coords[0]) / tf[0]))
                num_rows = int(np.ceil((mask_coords[1]-mask_coords[3]) / tf[4]))
                print("num cols : ", num_cols)
                print("num rows : ", num_rows)
                print("Mask : ", geom.bounds)
                #All images are padded to fit the geometry of the mask
                for i in range(len(files)) :
                    file = files[i]
                    rast = rasterio.open(os.path.join(path_in, file))
                    img = rast.read()
                    print("image shape : ", img.shape)

                    padded_data = np.ones((band_count[i], num_rows, num_cols), dtype=img.dtype)*255
                    
                    # Calculate the offset to pad the smaller raster
                    row_offset = int((geom.bounds[3] - rast.bounds.top) / abs(rast.res[0]))
                    col_offset = int((rast.bounds.left - geom.bounds[0]) / abs(rast.res[1]))
                    
                    # Define the window to copy the smaller raster into the padded data array
                    window = Window(col_offset, row_offset, img.shape[2], img.shape[1])

                    # Copy the smaller raster into the padded data array
                    padded_data[:, window.row_off : window.row_off+window.height, window.col_off : window.col_off+window.width] = img
                    # Update the metadata for the padded raster
                    target_tf = rasterio.Affine(list(rast.meta.copy()['transform'])[0], 0.0, mask_coords[0], 0.0, list(rast.meta.copy()['transform'])[4], mask_coords[3])
                    profile = rast.profile
                    profile.update(width=padded_data.shape[2], height=padded_data.shape[1], transform=target_tf,
                                    bounds=geom.bounds, count = band_count[i])
                    print("padded img shape : ", padded_data.shape)
                    out_path = os.path.join(path_in, f"{file.split('.')[0]}_temp.tif")
                    files[i] = f"{file.split('.')[0]}_temp.tif"
                    rast.close()
                    # Write the padded data to a new raster file
                    with rasterio.open(out_path, "w", **profile) as dst:
                        dst.write(padded_data)
                        dst.close()

            first_file = files[0]
            harmonize_crs(os.path.join(path_in, first_file), ref_filepath)
            path_out = os.path.join(out_dir_path, first_file.split('.')[0].replace("_temp", "") + f'_aligned_{corr_type}.tif')
            CR = call_arosics(os.path.join(path_in, first_file), ref_filepath, path_out=path_out, corr_type=corr_type, mp=mp, window_size=window_size, window_pos=window_pos, max_shift=max_shift, max_iter=max_iter, grid_res=grid_res, save_vector_plot=save_vector_plot, save_csv=save_csv)             
            list_CR.append(CR)
            for file in files[1:]:
                current_file_path = os.path.join(path_in, file)
                harmonize_crs(current_file_path, ref_filepath, check_ref=False)
                path_out = os.path.join(out_dir_path, file.split('.')[0].replace("_temp", "") + f'_aligned_{corr_type}.tif')
                CR = DESHIFTER(current_file_path, CR.coreg_info, path_out=path_out, fmt_out="GTIFF")
                CR.correct_shifts() 
                list_CR.append(CR)
        if rm_temp_files:
            for file in files:
                os.remove(os.path.join(path_in, file))
        
        return list_CR
    
    else:
        raise ValueError(f"The specified path '{path_in}' is not a file nor a directory.")


if __name__ == '__main__':
    
    complete_arosics_process(path_in = args.path_in,
                             ref_filepath = args.ref_filepath, 
                             out_dir_path = args.out_dir_path, 
                             corr_type = args.corr_type, 
                             mp = args.mp,
                             max_shift = args.max_shift,
                             max_iter = args.max_iter,
                             grid_res = args.grid_res,
                             window_pos = args.wp,
                             window_size = args.ws,
                             dynamic_corr = args.dynamic_corr,
                             apply_matrix = args.apply_matrix,
                             save_csv = args.save_csv,
                             save_vector_plot = args.save_plot
                             )
    
