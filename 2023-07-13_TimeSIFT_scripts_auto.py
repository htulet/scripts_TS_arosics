# Script for Time-SIFT multi-temporal images alignment
#
# This is python script for PhotoScan Pro, based on the following publication :
#
# D. Feurer, F. Vinatier, Joining multi-epoch archival aerial images in a single SfM block allows 3-D change detection with almost exclusively image information, ISPRS Journal of Photogrammetry and Remote Sensing, Volume 146, 2018, Pages 495-506, ISSN 0924-2716, https://doi.org/10.1016/j.isprsjprs.2018.10.016. (http://www.sciencedirect.com/science/article/pii/S0924271618302946)


"""
script manuel : scan.app.document
script auto : scan.Document(), keep the same thorough whole process
"""

import os
from os import path
import re
import numpy as np
import Metashape as scan
import time
import argparse
import shutil

#scan.License().activate('your_license_key')

parser = argparse.ArgumentParser()
parser.add_argument('--crs')
parser.add_argument('--pathDIR')
parser.add_argument('--out_dir_ortho')
parser.add_argument('--out_dir_dem')
parser.add_argument('--resol_ref')
parser.add_argument('--data_type')
args = parser.parse_args()


def str2bool(v):
    if v is None or isinstance(v, bool):
        return v
    if v.lower()=='none':
        return None
    elif v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def add_all_chunks(doc = scan.Document(), pathDIR=None):
    print(pathDIR)
    os.chdir(pathDIR)
    epochs = os.listdir(pathDIR)
    # we select only the non-empty subfolders 
    epochs = [ep for ep in epochs if not os.path.isfile(ep) and os.listdir(ep)]
    # We remove all existing chunks and add them one by one
    for chk in doc.chunks:
        doc.remove(chk)
    for ep in epochs: 
        add_TimeSIFT_chunk(doc, epoch_name=ep)


def add_TimeSIFT_chunk(doc = scan.Document(), epoch_name=""):
    if len([chk for chk in doc.chunks if re.search(epoch_name,chk.label) is not None])==0:
        doc.addChunk()
        chunk=doc.chunks[len(doc.chunks)-1]
        chunk.label=epoch_name
        [f for f in os.listdir(epoch_name) if os.path.isfile(os.path.join(epoch_name, f))]
        dirName=epoch_name
        listOfFiles = list()
        for (dirpath, dirnames, filenames) in os.walk(dirName):
            listOfFiles += [os.path.join(dirpath, file) for file in filenames]

        chunk.addPhotos(listOfFiles)
        for cam in chunk.cameras:
            cam.label = (str(chunk.label) + "_EPOCH_" + cam.label)



def add_all_MS_photos(doc = scan.Document(), pathDIR=None):
    print(pathDIR)
    os.chdir(pathDIR)
    for chk in doc.chunks:
        doc.remove(chk)
    doc.addChunk()
    chunk=doc.chunks[len(doc.chunks)-1]
    chunk.label="TimeSIFT"
            
    epochs=os.listdir(pathDIR)
    epochs = [ep for ep in epochs if not os.path.isfile(ep) and os.listdir(ep)]
    print(epochs)
    for i in range(len(epochs)):
        L = len(chunk.cameras)
        epoch_name = epochs[i]
        # We create one camera group per epoch
        chunk.addCameraGroup()
        cam_group = chunk.camera_groups[-1]
        cam_group.label = epoch_name
        
        dirName= os.path.join(pathDIR, epoch_name)
        listOfFiles = list()
        for (dirpath, dirnames, filenames) in os.walk(dirName):
            listOfFiles += [os.path.join(dirpath, file) for file in filenames]
        
        # We add all photos from the same epoch into the corresponding camera group 
        chunk.addPhotos(listOfFiles, layout = scan.MultiplaneLayout, group = i)  
        for cam in chunk.cameras[L:]:
            cam.label = (str(epoch_name) + "_EPOCH_" + cam.label)  



def merge_chunk_TimeSIFT(doc = scan.Document()):
    start_time = time.time()
    for chk_sel in doc.chunks:
        chunk_non_aligned = [chk.key for chk in doc.chunks if re.search("TimeSIFT", chk.label) is not None]
        chunk_non_aligned.append(chk_sel.key)
        doc.mergeChunks(chunks=chunk_non_aligned)
        TS_chunk = [chk for chk in doc.chunks if re.search("TimeSIFT", chk.label) is not None]
        doc.remove(TS_chunk)
        merged_chunk = doc.chunks[-1]
        merged_chunk.label = "TimeSIFT"
        doc.remove(chk_sel)
        for cam in [i for i in merged_chunk.cameras if re.search(chk_sel.label,i.label) is not None]:
            cam.transform = None
    print("Temps écoulé pour la fusion : ", time.time() - start_time) 

def align_TimeSIFT_chunk(doc = scan.Document()):
    start_time = time.time()
    TS_chunk=[chk for chk in doc.chunks if re.search("TimeSIFT",chk.label) is not None][0]
    TS_chunk.matchPhotos(downscale=1, generic_preselection=True, reference_preselection=True,
                      reference_preselection_mode=scan.ReferencePreselectionSource, keypoint_limit=100000,
                      tiepoint_limit=10000,keep_keypoints=False)
    # boucle pour ré-aligner les photos non-alignées
    nb_aligned_before = 0
    nb_aligned_after = 100
    while nb_aligned_after != nb_aligned_before:
        nb_aligned_before = len([i for i in TS_chunk.cameras if i.transform])
        TS_chunk.matchPhotos(downscale=1, generic_preselection=True, reference_preselection=True,
                          reference_preselection_mode=scan.ReferencePreselectionSource, keypoint_limit=200000,
                          tiepoint_limit=20000,keep_keypoints=False)
        TS_chunk.alignCameras()
        nb_aligned_after = len([i for i in TS_chunk.cameras if i.transform])
    aligned_cameras = [i for i in TS_chunk.cameras if i.transform]
    print("------",str(len(aligned_cameras)),"cameras aligned out of",len(TS_chunk.cameras),"----")
    TS_chunk.resetRegion()
    TS_chunk.optimizeCameras()
    TS_chunk.updateTransform()
    print("Temps écoulé pour l'alignement : ", time.time() - start_time)

def split_TimeSIFT_chunk(doc = scan.Document()):
    TS_chunk = [chk for chk in doc.chunks if (re.search("TimeSIFT", chk.label) is not None)][0]
    #TS_chunk_names=np.unique([cam.label.split("_EPOCH_")[0] for cam in TS_chunk.cameras])
    TS_chunk_names=np.unique([re.search(r'\d+', cam.label).group() for cam in TS_chunk.cameras])
    print(TS_chunk_names)
    for chk_name in TS_chunk_names:
        if len([chk for chk in doc.chunks if re.search(chk_name,chk.label) is not None])==0:
            NewChunk = TS_chunk.copy()
            NewChunk.label=chk_name
            #pattern = str(chk_name) + "_EPOCH_"
            pattern = str(chk_name)
            t = [pattern in cam.label for cam in NewChunk.cameras]
            list_cameras = [NewChunk.cameras[i] for i, x in enumerate(t) if not x]
            NewChunk.remove(list_cameras)


def merge_chunk_with_same_date(doc = scan.Document()):
    """
    Not useful anymore. Merging by date now done in split_TIMESift_chunk
    """
    dates = []
    Non_ts_chunks = [chk for chk in doc.chunks if re.search("TimeSIFT",chk.label) is None]
    for chk in Non_ts_chunks:
        chunk_date = chk.label[:8]
        if chunk_date not in dates :
            dates.append(chunk_date)
    print("Dates : ", dates)
    for date in dates:
        chunks_to_merge = [chk.key for chk in Non_ts_chunks if chk.label[:8]==date]
        doc.mergeChunks(chunks=chunks_to_merge)
        merged_chunk = doc.chunks[-1]
        merged_chunk.label = date
        for chk in [chk for chk in Non_ts_chunks if chk.label[:8]==date]:  
            doc.remove(chk)

    

def process_splited_TimeSIFT_chunks_one_by_one(doc = scan.Document(), pathDIR=None, out_dir_ortho = None, out_dir_DEM = None, site_name=None, resol_ref = None, crs = None):
    """
    Generate depth map, dense cloud, DEM and orthomosaic for one image. Always saves orthomosaic and saves DEM if specified
    """
    TS_chunks = [chk for chk in doc.chunks if (re.search("TimeSIFT", chk.label) is None)]
    TS_chunks = [chk for chk in TS_chunks if chk.enabled]
    TS_chunk_names=[chk.label for chk in TS_chunks]
    for chk in TS_chunks:
        start_time = time.time()
        NewChunk=chk
        NewChunk.buildDepthMaps(downscale=2, filter_mode=scan.AggressiveFiltering)
        t_depth_maps = time.time()
        print("Time to build depth map : ", t_depth_maps - start_time)
        NewChunk.buildPointCloud(point_colors=True)
        t_dense_cloud = time.time()
        print("Time to build dense cloud : ", t_dense_cloud - t_depth_maps)
        NewChunk.buildDem(source_data=scan.PointCloudData,resolution=resol_ref)
        t_DEM = time.time()
        print("Time to build DEM : ", t_DEM - t_dense_cloud)
        NewChunk.buildOrthomosaic(surface_data=scan.ElevationData,resolution=resol_ref)
        t_ortho = time.time()
        print("Time to build ortho : ", t_ortho - t_DEM)
        print("Total process time for the image : ", t_ortho - start_time)
        proj = scan.OrthoProjection()
        proj.type=scan.OrthoProjection.Type.Planar
        proj.crs=scan.CoordinateSystem(crs)
        img_compress=scan.ImageCompression
        img_compress.tiff_compression=scan.ImageCompression.TiffCompressionLZW
        img_compress.tiff_big = True
        #doc.save(os.path.join(out_dir_ortho, '_temp_.psx'))
        try :
            NewChunk.exportRaster(os.path.join(out_dir_ortho, f"{str(NewChunk.label)}{site_name}_ORTHO.tif"),source_data=scan.OrthomosaicData, image_format=scan.ImageFormatTIFF,
                                projection=proj, resolution=resol_ref,clip_to_boundary=True,save_alpha=False, split_in_blocks = False)
        except:
            NewChunk.exportRaster(os.path.join(out_dir_ortho, f"{str(NewChunk.label)}{site_name}_ORTHO.tif"),source_data=scan.OrthomosaicData, image_format=scan.ImageFormatTIFF,
                                    projection=proj, resolution=resol_ref,clip_to_boundary=True,save_alpha=False, split_in_blocks = True, block_width=10000, block_height=10000)

        if out_dir_DEM is not None:
            NewChunk.exportRaster(os.path.join(out_dir_DEM, f"{str(NewChunk.label)}{site_name}_DEM.tif"),source_data=scan.ElevationData, image_format=scan.ImageFormatTIFF,
                                projection=proj, resolution=resol_ref,clip_to_boundary=True, save_alpha=False)


#TODO : confirm whether or not DEMs and project are to be saved by default
def Time_SIFT_process(pathDIR,
                      out_dir_ortho, 
                      out_dir_DEM=None,      #""
                      out_dir_project=None,    #""
                      data_type="RGB",
                      resol_ref=0.05, 
                      crs="EPSG::32622", 
                      site_name = "",
                      calibrate_col = True,
                      sun_sensor = False,
                      doc = scan.Document(),
                      ):
    
    assert data_type in ["RGB", "MS"]
    #for file naming purposes
    if site_name != "":
       site_name = "_" + site_name

    calibrate_col = str2bool(calibrate_col)

    if not os.path.exists(out_dir_ortho):
        os.mkdir(out_dir_ortho)
    
    if out_dir_DEM is not None:  
        if out_dir_DEM == "" :
           out_dir_DEM = os.path.join(os.path.dirname(out_dir_ortho), "DEM")
        if not os.path.exists(out_dir_DEM):
           os.mkdir(out_dir_DEM)

    #TODO : store all times into log file, or add progress bars
    start_time = time.time()
    if data_type == "RGB" :
        add_all_chunks(doc, pathDIR = pathDIR)
        merge_chunk_TimeSIFT(doc)
        t_add_data = time.time()
        print("Temps écoulé pour le chargement et la fusion des photos : ", t_add_data - start_time)

    if sun_sensor and data_type=='MS':
        TS_chunk = [chk for chk in doc.chunks if (re.search("TimeSIFT", chk.label) is not None)][0]
        TS_chunk.calibrateReflectance(use_sun_sensor=True)

    elif data_type == "MS" :
        add_all_MS_photos(doc, pathDIR = pathDIR)
        t_add_data = time.time()
        print("Temps écoulé pour le chargement des photos : ", t_add_data - start_time)
        
    align_TimeSIFT_chunk(doc)
    t_align = time.time()
    print("Temps écoulé pour l'alignement : ", t_align - t_add_data)
    
    #The project needs to be saved before building DEMs and orthomosaics
    doc.save(os.path.join(out_dir_ortho, '_temp_.psx'))

    #Color calibration
    if calibrate_col and (data_type=='RGB' or 'MS'):
        TS_chunk = [chk for chk in doc.chunks if (re.search("TimeSIFT", chk.label) is not None)][0]
        TS_chunk.calibrateColors(scan.TiePointsData, white_balance=True)
    


    doc.save(os.path.join(out_dir_ortho, '_temp_.psx'))
    
    split_TimeSIFT_chunk(doc)
    #merge_chunk_with_same_date(doc)
    t_split = time.time()
    #print("Temps écoulé pour la division et regroupement par date : ", t_split - t_align)
    process_splited_TimeSIFT_chunks_one_by_one(doc, pathDIR = pathDIR, out_dir_ortho = out_dir_ortho, out_dir_DEM = out_dir_DEM, site_name = site_name, resol_ref = resol_ref, crs = crs)
    print("Temps écoulé pour le process final : ", time.time() - t_split)
    print("Temps écoulé pour la pipeline complète : ", time.time() - start_time)
    doc.save(os.path.join(out_dir_ortho, '_temp_.psx'))
    
    if out_dir_project is not None :
        if out_dir_project == "" :
            out_dir_project = out_dir_ortho
        if not os.path.exists(out_dir_project):
            os.mkdir(out_dir_project)
        doc.save(os.path.join(out_dir_project, f"Metashape_Project_{site_name}.psx"))
        
    #os.remove(os.path.join(out_dir_ortho, '_temp_.psx'))
    #shutil.rmtree(os.path.join(out_dir_ortho, 'temp.files'))
    
"""
try:

    #doc = scan.Document()
    #complete_process_RGB(doc, pathDIR="Y:\RGB", resol_ref=0.05, crs="EPSG::32622")
    print("args : ", args)
    Time_SIFT_process(pathDIR="Y:/RGB/RGB", out_dir_ortho= "Z:/shared/PhenOBS/Paracou/Metashape/RGB_Broad_Mosaics/4D_050324/ORTHO", out_dir_DEM= "Z:/shared/PhenOBS/Paracou/Metashape/RGB_Broad_Mosaics/4D_050324/DEM", 
                         data_type="RGB", resol_ref=0.05, crs="EPSG::32622")
    #complete_process_and_save(doc, pathDIR = args.path_in, resol_ref = args.resol_ref, crs = args.crs, doc)
    #complete_process_MS(doc, pathDIR="Y:\MS\P4M\batch_test", resol_ref=0.05, crs="EPSG::32622")
    #add_all_chunks(doc, "Y:\RGB")
    
    #doc.save(args.path_out)
    #doc.save("Z:/users/GaelleViennois/Phenobs_recalage/test.psx")
    

except Exception as e:
    print(f"An error occurred: {e}")
"""
