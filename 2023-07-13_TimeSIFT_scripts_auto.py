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

#scan.License().activate('E2TL6-7YYPS-9GUPV-YEGB8-NUPUE')

parser = argparse.ArgumentParser()
parser.add_argument('--crs')
parser.add_argument('--pathDIR')
parser.add_argument('--out_dir_ortho')
parser.add_argument('--out_dir_dem')
parser.add_argument('--resol_ref')
parser.add_argument('--data_type')
args = parser.parse_args()

def add_TimeSIFT_chunk(doc = scan.Document(), epoch_name="2019-05-23"):
    if epoch_name !="no":
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
    if pathDIR is None :
        pathDIR=scan.app.getExistingDirectory("select the folder where the different epochs folders are located")
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
        # Pour chaque epoch, on créé un groupe de cameras
        chunk.addCameraGroup()
        cam_group = chunk.camera_groups[-1]
        cam_group.label = epoch_name
        
        dirName= os.path.join(pathDIR, epoch_name)
        listOfFiles = list()
        for (dirpath, dirnames, filenames) in os.walk(dirName):
            listOfFiles += [os.path.join(dirpath, file) for file in filenames]
        
        # on ajoute toutes les photos liées à une même epoch dans le groupe de cameras correspondant
        chunk.addPhotos(listOfFiles, layout = scan.MultiplaneLayout, group = i)  
        for cam in chunk.cameras[L:]:
            cam.label = (str(epoch_name) + "_EPOCH_" + cam.label)  
            
def delete_GPS_exif_chunk(doc = scan.Document()):
    for chk in doc.chunks:
      if chk.enabled != True: 
        for cam in chk.cameras:
          cam.reference.location=None


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
    dates = []
    Non_ts_chunks = [chk for chk in doc.chunks if re.search("TimeSIFT",chk.label) is None]
    for chk in Non_ts_chunks:
      chunk_date = chk.label[:8]
      if chunk_date not in dates :
         dates.append(chunk_date)
    print("Dates : ", dates)
    for date in dates:
      chunks_to_merge = [chk.key for chk in Non_ts_chunks if chk.label[:8]==date]
      chk_sel = chunks_to_merge[0]
      doc.mergeChunks(chunks=chunks_to_merge)
      merged_chunk = doc.chunks[-1]
      merged_chunk.label = date
      for chk in [chk for chk in Non_ts_chunks if chk.label[:8]==date]:  
        doc.remove(chk)

    

def process_splited_TimeSIFT_chunks_one_by_one(doc = scan.Document(), pathDIR=None, out_dir_ortho = None, out_dir_DEM = None, counter = 0, site_name=None, resol_ref = None, crs = None):
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
    proj.crs=scan.CoordinateSystem(crs)#EPSG::32622
    img_compress=scan.ImageCompression
    img_compress.tiff_compression=scan.ImageCompression.TiffCompressionLZW
    img_compress.tiff_big = True
    doc.save(os.path.join(out_dir_ortho, 'temp.psx'))
    try :
        NewChunk.exportRaster(os.path.join(out_dir_ortho, f"ID{counter}{site_name}_{str(NewChunk.label)}_ORTHO.tif"),source_data=scan.OrthomosaicData, image_format=scan.ImageFormatTIFF,
                              projection=proj, resolution=resol_ref,clip_to_boundary=True,save_alpha=False, split_in_blocks = False)
    except:
        NewChunk.exportRaster(os.path.join(out_dir_ortho, f"ID{counter}{site_name}_{str(NewChunk.label)}_ORTHO.tif"),source_data=scan.OrthomosaicData, image_format=scan.ImageFormatTIFF,
                                  projection=proj, resolution=resol_ref,clip_to_boundary=True,save_alpha=False, split_in_blocks = True, block_width=10000, block_height=10000)
    NewChunk.exportRaster(os.path.join(out_dir_ortho, f"ID{counter}{site_name}_{str(NewChunk.label)}_ORTHO.tif"),source_data=scan.OrthomosaicData, image_format=scan.ImageFormatTIFF,
                          projection=proj, resolution=resol_ref,clip_to_boundary=True, save_alpha=False)
    if out_dir_DEM is not None:
        NewChunk.exportRaster(os.path.join(out_dir_DEM, f"ID{counter}{site_name}_{str(NewChunk.label)}_DEM.tif"),source_data=scan.ElevationData, image_format=scan.ImageFormatTIFF,
                              projection=proj, resolution=resol_ref,clip_to_boundary=True, save_alpha=False)


def add_all_chunks(doc = scan.Document(), pathDIR=None):
  print(pathDIR)
  os.chdir(pathDIR)
  epochs = os.listdir(pathDIR)
  # we select only the non-empty subfolders 
  epochs = [ep for ep in epochs if not os.path.isfile(ep) and os.listdir(ep) and ep!="20221108_P15" and ep!="20211018_Sud" and ep!="20210208_Nord"]
  print(epochs)
  # We remove all existing chunks and add them one by one
  for chk in doc.chunks:
    doc.remove(chk)
  for ep in epochs: 
    add_TimeSIFT_chunk(doc, epoch_name=ep)
  
def Time_SIFT_process(pathDIR,
                      dir_out, 
                      data_type="RGB",
                      resol_ref=0.05, 
                      crs="EPSG::32622", 
                      site_name = "",
                      save_MS_project = False,
                      calibrate_col = True,
                      doc = scan.Document(),
                      ):
    assert data_type in ["RGB", "MS"]

    #for file naming purposes
    if site_name != "":
       site_name = "_" + site_name
        
    # Initialize counter or load it if it exists
    counter_dir = os.path.join(dir_out, '.metadata')
    if not os.path.exists(counter_dir):
        os.mkdir(counter_dir)
    counter_file = os.path.join(counter_dir, '.project_counter.txt')
    counter = 0
    if os.path.exists(counter_file):
        with open(counter_file, 'r') as file:
            counter = int(file.read())
    # Increment the counter
    counter += 1
    # Create all required directories if they don't already exist   
    if not os.path.exists(dir_out):
        os.mkdir(dir_out)
    if not os.path.exists(os.path.join(dir_out, f"Project{counter}{site_name}")):
        os.mkdir(os.path.join(dir_out, f"Project{counter}{site_name}"))
    out_dir_ortho = os.path.join(dir_out, f"Project{counter}{site_name}", "ORTHO")
    out_dir_DEM = os.path.join(dir_out, f"Project{counter}{site_name}", "DEM")
    if not os.path.exists(out_dir_ortho):
        os.mkdir(out_dir_ortho)
    if not os.path.exists(out_dir_DEM):
        os.mkdir(out_dir_DEM)

    # Update the counter file
    with open(counter_file, 'w') as file:
        file.write(str(counter))

    start_time = time.time()
    if data_type == "RGB" :
        add_all_chunks(doc, pathDIR = pathDIR)
        merge_chunk_TimeSIFT(doc)
        t_add_data = time.time()
        print("Temps écoulé pour le chargement et la fusion des photos : ", t_add_data - start_time)
        #complete_process_RGB(doc=doc, pathDIR=pathDIR, resol_ref=resol_ref, crs=crs)
    elif data_type == "MS" :
        add_all_MS_photos(doc, pathDIR = pathDIR)
        t_add_data = time.time()
        print("Temps écoulé pour le chargement des photos : ", t_add_data - start_time)
        #complete_process_MS(doc=doc, pathDIR=pathDIR, resol_ref=resol_ref, crs=crs)
    

    
    align_TimeSIFT_chunk(doc)
    t_align = time.time()
    print("Temps écoulé pour l'alignement : ", t_align - t_add_data)
    
   #Color calibration
    if calibrate_col and data_type=='RGB':
        TS_chunk = [chk for chk in doc.chunks if (re.search("TimeSIFT", chk.label) is not None)][0]
        TS_chunk.calibrateColors(scan.PointCloudData, white_balance=False)
    
    #The project needs to be saved before building DEMs and orthomosaics
    doc.save(os.path.join(dir_out, f"Project{counter}{site_name}", 'temp.psx'))
    
    split_TimeSIFT_chunk(doc)
    merge_chunk_with_same_date(doc)
    t_split = time.time()
    print("Temps écoulé pour la division et regroupement par date : ", t_split - t_align)
    process_splited_TimeSIFT_chunks_one_by_one(doc, pathDIR = pathDIR, out_dir_ortho = out_dir_ortho, out_dir_DEM = out_dir_DEM, counter = counter, site_name = site_name, resol_ref = resol_ref, crs = crs)
    print("Temps écoulé pour le process final : ", time.time() - t_split)
    print("Temps écoulé pour la pipeline complète : ", time.time() - start_time)
    
    os.remove(os.path.join(dir_out, f"Project{counter}{site_name}", 'temp.psx'))
    #shutil.rmtree(os.path.join(dir_out, f"Project{counter}{site_name}", "temp.files"))
    
    if save_MS_project :
        doc.save(os.path.join(os.path.join(dir_out, f"Project{counter}{site_name}"), f"Project_{counter}_{site_name}.psx"))

def temp_func_2(doc = scan.Document(), pathDIR=None, out_dir_ortho = None, out_dir_DEM = None, counter = 7, site_name=None, resol_ref = 0.05, crs = None):     #resume progress if project saved right after alignement
    if out_dir_DEM is not None:  
        if out_dir_DEM == "" :
           out_dir_DEM = os.path.join(os.path.dirname(out_dir_ortho), "DEM")
        if not os.path.exists(out_dir_DEM):
           os.mkdir(out_dir_DEM)
    doc.open(os.path.join(out_dir_ortho, 'temp.psx'))
    split_TimeSIFT_chunk(doc)
    #merge_chunk_with_same_date(doc)
    t_split = time.time()
    doc.save(os.path.join(out_dir_ortho, 'temp.psx'))
    #print("Temps écoulé pour la division et regroupement par date : ", t_split - t_align)
    process_splited_TimeSIFT_chunks_one_by_one(doc, pathDIR = pathDIR, out_dir_ortho = out_dir_ortho, out_dir_DEM = out_dir_DEM, counter = counter, site_name = site_name, resol_ref = resol_ref, crs = crs)
    print("Temps écoulé pour le process final : ", time.time() - t_split)
    

def temp_func_3(doc = scan.Document(), pathDIR=None, out_dir_ortho = None, out_dir_DEM = None, counter = 5, site_name=None, resol_ref = 0.05, crs = None):   #test split MS
    if not os.path.exists(out_dir_ortho):
        os.mkdir(out_dir_ortho)
    if out_dir_DEM is not None:  
        if out_dir_DEM == "" :
           out_dir_DEM = os.path.join(os.path.dirname(out_dir_ortho), "DEM")
        if not os.path.exists(out_dir_DEM):
           os.mkdir(out_dir_DEM)
    doc.open(pathDIR)
    TS_chunks = [chk for chk in doc.chunks if (re.search("TimeSIFT", chk.label) is None)][:-1]
    TS_chunks = [chk for chk in TS_chunks if chk.enabled]
    TS_chunk_names=[chk.label for chk in TS_chunks]
    print(TS_chunk_names)
    for chk in TS_chunks:
        NewChunk=chk
        proj = scan.OrthoProjection()
        proj.type=scan.OrthoProjection.Type.Planar
        proj.crs=scan.CoordinateSystem(crs)#EPSG::32622
        img_compress=scan.ImageCompression
        img_compress.tiff_compression=scan.ImageCompression.TiffCompressionLZW
        img_compress.tiff_big = True
        #doc.save(os.path.join(out_dir_ortho, 'temp.psx'))
        NewChunk.exportRaster(os.path.join(out_dir_ortho, f"ID{counter}{site_name}_{str(NewChunk.label)}_ORTHO.tif"),source_data=scan.OrthomosaicData, image_format=scan.ImageFormatTIFF,
                          projection=proj, resolution=resol_ref,clip_to_boundary=True,save_alpha=False, split_in_blocks = True, block_width=12000, block_height=12000)
        

def Time_SIFT_process_v2(pathDIR,
                         out_dir_ortho, 
                         out_dir_DEM=None,      #""
                         out_dir_project=None,    #""
                         data_type="RGB",
                         resol_ref=0.05, 
                         crs="EPSG::32622", 
                         site_name = "",
                         calibrate_col = True,
                         doc = scan.Document(),
                          ):
    assert data_type in ["RGB", "MS"]
    #for file naming purposes
    if site_name != "":
       site_name = "_" + site_name

    if not os.path.exists(out_dir_ortho):
        os.mkdir(out_dir_ortho)
    
    if out_dir_DEM is not None:  
        if out_dir_DEM == "" :
           out_dir_DEM = os.path.join(os.path.dirname(out_dir_ortho), "DEM")
        if not os.path.exists(out_dir_DEM):
           os.mkdir(out_dir_DEM)
      
    # Initialize counter or load it if it exists
    counter_dir = os.path.join(out_dir_ortho, '.metadata')
    if not os.path.exists(counter_dir):
        os.mkdir(counter_dir)
    
    counter_file = os.path.join(counter_dir, '.project_counter.txt') 
    counter = 0
    if os.path.exists(counter_file):
        with open(counter_file, 'r') as file:
            counter = int(file.read())
    # Increment the counter
    counter += 1
    # Update the counter file
    with open(counter_file, 'w') as file:
        file.write(str(counter))
    start_time = time.time()
    if data_type == "RGB" :
        add_all_chunks(doc, pathDIR = pathDIR)
        merge_chunk_TimeSIFT(doc)
        t_add_data = time.time()
        print("Temps écoulé pour le chargement et la fusion des photos : ", t_add_data - start_time)
        
    elif data_type == "MS" :
        add_all_MS_photos(doc, pathDIR = pathDIR)
        t_add_data = time.time()
        print("Temps écoulé pour le chargement des photos : ", t_add_data - start_time)
        
    align_TimeSIFT_chunk(doc)
    t_align = time.time()
    print("Temps écoulé pour l'alignement : ", t_align - t_add_data)
    
    #The project needs to be saved before building DEMs and orthomosaics
    doc.save(os.path.join(out_dir_ortho, 'temp.psx'))

    #Color calibration
    if calibrate_col and data_type=='RGB':
        TS_chunk = [chk for chk in doc.chunks if (re.search("TimeSIFT", chk.label) is not None)][0]
        TS_chunk.calibrateColors(scan.TiePointsData, white_balance=True)

    doc.save(os.path.join(out_dir_ortho, 'temp.psx'))
    
    split_TimeSIFT_chunk(doc)
    #merge_chunk_with_same_date(doc)
    t_split = time.time()
    #print("Temps écoulé pour la division et regroupement par date : ", t_split - t_align)
    process_splited_TimeSIFT_chunks_one_by_one(doc, pathDIR = pathDIR, out_dir_ortho = out_dir_ortho, out_dir_DEM = out_dir_DEM, counter = counter, site_name = site_name, resol_ref = resol_ref, crs = crs)
    print("Temps écoulé pour le process final : ", time.time() - t_split)
    print("Temps écoulé pour la pipeline complète : ", time.time() - start_time)
    doc.save(os.path.join(out_dir_ortho, 'temp.psx'))
    
    if out_dir_project is not None :
        if out_dir_project == "" :
            out_dir_project = out_dir_ortho
        if not os.path.exists(out_dir_project):
            os.mkdir(out_dir_project)
        doc.save(os.path.join(os.path.join(out_dir_project, f"Project_{counter}_{site_name}.psx")))
        
    #os.remove(os.path.join(out_dir_ortho, 'temp.psx'))
    #shutil.rmtree(os.path.join(out_dir_ortho, 'temp.files'))
    

try:

    #doc = scan.Document()
    #complete_process_RGB(doc, pathDIR="Y:\RGB", resol_ref=0.05, crs="EPSG::32622")
    print("args : ", args)
    Time_SIFT_process_v2(pathDIR="Y:/RGB/RGB", out_dir_ortho= "Z:/shared/PhenOBS/Paracou/Metashape/RGB_Broad_Mosaics/4D_050324/ORTHO", out_dir_DEM= "Z:/shared/PhenOBS/Paracou/Metashape/RGB_Broad_Mosaics/4D_050324/DEM", 
                         data_type="RGB", resol_ref=0.05, crs="EPSG::32622")
    #complete_process_and_save(doc, pathDIR = args.path_in, resol_ref = args.resol_ref, crs = args.crs, doc)
    #complete_process_MS(doc, pathDIR="Y:\MS\P4M\batch_test", resol_ref=0.05, crs="EPSG::32622")
    #add_all_chunks(doc, "Y:\RGB")
    
    #doc.save(args.path_out)
    #doc.save("Z:/users/GaelleViennois/Phenobs_recalage/test.psx")
    

except Exception as e:
    print(f"An error occurred: {e}")

