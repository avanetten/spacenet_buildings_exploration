#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: avanetten
"""

'''
Transform SpaceNet geojson buidling labels data into raster masks.
Download data via:
aws s3api get-object --bucket spacenet-dataset \
--key AOI_1_Rio/processedData/processedBuildingLabels.tar.gz \
--request-payer requester processedBuildingLabels.tar.gz

Download spacenet utilities from:
https://github.com/SpaceNetChallenge/utilities/tree/master/python/spaceNet

For further details, see:
https://medium.com/the-downlinq/getting-started-with-spacenet-data-827fd2ec9f53
'''

from __future__ import print_function
from matplotlib.collections import PatchCollection
from osgeo import gdal, osr, ogr, gdalnumeric
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np
import shutil
import json
import glob
import sys
import os


####################
# EDIT THESE PATHS
# download spacenet utilities from:
#   https://github.com/SpaceNetChallenge/utilities/tree/master/python/spaceNet 
path_to_spacenet_utils = '/path_to_spacenet_utils'
# Set data dir
spacenet_data_dir = '/path_to_spacenet_data'
# This is the directory where this script is located
spacenet_explore_dir = os.path.dirname(os.path.realpath(__file__))
# exclore N images in 3band data
N_ims = 15
####################

# import packages
sys.path.extend([path_to_spacenet_utils])
from spaceNetUtilities import geoTools as gT

       
###############################################################################    
def geojson_to_pixel_arr(raster_file, geojson_file, pixel_ints=True,
                       verbose=False):
    '''
    Tranform geojson file into array of points in pixel (and latlon) coords
    pixel_ints = 1 sets pixel coords as integers
    '''
    
    # load geojson file
    with open(geojson_file) as f:
        geojson_data = json.load(f)

    # load raster file and get geo transforms
    src_raster = gdal.Open(raster_file)
    targetsr = osr.SpatialReference()
    targetsr.ImportFromWkt(src_raster.GetProjectionRef())
        
    geom_transform = src_raster.GetGeoTransform()
    if verbose:
        print ("geom_transform:", geom_transform)
    
    # get latlon coords
    latlons = []
    types = []
    for feature in geojson_data['features']:
        coords_tmp = feature['geometry']['coordinates'][0]
        type_tmp = feature['geometry']['type']
        if verbose: 
            print ("features:", feature.keys())
            print ("geometry:features:", feature['geometry'].keys())

            #print "feature['geometry']['coordinates'][0]", z
        latlons.append(coords_tmp)
        types.append(type_tmp)
        #print feature['geometry']['type']
    
    # convert latlons to pixel coords
    pixel_coords = []
    latlon_coords = []
    for i, (poly_type, poly0) in enumerate(zip(types, latlons)):
        
        if poly_type.upper() == 'MULTIPOLYGON':
            #print "oops, multipolygon"
            for poly in poly0:
                poly=np.array(poly)
                if verbose:
                    print ("poly.shape:", poly.shape)
                    
                # account for nested arrays
                if len(poly.shape) == 3 and poly.shape[0] == 1:
                    poly = poly[0]
                    
                poly_list_pix = []
                poly_list_latlon = []
                if verbose: 
                    print ("poly", poly)
                for coord in poly:
                    if verbose: 
                        print ("coord:", coord)
                    lon, lat, z = coord 
                    px, py = gT.latlon2pixel(lat, lon, input_raster=src_raster, 
                                         targetsr=targetsr, 
                                         geom_transform=geom_transform)
                    poly_list_pix.append([px, py])
                    if verbose:
                        print ("px, py", px, py)
                    poly_list_latlon.append([lat, lon])
                
                if pixel_ints:
                    ptmp = np.rint(poly_list_pix).astype(int)
                else:
                    ptmp = poly_list_pix
                pixel_coords.append(ptmp)
                latlon_coords.append(poly_list_latlon)            

        elif poly_type.upper() == 'POLYGON':
            poly=np.array(poly0)
            if verbose:
                print ("poly.shape:", poly.shape)
                
            # account for nested arrays
            if len(poly.shape) == 3 and poly.shape[0] == 1:
                poly = poly[0]
                
            poly_list_pix = []
            poly_list_latlon = []
            if verbose: 
                print ("poly", poly)
            for coord in poly:
                if verbose: 
                    print ("coord:", coord)
                lon, lat, z = coord 
                px, py = gT.latlon2pixel(lat, lon, input_raster=src_raster, 
                                     targetsr=targetsr, 
                                     geom_transform=geom_transform)
                poly_list_pix.append([px, py])
                if verbose:
                    print ("px, py", px, py)
                poly_list_latlon.append([lat, lon])
            
            if pixel_ints:
                ptmp = np.rint(poly_list_pix).astype(int)
            else:
                ptmp = poly_list_pix
            pixel_coords.append(ptmp)
            latlon_coords.append(poly_list_latlon)
            
        elif poly_type.upper() == 'POINT':
            print ("Skipping shape type: POINT in geojson_to_pixel_arr()")
            continue
        
        else:
            print ("Unknown shape type:", poly_type, " in geojson_to_pixel_arr()")
            return
            
    return pixel_coords, latlon_coords

###############################################################################
def create_building_mask(rasterSrc, vectorSrc, npDistFileName='', 
                            noDataValue=0, burn_values=1):

    '''
    Create building mask for rasterSrc,
    Similar to labeltools/createNPPixArray() in spacenet utilities
    '''
    
    ## open source vector file that truth data
    source_ds = ogr.Open(vectorSrc)
    source_layer = source_ds.GetLayer()

    ## extract data from src Raster File to be emulated
    ## open raster file that is to be emulated
    srcRas_ds = gdal.Open(rasterSrc)
    cols = srcRas_ds.RasterXSize
    rows = srcRas_ds.RasterYSize

    ## create First raster memory layer, units are pixels
    # Change output to geotiff instead of memory 
    memdrv = gdal.GetDriverByName('GTiff') 
    dst_ds = memdrv.Create(npDistFileName, cols, rows, 1, gdal.GDT_Byte, 
                           options=['COMPRESS=LZW'])
    dst_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
    dst_ds.SetProjection(srcRas_ds.GetProjection())
    band = dst_ds.GetRasterBand(1)
    band.SetNoDataValue(noDataValue)    
    gdal.RasterizeLayer(dst_ds, [1], source_layer, burn_values=[burn_values])
    dst_ds = 0
    
    return 

###############################################################################
def create_dist_map(rasterSrc, vectorSrc, npDistFileName='', 
                           noDataValue=0, burn_values=1, 
                           dist_mult=1, vmax_dist=64):

    '''
    Create building signed distance transform from Yuan 2016 
    (https://arxiv.org/pdf/1602.06564v1.pdf).
    vmax_dist: absolute value of maximum distance (meters) from building edge
    Adapged from createNPPixArray in labeltools
    '''
    
    ## open source vector file that truth data
    source_ds = ogr.Open(vectorSrc)
    source_layer = source_ds.GetLayer()

    ## extract data from src Raster File to be emulated
    ## open raster file that is to be emulated
    srcRas_ds = gdal.Open(rasterSrc)
    cols = srcRas_ds.RasterXSize
    rows = srcRas_ds.RasterYSize

    geoTrans, poly, ulX, ulY, lrX, lrY = gT.getRasterExtent(srcRas_ds)
    transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs \
                                                = gT.createUTMTransform(poly)
    line = ogr.Geometry(ogr.wkbLineString)
    line.AddPoint(geoTrans[0], geoTrans[3])
    line.AddPoint(geoTrans[0]+geoTrans[1], geoTrans[3])

    line.Transform(transform_WGS84_To_UTM)
    metersIndex = line.Length()

    memdrv = gdal.GetDriverByName('MEM')
    dst_ds = memdrv.Create('', cols, rows, 1, gdal.GDT_Byte)
    dst_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
    dst_ds.SetProjection(srcRas_ds.GetProjection())
    band = dst_ds.GetRasterBand(1)
    band.SetNoDataValue(noDataValue)

    gdal.RasterizeLayer(dst_ds, [1], source_layer, burn_values=[burn_values])
    srcBand = dst_ds.GetRasterBand(1)

    memdrv2 = gdal.GetDriverByName('MEM')
    prox_ds = memdrv2.Create('', cols, rows, 1, gdal.GDT_Int16)
    prox_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
    prox_ds.SetProjection(srcRas_ds.GetProjection())
    proxBand = prox_ds.GetRasterBand(1)
    proxBand.SetNoDataValue(noDataValue)

    opt_string = 'NODATA='+str(noDataValue)
    options = [opt_string]

    gdal.ComputeProximity(srcBand, proxBand, options)

    memdrv3 = gdal.GetDriverByName('MEM')
    proxIn_ds = memdrv3.Create('', cols, rows, 1, gdal.GDT_Int16)
    proxIn_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
    proxIn_ds.SetProjection(srcRas_ds.GetProjection())
    proxInBand = proxIn_ds.GetRasterBand(1)
    proxInBand.SetNoDataValue(noDataValue)
    opt_string2 = 'VALUES='+str(noDataValue)
    options = [opt_string, opt_string2]
    #options = ['NODATA=0', 'VALUES=0']

    gdal.ComputeProximity(srcBand, proxInBand, options)

    proxIn = gdalnumeric.BandReadAsArray(proxInBand)
    proxOut = gdalnumeric.BandReadAsArray(proxBand)

    proxTotal = proxIn.astype(float) - proxOut.astype(float)
    proxTotal = proxTotal*metersIndex
    proxTotal *= dist_mult

    # clip array
    proxTotal = np.clip(proxTotal, -1*vmax_dist, 1*vmax_dist)

    if npDistFileName != '':
        # save as numpy file since some values will be negative
        np.save(npDistFileName, proxTotal)
        #cv2.imwrite(npDistFileName, proxTotal)

    #return proxTotal
    return

###############################################################################
def plot_truth_coords(input_image, pixel_coords,   
                  figsize=(8,8), plot_name='',
                  add_title=False, poly_face_color='orange', 
                  poly_edge_color='red', poly_nofill_color='blue', cmap='bwr'):
    '''Plot ground truth coordinaates, pixel_coords should be a numpy array'''
    
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(2*figsize[0], figsize[1]))
    
    if add_title:
        suptitle = fig.suptitle(plot_name.split('/')[-1], fontsize='large')
    
    # create patches
    patches = []
    patches_nofill = []
    if len(pixel_coords) > 0:
        # get patches    
        for coord in pixel_coords:
            patches_nofill.append(Polygon(coord, facecolor=poly_nofill_color, 
                                          edgecolor=poly_edge_color, lw=3))
            patches.append(Polygon(coord, edgecolor=poly_edge_color, fill=True, 
                                   facecolor=poly_face_color))
        p0 = PatchCollection(patches, alpha=0.25, match_original=True)
        #p1 = PatchCollection(patches, alpha=0.75, match_original=True)
        p2 = PatchCollection(patches_nofill, alpha=0.75, match_original=True)
                   
    # ax0: raw image
    ax0.imshow(input_image)
    if len(patches) > 0:
        ax0.add_collection(p0)
    ax0.set_title('Input Image + Ground Truth Buildings') 
    
    # truth polygons
    zero_arr = np.zeros(input_image.shape[:2])
    # set background to white?
    #zero_arr[zero_arr == 0.0] = np.nan
    ax1.imshow(zero_arr, cmap=cmap)
    if len(patches) > 0:
        ax1.add_collection(p2)
    ax1.set_title('Ground Truth Building Polygons')
        
    #plt.axis('off')
    plt.tight_layout()
    if add_title:
        suptitle.set_y(0.95)
        fig.subplots_adjust(top=0.96)
    #plt.show()
 
    if len(plot_name) > 0:
        plt.savefig(plot_name)
    
    return patches, patches_nofill
    
                     
###############################################################################
def plot_building_mask(input_image, pixel_coords, mask_image,   
                  figsize=(8,8), plot_name='',
                  add_title=False, poly_face_color='orange', 
                  poly_edge_color='red', poly_nofill_color='blue', cmap='bwr'):


    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, 
                                        figsize=(3*figsize[0], figsize[1]))
    
    if add_title:
        suptitle = fig.suptitle(plot_name.split('/')[-1], fontsize='large')

    # create patches
    patches = []
    patches_nofill = []
    if len(pixel_coords) > 0:
        # get patches    
        for coord in pixel_coords:
            patches_nofill.append(Polygon(coord, facecolor=poly_nofill_color, 
                                          edgecolor=poly_edge_color, lw=3))
            patches.append(Polygon(coord, edgecolor=poly_edge_color, fill=True, 
                                   facecolor=poly_face_color))
        p0 = PatchCollection(patches, alpha=0.25, match_original=True)
        p1 = PatchCollection(patches_nofill, alpha=0.75, match_original=True)
        
    #if len(patches) > 0:
    #    p0 = PatchCollection(patches, alpha=0.25, match_original=True)
    #    #p1 = PatchCollection(patches, alpha=0.75, match_original=True)
    #    p1 = PatchCollection(patches_nofill, alpha=0.75, match_original=True)                   
 
    # ax0: raw image
    ax0.imshow(input_image)
    if len(patches) > 0:
        ax0.add_collection(p0)
    ax0.set_title('Input Image + Ground Truth Buildings') 
    
    # truth polygons
    zero_arr = np.zeros(input_image.shape[:2])
    # set background to white?
    #zero_arr[zero_arr == 0.0] = np.nan
    ax1.imshow(zero_arr, cmap=cmap)
    if len(patches) > 0:
        ax1.add_collection(p1)
    ax1.set_title('Ground Truth Building Polygons')
        
    # old method of truth, with mask
    ## ax0: raw image
    #ax0.imshow(input_image)
    ## ground truth
    ## set zeros to nan
    #palette = plt.cm.gray
    #palette.set_over('orange', 1.0)
    #z = mask_image.astype(float)
    #z[z==0] = np.nan
    #ax0.imshow(z, cmap=palette, alpha=0.25, 
    #        norm=matplotlib.colors.Normalize(vmin=0.5, vmax=0.9, clip=False))
    #ax0.set_title('Input Image + Ground Truth Buildings') 
   
    # mask
    ax2.imshow(mask_image, cmap=cmap)
    # truth polygons?
    #if len(patches) > 0:
    #    ax1.add_collection(p1)
    ax2.set_title('Ground Truth Building Mask')    
          
    #plt.axis('off')
    plt.tight_layout()
    if add_title:
        suptitle.set_y(0.95)
        fig.subplots_adjust(top=0.96)
    #plt.show()
 
    if len(plot_name) > 0:
        plt.savefig(plot_name)
    
    return

###############################################################################
def plot_dist_transform(input_image, pixel_coords, dist_image, 
                        figsize=(8,8), plot_name='', add_title=False, 
                        colorbar=True,
                        poly_face_color='orange', poly_edge_color='red', 
                        poly_nofill_color='blue', cmap='bwr'):
    '''Explore distance transform'''

    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, 
                                        figsize=(3*figsize[0], figsize[1]))

    #fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(2*figsize[0], figsize[1]))
    mind, maxd = np.round(np.min(dist_image),2), np.round(np.max(dist_image),2)
    
    if add_title:
        suptitle = fig.suptitle(plot_name.split('/')[-1], fontsize='large')

    # create patches
    patches = []
    patches_nofill = []
    if len(pixel_coords) > 0:
        # get patches    
        for coord in pixel_coords:
            patches_nofill.append(Polygon(coord, facecolor=poly_nofill_color, 
                                          edgecolor=poly_edge_color, lw=3))
            patches.append(Polygon(coord, edgecolor=poly_edge_color, fill=True, 
                                   facecolor=poly_face_color))
        p0 = PatchCollection(patches, alpha=0.25, match_original=True)
        p1 = PatchCollection(patches, alpha=0.75, match_original=True)
        #p2 = PatchCollection(patches_nofill, alpha=0.75, match_original=True)
        
    #if len(patches) > 0:
    #    p0 = PatchCollection(patches, alpha=0.25, match_original=True)
    #    p1 = PatchCollection(patches, alpha=0.75, match_original=True)
                   
 
    # ax0: raw image
    ax0.imshow(input_image)
    if len(patches) > 0:
        ax0.add_collection(p0)
    ax0.set_title('Input Image + Ground Truth Buildings') 
    
    ## truth polygons
    #zero_arr = np.zeros(input_image.shape[:2])
    ## set background to white?
    ##zero_arr[zero_arr == 0.0] = np.nan
    #ax1.imshow(zero_arr, cmap=cmap)
    #if len(patches) > 0:
    #    ax1.add_collection(p1)
    #ax1.set_title('Ground Truth Building Outlines')
    
    # transform
    cbar_pointer = ax1.imshow(dist_image)
    dist_suffix = " (min=" + str(mind) + ", max=" + str(maxd) + ")"
    ax1.set_title("Yuan 2016 Distance Transform" + dist_suffix)
    
    # overlay buildings on distance transform
    ax2.imshow(dist_image)
    # truth polygons
    if len(patches) > 0:
        ax2.add_collection(p1)
    # truth mask
    #ax2.imshow(z, cmap=palette, alpha=0.5, 
    #     norm=matplotlib.colors.Normalize(vmin=0.5, vmax=0.9, clip=False))
    ax2.set_title("Ground Truth Polygons Overlaid on Distance Transform")
    
    if colorbar:
        #from mpl_toolkits.axes_grid1 import make_axes_locatable
        #divider = make_axes_locatable(ax2)
        #cax = divider.append_axes('right', size='5%', pad=0.05)
        #fig.colorbar(cbar_pointer, cax=cax, orientation='vertical')
        left, bottom, width, height = [0.38, 0.85, 0.24, 0.03]
        cax = fig.add_axes([left, bottom, width, height])
        fig.colorbar(cbar_pointer, cax=cax, orientation='horizontal')

    #plt.axis('off')
    plt.tight_layout()
    if add_title:
        suptitle.set_y(0.95)
        fig.subplots_adjust(top=0.96)
    #plt.show()
 
    if len(plot_name) > 0:
        plt.savefig(plot_name)
    
    return
 
###############################################################################
def plot_all_transforms(input_image, pixel_coords, mask_image, dist_image, 
                        figsize=(8,8), plot_name='', add_global_title=False, 
                        colorbar=False, add_titles=False,
                        poly_face_color='orange', poly_edge_color='red', 
                        poly_nofill_color='blue', cmap='bwr'):
    '''Explore all transforms'''

    fig, (ax0, ax1, ax2, ax3) = plt.subplots(1, 4, 
                                        figsize=(4*figsize[0], figsize[1]))

    #fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(2*figsize[0], figsize[1]))
    
    if add_global_title:
        suptitle = fig.suptitle(plot_name.split('/')[-1], fontsize='large')

    # create patches
    patches = []
    patches_nofill = []
    if len(pixel_coords) > 0:
        # get patches    
        for coord in pixel_coords:
            patches_nofill.append(Polygon(coord, facecolor=poly_nofill_color, 
                                          edgecolor=poly_edge_color, lw=3))
            patches.append(Polygon(coord, edgecolor=poly_edge_color, fill=True, 
                                   facecolor=poly_face_color))
        p0 = PatchCollection(patches, alpha=0.25, match_original=True)
        #p1 = PatchCollection(patches, alpha=0.75, match_original=True)
        p2 = PatchCollection(patches_nofill, alpha=0.75, match_original=True)
        
    #if len(patches) > 0:
    #    p0 = PatchCollection(patches, alpha=0.25, match_original=True)
    #    p1 = PatchCollection(patches, alpha=0.75, match_original=True)
                   
 
    # ax0: raw image
    ax0.imshow(input_image)
    if len(patches) > 0:
        ax0.add_collection(p0)
    if add_titles:
        ax0.set_title('Input Image + Ground Truth Buildings') 

    # truth polygons
    zero_arr = np.zeros(input_image.shape[:2])
    # set background to white?
    #zero_arr[zero_arr == 0.0] = np.nan
    ax1.imshow(zero_arr, cmap=cmap)
    if len(patches) > 0:
        ax1.add_collection(p2)
    if add_titles:
        ax1.set_title('Ground Truth Building Polygons')        

    # mask
    ax2.imshow(mask_image, cmap=cmap)
    # truth polygons?
    #if len(patches) > 0:
    #    ax1.add_collection(p1)
    if add_titles:
        ax2.set_title('Ground Truth Building Mask')    

    # distance transform
    cbar_pointer = ax3.imshow(dist_image)
    # overlay buildings on distance transform? 
    #if len(patches) > 0:
    #    ax3.add_collection(p1)
    if add_titles:
        #mind, maxd = np.round(np.min(dist_image),2), \
        #                                   np.round(np.max(dist_image),2)
        #dist_suffix = ""#" (min=" + str(mind) + ", max=" + str(maxd) + ")"
        #ax3.set_title("Yuan 2016 Distance Transform" + dist_suffix)
        ax3.set_title("Ground Truth Polygons Overlaid on Distance Transform")
    
    if colorbar:
        #from mpl_toolkits.axes_grid1 import make_axes_locatable
        #divider = make_axes_locatable(ax2)
        #cax = divider.append_axes('right', size='5%', pad=0.05)
        #fig.colorbar(cbar_pointer, cax=cax, orientation='vertical')
        left, bottom, width, height = [0.38, 0.85, 0.24, 0.03]
        cax = fig.add_axes([left, bottom, width, height])
        fig.colorbar(cbar_pointer, cax=cax, orientation='horizontal')

    #plt.axis('off')
    plt.tight_layout()
    if add_global_title:
        suptitle.set_y(0.95)
        fig.subplots_adjust(top=0.96)
    #plt.show()
 
    if len(plot_name) > 0:
        plt.savefig(plot_name)
    
    return    
    
    
#%%    
###############################################################################
def main():    

    imDir = os.path.join(spacenet_data_dir, '3band')
    vecDir = os.path.join(spacenet_data_dir, 'vectorData/geoJson')
    imDir_out = os.path.join(spacenet_explore_dir, '3band')

    ground_truth_patches = []
    pos_val, pos_val_vis = 1, 255
 
    ########################
    # Create directories

    #coordsDir = spacenet_explore_dir + 'pixel_coords_mask/'
    coords_demo_dir = os.path.join(spacenet_explore_dir, 'pixel_coords_demo')

    maskDir = os.path.join(spacenet_explore_dir, 'building_mask')
    maskDir_vis = os.path.join(spacenet_explore_dir, 'building_mask_vis')
    mask_demo_dir = os.path.join(spacenet_explore_dir, 'mask_demo')

    distDir = os.path.join(spacenet_explore_dir, 'distance_trans')
    dist_demo_dir = os.path.join(spacenet_explore_dir, 'distance_trans_demo')
    
    all_demo_dir = os.path.join(spacenet_explore_dir, 'all_demo')

    # make dirs
    for p in [imDir_out, coords_demo_dir, maskDir, maskDir_vis, mask_demo_dir,
              distDir, dist_demo_dir, all_demo_dir]:
        if not os.path.exists(p):
            os.mkdir(p)

    # get input images and copy to working directory
    rasterList = glob.glob(os.path.join(imDir, '*.tif'))[10:10+N_ims]   
    for im_tmp in rasterList:
        shutil.copy(im_tmp, imDir_out)
            
    # Create masks and demo images
    pixel_coords_list = []
    for i,rasterSrc in enumerate(rasterList):
        
        print (i, "Evaluating", rasterSrc)

        input_image = plt.imread(rasterSrc) # cv2.imread(rasterSrc, 1)
        
         # get name root
        name_root0 = rasterSrc.split('/')[-1].split('.')[0]
        # remove 3band or 8band prefix
        name_root = name_root0[6:]
        vectorSrc = os.path.join(vecDir, name_root + '_Geo.geojson')
        maskSrc = os.path.join(maskDir, name_root0 + '.tif')
        
        ####################################################
        # pixel coords and ground truth patches
        pixel_coords, latlon_coords = \
            geojson_to_pixel_arr(rasterSrc, vectorSrc, 
                                                      pixel_ints=True,
                                                      verbose=False)
        pixel_coords_list.append(pixel_coords)
       
        plot_name = os.path.join(coords_demo_dir, name_root + '.png')
        patch_collection, patch_coll_nofill = \
            plot_truth_coords(input_image, pixel_coords,   
                  figsize=(8,8), plot_name=plot_name,
                  add_title=False)
        ground_truth_patches.append(patch_collection)
        plt.close('all')
        ####################################################
        
        ####################################################
        #building mask
        outfile = os.path.join(maskDir, name_root0 + '.tif')
        outfile_vis = os.path.join(maskDir_vis, name_root0 + '.tif')
    
        # create mask from 0-1 and mask from 0-255 (for visual inspection)
        create_building_mask(rasterSrc, vectorSrc, 
                                                  npDistFileName=outfile, 
                                                  burn_values=pos_val)
        create_building_mask(rasterSrc, vectorSrc, 
                                                  npDistFileName=outfile_vis, 
                                                  burn_values=pos_val_vis)
        
        plot_name = os.path.join(mask_demo_dir, name_root + '.png')
        mask_image = plt.imread(outfile)    # cv2.imread(outfile, 0)
        plot_building_mask(input_image, pixel_coords,
                           mask_image,
                           figsize=(8,8), plot_name=plot_name,
                           add_title=False) 
        plt.close('all')
        ####################################################   
        
        ####################################################
        # signed distance transform
        # remove 3band or 8band prefix
        outfile = os.path.join(distDir, name_root0 + '.npy')#'.tif'    
        create_dist_map(rasterSrc, vectorSrc, 
                                        npDistFileName=outfile, 
                                        noDataValue=0, burn_values=pos_val, 
                                        dist_mult=1, vmax_dist=64)
        # plot
        #plot_name = os.path.join(dist_demo_dir + name_root, '_no_colorbar.png')
        plot_name = os.path.join(dist_demo_dir, name_root + '.png')
        mask_image = plt.imread(maskSrc)    # cv2.imread(maskSrc, 0)
        dist_image = np.load(outfile)
        plot_dist_transform(input_image, pixel_coords, 
                                                dist_image, figsize=(8,8),
                                                plot_name=plot_name, 
                                                add_title=False,
                                                colorbar=True)#False)
        plt.close('all')
        ####################################################

        ####################################################
        # plot all transforms
        plot_name = os.path.join(all_demo_dir, name_root + '.png')#+ '_titles.png'
        mask_image = plt.imread(maskSrc)    # cv2.imread(maskSrc, 0)
        dist_image = np.load(outfile)
        plot_all_transforms(input_image, pixel_coords, mask_image, dist_image, 
                        figsize=(8,8), plot_name=plot_name, add_global_title=False, 
                        colorbar=False, 
                        add_titles=False,#True,
                        poly_face_color='orange', poly_edge_color='red', 
                        poly_nofill_color='blue', cmap='bwr')        
        plt.close('all')
        ####################################################


###############################################################################    
if __name__ == '__main__':
    main()              
