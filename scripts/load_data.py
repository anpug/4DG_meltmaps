#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 15:02:31 2023

@author: annpu
"""

import xarray as xr
import numpy as np
import pyproj 
from osgeo import gdal
import pickle

def preprocess_func(RCM):
    RCM.coords['lon_2'] = (RCM.coords['lon_2'] + 180) % 360 - 180
    return RCM

def tranform_coord(lon, lat): 
    transformer = pyproj.Transformer.from_crs('EPSG:4326','EPSG:3413')  
    polar_lon,polar_lat = transformer.transform(lat,lon)
    return polar_lon, polar_lat


def load_RCM(home_dir, RCM_name, year, season = None): 
    '''
    Function that loads in RCM from one netcdf-file 
    
    Input: RCM_name: Name of RCM. For HIRHAM I have merged daily files into one using cdo mergetime.
           year: year 
           season: Due to the large size of RACMO you can only load in 3 month at a time and you 
                   need to sepcify which 3 month to import. 
    
    Output: polar_lat, polar_lon; Coordinates in EPSG:3413. 
            melt_data: Data from netcdf file. Full year for MAR and HIRHAM and only 3 month for RACMO.
    '''

    if RCM_name == 'RACMO': 
        # Import Racmo, due to large dataset the annual files are devided into 4:
        if season == 'JFM':
            RCM = xr.open_dataset(home_dir + f'/RACMO2.3/Daily-1km/snowmelt.{year}_JFM.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc' )
        elif season == 'AMJ': 
            RCM = xr.open_dataset(home_dir + f'/RACMO2.3/Daily-1km/snowmelt.{year}_AMJ.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc' )
        elif season == 'JAS': 
            RCM = xr.open_dataset(home_dir + f'/RACMO2.3/Daily-1km/snowmelt.{year}_JAS.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc' )            
        elif season == 'OND':
            RCM = xr.open_dataset(home_dir + f'/RACMO2.3/Daily-1km/snowmelt.{year}_OND.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc' )
        
        # extract data as a array
        melt_data = RCM.snowmeltcorr.data
        
        # convert into m from km:
        polar_lat = RCM.LAT.data*1000
        polar_lon = RCM.LON.data*1000
        
    elif RCM_name == 'MAR': 
        # import nc-file:
        RCM= xr.open_dataset(home_dir + f'/marv3.12/MARv3.12.1-10km-daily-ERA5-{year}.nc')
        
        # Get meltdata as a array:
        melt_data = RCM.ME.data.squeeze()
        
        # Convert into m from km: 
        polar_x = RCM.x.data*1000
        polar_y = RCM.y.data*1000
        
        polar_lon,polar_lat = np.meshgrid(polar_x, polar_y)
        
    elif RCM_name == 'HIRHAM': 
        # Import nc-file: 
        RCM =  xr.open_dataset(home_dir + f'/HIRHAM/{year}/HIRHAM_{year}.nc')
        melt_data = RCM.snmel.data.squeeze()
        
        # Get latlon coordinates:
        lon = RCM.lon.data 
        lat = RCM.lat.data 
        
        # Transform into EPSG 3413 using pyproj: 
        polar_lon,polar_lat = tranform_coord(lon, lat)

    elif RCM_name == 'CARRA': 
        # Load in a full year of CARRA data using with lon from -180,180 degrees: 
        RCM = xr.open_mfdataset(home_dir + f'/CARRA/Daily2D_{year}*.nc', 
                                concat_dim = 'time' , combine = 'nested', 
                                preprocess=preprocess_func, parallel = True, data_vars = ['snmel'])

        # Extract lon og lat and reproject to North Polar Stereographic (EPSG:3431)
        lon = RCM.lon_2.data 
        lat = RCM.lat_2.data 
        polar_lon,polar_lat = tranform_coord(lon, lat)

        # Convert melt data to array: 
        melt_data = np.array(RCM.snmel.data)
        
    return polar_lat, polar_lon, melt_data

def load_yearly_RCM(RCM_name, year): 
    
    '''
    Function that loads a specific year of the RCM. Since RACMO are devided into seasonal files each file 
    needs to loaded sperately and combined. The rest of the RCMs are loaded from one file only. 
    
    Input: RCM_name: Name of RCM. For HIRHAM I have merged daily files into one using cdo mergetime.
           year: year 
    
    Output: polar_lat, polar_lon; Coordinates in EPSG:3413. 
            melt_data: Annual melt data from specified RCMS as an array. 
    
    '''   
    if RCM_name == 'RACMO': # Only for RACMO since files are too big to be loaded in at once.
        # Now we import a year and one RCM for the full year: (always on the same grid no matter the season)
        polar_lat, polar_lon, melt_dataJFM = load_RCM(RCM_name, year, season = 'JFM')
        _, _, melt_dataAMJ = load_RCM(RCM_name, year, season = 'AMJ')
        _, _, melt_dataJAS = load_RCM(RCM_name, year, season = 'JAS')
        _, _, melt_dataOND = load_RCM(RCM_name, year, season = 'OND')
        
        # Combine into one full year:
        melt_data = np.vstack((melt_dataJFM, melt_dataAMJ, melt_dataJAS, melt_dataOND))
    else: # All other RCMs:
        polar_lat, polar_lon, melt_data = load_RCM(RCM_name, year)
        
    return polar_lat, polar_lon, melt_data

def open_pickle(file_path): 
    f = open(file_path, 'rb')
    return pickle.load(f)

def read_tif(file_path):
    tif = gdal.Open(file_path)
    band= tif.GetRasterBand(1)
    arr = band.ReadAsArray().astype(np.float64)
    arr[arr == 65535] = np.nan
    return arr 

def write_geotiff(filename, arr, in_ds):
    if arr.dtype == np.float32:
        arr_type = gdal.GDT_Float32
    else:
        arr_type = gdal.GDT_Int32

    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(filename, arr.shape[1], arr.shape[0], 1, arr_type)
    out_ds.SetProjection(in_ds.GetProjection())
    out_ds.SetGeoTransform(in_ds.GetGeoTransform())
    band = out_ds.GetRasterBand(1)
    band.WriteArray(arr)
    band.FlushCache()
    band.ComputeStatistics(False)


def icemask_from_tif(file_path):
    '''
    Opens a meltmap (geotif-file) and convert it into a ice mask. Since the extent of the meltmaps are always 
    the same a random melt scene can be used to generate the icemask. Bedrock and ocean are set to np.nan whereas 
    ice is set to 1. 
    
    Input: file_path: Path to meltmap
    Output: icemask: a mask consising of either 1 or nan. 
    '''
    
    icemask = read_tif(file_path)
    # Ocean and bedrock is set to nan. 
    icemask[icemask >= 5] = np.nan
    # Different classes of melt is set to 1 
    icemask[icemask <= 4] = 1
    return icemask 


