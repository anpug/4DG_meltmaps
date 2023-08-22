#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 15:02:31 2023

@author: annpu
"""

import xarray as xr
import numpy as np
import pyproj 
from scipy.interpolate import griddata


def transform_coords(lon, lat): 
    transformer = pyproj.Transformer.from_crs('EPSG:4326','EPSG:3413')  
    polar_lon,polar_lat = transformer.transform(lat,lon)
    return polar_lon, polar_lat

def _preprocess(RCM):
    RCM.coords['lon_2'] = (RCM.coords['lon_2'] + 180) % 360 - 180
    return RCM
    
    return RCM 

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
        
        # Trnasform coordinates to polar stereografic: 
        polar_lon,polar_lat = transform_coords(lon, lat)
        
    elif RCM_name == 'CARRA': 
        RCM = xr.open_mfdataset(home_dir + f'/CARRA/Daily2D_{year}*.nc', 
                                concat_dim = 'time' , 
                                combine = 'nested', 
                                preprocess=_preprocess, 
                                parallel = True, 
                                data_vars = ['snmel'])
        
        # For some reason they are called lon_2 and lat_2:
        lon = RCM.lon_2.data 
        lat = RCM.lat_2.data 

        # Trnasform coordinates to polar stereografic: 
        polar_lon,polar_lat = transform_coords(lon, lat)
        melt_data = np.array(RCM.snmel)
        
    return polar_lat, polar_lon, melt_data

def load_yearly_RCM(home_dir,RCM_name, year): 
    
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
        polar_lat, polar_lon, melt_dataJFM = load_RCM(home_dir,RCM_name, year, season = 'JFM')
        _, _, melt_dataAMJ = load_RCM(home_dir,RCM_name, year, season = 'AMJ')
        _, _, melt_dataJAS = load_RCM(home_dir,RCM_name, year, season = 'JAS')
        _, _, melt_dataOND = load_RCM(home_dir,RCM_name, year, season = 'OND')
        
        # Combine into one full year:
        melt_data = np.vstack((melt_dataJFM, melt_dataAMJ, melt_dataJAS, melt_dataOND))
    else: # All other RCMs:
        polar_lat, polar_lon, melt_data = load_RCM(home_dir,RCM_name, year)
        
    return polar_lat, polar_lon, melt_data


def grid_data(lat_ori, lon_ori, var_ori, ref_lat, ref_lon, method_name):
    
    '''
    Function used to interpolate/regrid melt data to the ASCAT grid (in this case). Important that the original 
    grid and new grid is in the same projection (in this case EPSG:3413)
    
    Input: lat_ori, lon_ori: coordinates of the original grid (2D arrays)
           var_ori: Data to be regridded on the original grid (2D arrays).
           ref_lon, ref_lat: 
           method: 
    
    Output: (melt)Data on the new grid. 
    '''
   
    grid_var = griddata((lon_ori.flatten(), lat_ori.flatten()), var_ori.flatten(), (ref_lon, ref_lat), method=method_name)

    return grid_var
