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
import geopandas as gpd
from scipy.spatial import cKDTree
import pandas as pd
import rasterio as rio
import rioxarray as rxr



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
        
        # Variables not to be included in the data load:
        varlist = ['ahfl','ahfs','dlwrad','dswrad','evspsbl','rainfall','snfall','grndhflx','tsl2','gld','rogl',
                   'tradsm','sradsm','albedom','sn','supimp','tsl','rfrz']
        # Import nc-file: 
        #RCM =  xr.open_mfdataset(home_dir + f'/HIRHAM_2/Daily2D_{year}*.nc', 
        #                        parallel=True, drop_variables = varlist,
        #                        concat_dim='time', combine='nested')
        #RCM =  xr.open_mfdataset(home_dir + f'/HIRHAM/{year}/Daily2D_{year}*.nc', 
        #                        parallel=True, drop_variables = varlist,
        #                        concat_dim='time', combine='nested')
        # HIRHAM-ERA5_v2: 
        RCM =  xr.open_mfdataset(home_dir + f'/Daily2D_{year}*.nc', 
                                parallel=True, drop_variables = varlist,
                                concat_dim='time', combine='nested')
        
        melt_data = RCM.tas.data.squeeze()
        
        # Get latlon coordinates:
        lon = RCM.lon.data 
        lat = RCM.lat.data 
        
        # Transform into EPSG 3413 using pyproj: 
        polar_lon,polar_lat = tranform_coord(lon, lat)

        
    return polar_lat, polar_lon, melt_data

def load_yearly_RCM(home_dir, RCM_name, year): 
    
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
        polar_lat, polar_lon, melt_dataJFM = load_RCM(home_dir, RCM_name, year, season = 'JFM')
        _, _, melt_dataAMJ = load_RCM(home_dir, RCM_name, year, season = 'AMJ')
        _, _, melt_dataJAS = load_RCM(home_dir, RCM_name, year, season = 'JAS')
        _, _, melt_dataOND = load_RCM(home_dir, RCM_name, year, season = 'OND')
        
        # Combine into one full year:
        melt_data = np.vstack((melt_dataJFM, melt_dataAMJ, melt_dataJAS, melt_dataOND))
    else: # All other RCMs:
        polar_lat, polar_lon, melt_data = load_RCM(home_dir, RCM_name, year)
        
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
    band.SetNoDataValue(np.nan)  # Set np.nan as nodata
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


def load_ASCAT_grid(home_dir): 
    # Open nc file with the desired grid; (Made using gdal_translate)
    tif_grid = xr.open_dataset(home_dir + '/tif_grid.nc')

    # Generate ASCAT grid: 
    tif_x = tif_grid.x
    tif_y = tif_grid.y
    tif_lon, tif_lat = np.meshgrid(tif_x,tif_y)
    return tif_lat, tif_lon


def load_dmidata(home_dir, year, rcm_name): 
    varlist = ['ahfl','ahfs','dlwrad','dswrad','evspsbl','rainfall','snfall','grndhflx','tsl2','gld','rogl',
           'tradsm','sradsm','albedom','supimp','tsl','rfrz']
    if rcm_name == 'HIRHAM_ERA5': 
            # Open and concatenate netCDF files into one dataset
            #RCM = xr.open_mfdataset(home_dir+f'/HIRHAM_2/Daily2D_{year}*.nc',
            #                    parallel=True, drop_variables = varlist,
            #                    concat_dim='time', combine='nested')
            RCM =  xr.open_mfdataset(home_dir + f'/HIRHAM5-ERA5_v2/Daily2D_{year}*.nc', 
                                parallel=True, drop_variables = varlist,
                                concat_dim='time', combine='nested')
    
    elif rcm_name == 'HIRHAM_ERAI':
        if year >= 2015: 
            RCM =  xr.open_mfdataset(home_dir + f'/HIRHAM/{year}/HIRHAM_{year}.nc', drop_variables = varlist)
        
        else:
            RCM = xr.open_mfdataset(home_dir+f'/HIRHAM/{year}/Daily2D_{year}*.nc',
                                parallel=True, drop_variables = varlist,
                                concat_dim='time', combine='nested') 
    
    return RCM
        
def sample_xarray_points(RCM, aws_df): 
    RCM_xy = RCM.sel(x = np.array(aws_df.geometry.x), y = np.array(aws_df.geometry.y), method = 'nearest')

    # Nu skal det så ind i en dataframe 
    data4df = np.zeros((np.size(RCM_xy.time), len(aws_df)))
    for i in range(len(aws_df)):
        # Set up data for dataframe:
        data4df[:, i] = np.array(RCM_xy.isel(x=i, y=i).ME).squeeze()

    #Put into dataframe:
    sampled_df = pd.DataFrame(data = data4df, 
                index = RCM_xy.time.dt.floor('d'),
                columns = list( aws_df.index))

    return sampled_df

# Function that loads in RCM from one netcdf-file in a given year and season. The functions then samples the data in the AWS locations and only returns the melt data at these locations.
def get_melt_in_points(home_dir, date_range, aws_df, rcm_name): 
    # Create a dataframe: 
    AWS_RCM = pd.DataFrame(columns = list(aws_df.index))

    # Since MAR and RACMO data structure are in the data structure the method is very similar:
    if rcm_name == 'MAR' or rcm_name == 'RACMO' or rcm_name == 'RACMO_s': 

        for year in range(date_range[0], date_range[1]+1):
            # Load RCM:
            
            if rcm_name == 'MAR':
                RCM = xr.open_dataset(home_dir + f'/marv3.12/MARv3.12.1-10km-daily-ERA5-{year}.nc')
                # For MAR we need to correct the xy to m:
                RCM.coords['x'] = RCM.coords['x']*1000
                RCM.coords['y'] = RCM.coords['y']*1000
                # Rename TIME to time:
                RCM = RCM.rename_vars({'TIME': 'time'})
                sampled_df = sample_xarray_points(RCM, aws_df)

            elif rcm_name == 'RACMO':
                seasons = ['AMJ', 'JAS', 'OND']
                RCM = xr.open_dataset(home_dir + f'/RACMO2.3/Daily-1km/snowmelt.{year}_JFM.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc')
                # Rename snowmeltcorr to ME: 
                RCM = RCM.rename_vars({'snowmeltcorr': 'ME'})
                sampled_df = sample_xarray_points(RCM, aws_df)
                for s in seasons: # load in the rest of the seasons:
                    RCM = xr.open_dataset(home_dir + f'/RACMO2.3/Daily-1km/snowmelt.{year}_{s}.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc')
                    # Rename snowmeltcorr to ME: 
                    RCM = RCM.rename_vars({'snowmeltcorr': 'ME'})
                    sampled_df_s = sample_xarray_points(RCM, aws_df)
                    sampled_df = pd.concat([sampled_df, sampled_df_s])

            elif rcm_name == 'RACMO_s': # only the JAS t2m data
                
                RCM = xr.open_dataset(home_dir + f'/RACMO2.3/t2m/t2m.{year}_JAS.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc')
                # remane t2m to TT: 
                RCM = RCM.rename_vars({'t2m': 'TT'})
                # Convert from km to m: 
                #RCM.coords['x'] = RCM.coords['x']*1000
                #RCM.coords['y'] = RCM.coords['y']*1000  
                # Sample in points:
                sampled_df = sample_xarray_points(RCM, aws_df)
                

            # concat onto dataframe:
            AWS_RCM = pd.concat([AWS_RCM, sampled_df])
            print(f'Dataframe updated for {year}\n')

    elif rcm_name == 'HIRHAM_ERA5' or rcm_name == 'HIRHAM_ERAI': 
       
        # Load HIRHAM for the first year:
        RCM = load_dmidata(home_dir, date_range[0], rcm_name)
        
        print('Loaded data for:', date_range[0])
        # Convert to the correct time format: (datetime64)
        time_index = np.arange(np.datetime64(f'{date_range[0]}-01-01'), 
                               np.datetime64(f'{date_range[0]+1}-01-01'), 
                               np.timedelta64(1, 'D'))

        # Flatten and stack latitude and longitude values to create a 2D grid
        lon = np.array(RCM.lon.data).flatten()
        lat = np.array(RCM.lat.data).flatten()
        lat_lon_grid = np.column_stack((lat, lon))

        # Build a KD-tree for efficient nearest neighbor searches
        tree = cKDTree(lat_lon_grid)
        
        # Extract AWS coordinates and find their nearest neighbors in the grid
        aws_coords = np.column_stack([aws_df.lat, aws_df.lon])
        dist , closest_indices = tree.query(aws_coords, k=1)

        print('Nearest coordinates found\n')
        # Reshape melt data for compatibility with the lat-lon grid
        melt_data = np.array(RCM.snmel.data).reshape(RCM.snmel.data.shape[0], -1)
        
        # Extract melt data for the nearest grid points
        melt_data_aws = melt_data[:, closest_indices]

        # Create a DataFrame with the melt data
        sampled_df = pd.DataFrame(data=melt_data_aws, 
                                  index=time_index,
                                  columns=list(aws_df.index))
        
        AWS_RCM = pd.concat([AWS_RCM, sampled_df])

        # Now sample the data in the remaining years: 
        for year in range(date_range[0]+1, date_range[1]+1):
            print(year)
            # Open files and combine into an annual dataset:
            RCM = load_dmidata(home_dir, year, rcm_name)
            print('data loaded')
            # Reshape melt data for compatibility with the lat-lon grid
            melt_data = np.array(RCM.snmel.data).reshape(RCM.snmel.data.shape[0], -1)
            
            # Extract melt data for the nearest grid points
            melt_data_aws = melt_data[:, closest_indices]

            time_index = np.arange(np.datetime64(f'{year}-01-01'), 
                                           np.datetime64(f'{year+1}-01-01'), 
                                           np.timedelta64(1, 'D'))
            
            # Create a DataFrame with the melt data
            sampled_df = pd.DataFrame(data=melt_data_aws, 
                                  index=time_index,
                                  columns=list(aws_df.index))
            AWS_RCM = pd.concat([AWS_RCM, sampled_df])

    return AWS_RCM


def load_AWS_data(id): 
    old_names = dict({'CP1': 'CrawfordPoint1', 'DY2': 'DYE-2', 'HUM': 'Humboldt', 'NAE': 'NASA-E', 'NAU': 'NASA-U', 
                  'NEM': 'NEEM', 'NSE': 'NASA-SE','SDM': 'SouthDome', 'SDL': 'Saddle', 'TUN': 'Tunu-N', 
                  'SWC': 'SwissCamp','Summit': 'Summit', 'PetermannELA': 'PetermannELA', 'GITS': 'GITS'})

    # If the AWS station is part of the old gc-net when we want to open the historical data:
    if id in ['CP1', 'DY2', 'HUM', 'NAE', 'NAU', 'NEM', 'NSA', 'NSE','SDM', 'SDL', 'TUN', 'SWC', 'Summit', 'PetermannELA', 'GITS']:
        id = old_names[id]
        # Load in AWS data:
        AWS_data = pd.read_csv(home_dir + f'/AWS/combined/historical/{id}_daily.csv', header = 14,low_memory=False)
        AWS_data= AWS_data.iloc[10:] # Remove first 10 rows where there no data. 
        AWS_data = AWS_data.astype({'TA1': 'float', 'TA2': 'float', 'TA3': 'float', 'TA4': 'float'}) # Convert to float
        
        # Compute mean of TA1, TA2, TA3, TA4 and store in a new column named t_u
        AWS_data['t_u'] = AWS_data[['TA1', 'TA2', 'TA3', 'TA4']].mean(axis=1)

        AWS_data = AWS_data.rename(columns = {AWS_data.columns[0]: 'time'}) # rename the first column to time
        #AWS_data.time = AWS_data.time.tz_localize(None) # Remove timezone info


    else: 
        # Load in AWS that part of PROMICE ogrinally:
        AWS_data = pd.read_csv(home_dir + f'/AWS/combined/day/{id}_day.csv') # Load in AWS data
    
    # Convert to datetime and set time as index:
    AWS_data['time'] = pd.to_datetime(AWS_data['time']) # Convert to datetime
    AWS_data = AWS_data.set_index('time') # Set time as index
    AWS_data.index = AWS_data.index.tz_localize(None) # remove tz 
    # select only t_u: 
    return AWS_data[['t_u']]


def load_aws_locations(home_dir): 
    # Load in AWS stations and mask applied to select AWS stations: 
    aws_loc = pd.read_csv(home_dir + '/AWS/combined/AWS_latest_locations.csv')
    ASCAT_mask =gpd.read_file(home_dir + '/QGIS/ASCAT_mask.geojson') 
    # set station ID as index and remove stations where we don't have datafiles from:
    aws_loc = aws_loc.set_index('stid')
    aws_loc = aws_loc.drop(['QAS_Lv3','ZAK_L','WEG_L', 'KPC_Uv3', 'KPC_Lv3']) # remove nan entry and st where I dont have data. 

    # Convert to geopandas dataframe and reproject to EPSG:3413:
    gdf = gpd.GeoDataFrame(aws_loc, geometry=gpd.points_from_xy(aws_loc.lon, aws_loc.lat), crs="EPSG:4326").to_crs('EPSG:3413')

    # Clip to ASCAT mask
    clipped_aws_loc = gpd.clip(gdf, ASCAT_mask) 
    return clipped_aws_loc

def make_data_mask(home_dir):
    # Load snowline and ascat scene to compte ascat mask:
    snow_line = gpd.read_file(home_dir + '/snowline/2020_snowline_epsg3413.geojson')
    ascat_scene = rxr.open_rasterio(home_dir + '/ASCAT/2020/meltmap_v02_2020-01-07.tif')

    # Clip the ascat scene to the snowline:
    masked_data = ascat_scene.rio.clip(snow_line.geometry, snow_line.crs, drop=False)

    # convert the masked data to a numpy array dtype float:
    mask= masked_data.to_numpy().astype('float64').squeeze()
    # Where masked array is 0, set to nan:
    mask[mask == 0] = 1
    # where masked array is over zero, set to 1:
    mask[mask > 1] = np.nan

    # convert the ascat scene to a numpy array dtype float to get the unmasked data:
    unmasked = ascat_scene.data.astype('float64').squeeze()
    unmasked[unmasked == 0] = 1
    # where masked array is over zero, set to 1:
    unmasked[unmasked > 1] = np.nan

    return mask, unmasked