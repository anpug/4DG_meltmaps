#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 11:38:49 2023

@author: annpu
"""

from scipy.interpolate import griddata
import numpy as np
import cartopy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
import matplotlib.colors as colors
import seaborn as sns



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

def number_of_meltdays(rcm, threshold): 
    '''
    Compute the number of meltdays in a year for a RCM. Meltdays are thresholded by a magnitude threshold only. 
    '''    
    melt = np.zeros(rcm.shape)
    melt[rcm >= threshold] = 1
    sum_meltdays = np.sum(melt, axis = 0)
    sum_meltdays[sum_meltdays == 0] = np.nan
    return sum_meltdays 

def compute_melt_extent(rcm, threshold): 
    '''
    Compute the number of meltdays in a year for a RCM. Meltdays are thresholded by a magnitude threshold only. 
    '''    
    melt = np.zeros(rcm.shape)
    #melt = np.full([365, 540, 299], np.nan)
    melt[rcm >= threshold] = 1
    melt_extent = np.sum(melt.reshape(melt.shape[0], -1), axis = 1)
    return melt_extent/59959*100


def compute_max_melt_extent(melt_days, ascat_mask):
    melt_days[(melt_days.astype(str)==str(np.nan)) & (ascat_mask.astype(str) != str(np.nan)) ] = 0
    # Compute extent:
    extent = np.zeros(melt_days.shape)
    extent[melt_days*ascat_mask>0] = 1
    return np.sum(extent)/66427*100


# For consistent plotting on Greenland
def plot_melt_volume(data, ax, transform = cartopy.crs.epsg('3413'), mask = None, 
                   vmin = 0, vmax = 100, 
                   cbar = False, cbar_label = None, cmap = 'magma'): 
    grid = xr.open_dataset(home_dir + '/tif_grid.nc')
    grid_X, grid_Y = np.meshgrid(grid.x,grid.y)
    
    if mask is not None: 
        # Mask data: _r 
        data = data * mask 
    
    plot = ax.scatter(grid_X, grid_Y, c = np.flip(data+0.00001, axis =0 ), 
               transform = transform, s = 0.35, cmap = cmap, 
               #vmin = vmin, vmax = vmax, 
               norm = colors.LogNorm(vmin = 0.01, vmax = 30000 )) # ,vmin= 0,vmax = 2000
    #ax.contour(grid_X, grid_Y, np.flip(data+0.00001, axis =0 ), transform = transform, 
    #           levels = 100, norm = colors.LogNorm( ), colors = 'black')
    #ax.add_feature(cartopy.feature.LAND,color = 'tan', alpha = 0.7)
    ax.add_feature(cartopy.feature.LAND,color = '#C6C7C8', alpha = 0.7)
    
    ax.coastlines(resolution='50m', linewidth = 0.05, color = 'black')
    ax.gridlines(linewidth=0.8, linestyle = ':')

    if cbar:
        plt.colorbar(plot , ax = ax , label = cbar_label)
    return plot

# For consistent plotting on Greenland
def plot_greenland(data, ax, transform = cartopy.crs.epsg('3413'), mask = None, 
                   vmin = 0, vmax = 100, 
                   cbar = False, cbar_label = None, cmap = sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True), text_label = None): 
    if mask is not None: 
        # Mask data: 
        data = data * mask 
    grid = xr.open_dataset(home_dir + '/tif_grid.nc')
    grid_X, grid_Y = np.meshgrid(grid.x,grid.y)
    plot = ax.scatter(grid_X, grid_Y, c = np.flip(data, axis =0 ), 
               transform = transform, s = 0.5, cmap = cmap, 
               vmin = vmin, vmax = vmax)
    ax.add_feature(cartopy.feature.LAND,color = '#C6C7C8', alpha = 0.7)
    ax.coastlines(resolution='10m', linewidth = 0.05, color = 'black')
    ax.gridlines(linewidth=0.3)

    if cbar:
        plt.colorbar(plot , ax = ax , label = cbar_label)


    if text_label is not None: 
        ax.text(0.99, 0.01, text_label, transform=ax.transAxes, fontsize=16,  va='bottom', ha='right')

    return plot