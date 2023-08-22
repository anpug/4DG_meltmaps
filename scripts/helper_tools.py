#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 11:38:49 2023

@author: annpu
"""

import numpy as np 



def number_of_meltdays(rcm, threshold): 
    '''
    Compute the number of meltdays in a year for a RCM. Meltdays are thresholded by a magnitude threshold only. 
    '''    
    melt = np.zeros(rcm.shape)
    #melt = np.full([365, 540, 299], np.nan)
    melt[rcm >= threshold] = 1
    sum_meltdays = np.sum(melt, axis = 0)
    sum_meltdays[sum_meltdays == 0] = np.nan
    return sum_meltdays 


        
def rolling_window(a, size):
    '''
    Takes an array (1D or 2D) and compute an array using a rolling window (2D or 3D array). The rollig array is 
    created along the first axis (i.e. time). 
    
    Input: a: The original array that the rolling array needs to be computed from. 
           size: size of the rolling window. 
    
    Output: The rolling window array. 
    '''
    
    shape = a.shape[:-1] + (a.shape[-1] - size + 1, size)
    strides = a.strides + (a.strides[-1],)
    
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def find_first_last(arr):
    '''
    Takes a 3D boolean array and finds the first and last occurence of True along the zero axis (time in this case).
    
    Input: arr: 3D boolean array
    Output: first_occurrence: The index of the first True occurence. 
            last_occurence: The index of the last True occurance
    '''
    
    # Find the index of the first occurence:  
    first_occurrence = np.argmax(arr, axis=0)
    # Flip array around and find the last occurence. Since the index correspond to the flipped array, substract 
    # with the size of array: 
    last_occurrence = arr.shape[0] - np.argmax(np.flip(arr, axis=0), axis=0) - 1
    return first_occurrence.astype(np.float32),  last_occurrence.astype(np.float32)


def melt_season(array, thresholds, length):
    
    '''
    Function used to compute the start and end of the melt season. The melt season is defined using a threshold 
    for the magnitude (either mmWE/day for RCM or class threshold for ASCAT) and a length of consective melt 
    threshold. 
    
    Input: array: a year of either daily RCM active melt or ASCAT melt classification (3D array)
           threshold: lower and upper magnitude threshold (the upper threshold for RCM should be very large as
                      there should be no upper threshold)
           length: Length threshold of how many days of consective melt defines the start and end of melt season. 
    Output; Start and end of the melt season for each pixel. 
    '''
    
    outarray = np.zeros_like(array).astype('bool') # Zero array 
    subarray1 = np.ones(length) * thresholds[0] # Lower magnitude threhold
    subarray2 = np.ones(length) * thresholds[1] # Upper magnitude threhold
    # Compute boolean for length and magnitude thresholds is statisfied: 
    index1 = np.all(rolling_window(array, length) >= subarray1, axis=2) # Lower magnitude threholding 
    index2 = np.all(rolling_window(array, length) <= subarray2, axis=2) # Upper magnitude threholding 
    index = index1 & index2 # Both magnitude thresholds 
    # Correct indexing
    outarray[:, :index.shape[1]] = index
    # Find start and end: 
    start, end = find_first_last(outarray.T)
    start[start == 0] = np.nan # Set no start of melt season to nan 
    start = start.reshape((540,299)) # Correct grid shape 
    end[end >= 364] = np.nan # Set no end of melt season to nan 
    end = end.reshape((540,299)) # Correct grid shape 
    end = end + length - 1 # correct indexing since it has to be the last day where the thresholding was satisfied
    return start, end


def rolling_window_rcm(array, threshold, length):
    outarray = np.zeros_like(array).astype('bool')
    subarray = np.ones(length) * threshold
    index = np.any(rolling_window(array, length) >= subarray, axis=2)
    outarray[:, :index.shape[1]] = index
    return outarray.T

def window_comparison(ascat, rcm, threshold, length): 
    
    '''
    Compute find the ratio (procentage) between days with melt/refreezing/no melt in both ASCAT and RCM 
    and total number of days with melt/refreezing/no melt in ASCAT. It is possible to specify a window in which 
    the RCM "looks into the future" using the length variable. If length = 1, no window is applied. 
    
    Input: ascat: Boolean-like that specifies if selected ascat claasification is satisfied. 
           rcm: RCM with daily active melt. 
           threshold: Threshold for when a day is clssified as a melt day in the RCM. 
           length: The number of days in the future included in the RCM used to account for the preporsseing 
                   window of ASCAT. If length = 1, no window is applied. 
    
    Output: The melt/refreezing/no melt ratio for each pixel in the ascat grid. Pixels where melt/refreezing/no melt
            is not detecetd by ASCAT is set to nan. 
    '''
    
    # Computed boolean array melt occurs over a threshold inside a window of specified length: 
    window_rcm = rolling_window_rcm(rcm.reshape((rcm.shape[0], -1)).T,threshold, length)
    melt_rcm = window_rcm.reshape(rcm.shape) # reshape into original dimentions to match ASCAT 
    
    # Days where 
    match = np.equal(ascat.astype('bool'), melt_rcm).astype('float')
    # Exclude indexs where ASCAT is False: 
    match[ascat == 0] = np.nan
    # Sum days: 
    sum_match = np.nansum(match, axis = 0)
    sum_ascat = np.nansum(ascat, axis = 0)
    return np.divide(sum_match,sum_ascat)*100 # Computed ratio. 