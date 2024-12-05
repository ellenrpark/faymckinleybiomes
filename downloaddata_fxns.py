#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:09:47 2024

@author: epark
"""

import os
import xarray as xr
import gzip
import requests
import pandas as pd
import numpy as np
from scipy.interpolate import griddata

# Functions to download data for fay and mckinley biomes

def GetHadISSTData(datadir):
    
    # Data downloaded from:
    # https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html
    
    url_list = ['https://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST_ice.nc.gz',
                'https://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST_sst.nc.gz']
    
    # with gzip.open(requests.get(sst_url).content) as fp:
    #     ds = xr.open_dataset(fp)
        
    for url in url_list:
        # Dowload .gz file
        filename = url.split("/")[-1]
        
        print('Downloading and extracting: '+filename[:-3])
        filepath = datadir+'raw/'+filename
        with open(filepath, "wb") as f:
            r = requests.get(url)
            f.write(r.content)
        
        # Unzip .gz file
        with gzip.open(filepath) as fp:
            ds = xr.open_dataset(fp)
            ds.to_netcdf(filepath[:-3])
            ds.close()
            
        # Delete .gz file
        os.remove(filepath)
        
    return

def GetArgoMLD(datadir, updatestr='04142022'):
    
    # Data downloaded from:
    # https://mixedlayer.ucsd.edu/
    
    url = 'https://mixedlayer.ucsd.edu/data/Argo_mixedlayers_monthlyclim_'+updatestr+'.nc'
    
    filename = 'Argo_mixedlayers.nc'
    filepath = datadir+'raw/'+filename
    
    print('Downloading '+filename)
    data = xr.open_dataset(requests.get(url).content)
    data.to_netcdf(filepath)
    
    
    
    return

def GetMODIS(datadir, download_filename = 'AQUA_MODIS_L3m_CHL_list.csv'):
    
    # Make file to save chlorophyll 
    if os.path.exists(datadir+'raw/chl/') == False:
        os.mkdir(datadir+'raw/chl/')
        
    # Get 1ºx1º lat-lon grid from sst
    sst = xr.open_dataset(datadir+'raw/HadISST_sst.nc')
    sst.close()
    
    lon_gridvalues = sst.longitude.values
    lat_gridvalues = sst.latitude.values


    XX, YY = np.meshgrid(lon_gridvalues, lat_gridvalues)
    X = XX.flatten()
    Y = YY.flatten()
        
    # Need csv of necessary files
    flist = pd.read_csv(download_filename, header=None).loc[:,0].values
    
    for url in flist:
        
        filename = url.split('/')[-1]
        # Download monthly chlorophyll
        filepath = datadir+'raw/chl/' + filename
        
        print('Downloading '+filename)
        open(filepath, 'wb').write(requests.get(url).content)
        data = xr.open_dataset(filepath)
        data.close()
        
        # Re-gridd to 1º x 1º grid
        xx, yy = np.meshgrid(data.lon.values, data.lat.values)
        points = (xx.flatten(), yy.flatten())
        chl = griddata(points, data.chlor_a.values.flatten(), (X, Y), method='nearest').reshape(XX.shape)
        chl = chl.reshape(XX.shape)
        
        # Save data
        chl_ds = xr.Dataset(
                            data_vars=dict(
                                chlor_a=(["x", "y"], chl),
                            ),
                            coords=dict(
                                lat=(["x"], lat_gridvalues),
                                lon=(["y"], lon_gridvalues),
                            ),
                            attrs=dict(description="MODIS AQUA chlor_a downsampled to 1ºx1º from 9 km x 9km using griddata nearest neighbor"),
                            )
        outfname = filepath[:-6]+'1deg.nc'
        chl_ds.to_netcdf(outfname)
        
        # Delete initial downloaded data
        os.remove(filepath)
        
    
    return

def DownloadData(datadir):
    GetHadISSTData(datadir)
    GetArgoMLD(datadir)
    GetMODIS(datadir)
    return
    
