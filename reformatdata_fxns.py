#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 12:43:15 2024

@author: epark
"""

import numpy as np
import xarray as xr
import glob 
import pandas as pd



def ReformatIce(datadir, end_year):
    # Load data
    datapath = datadir+'raw/'
    yr_range = np.arange(2010, end_year+1).astype('int')

    data = xr.open_dataset(datapath+'HadISST_ice.nc')
    data.close()
    
    yrs =  np.array([pd.Timestamp(x.astype('datetime64[s]')).year for x in data.time.values])
    mns =  np.array([pd.Timestamp(x.astype('datetime64[s]')).month for x in data.time.values])
    
    # data_time = xr.open_dataset(InDir+'HadISST_ice.nc')
    # yrs =  np.array([pd.Timestamp(x.astype('datetime64[s]')).year for x in data_time.time.values])
    # mns =  np.array([pd.Timestamp(x.astype('datetime64[s]')).month for x in data_time.time.values])
    # data_time.close()
    
    ice = np.zeros((yr_range.shape[0],data.latitude.values.shape[0], data.longitude.values.shape[0]))* np.NaN
    
    for i in np.arange(yr_range.shape[0]):
        
        yr = yr_range[i]
        
        inds= np.where(yrs==yr)[0]
        
        sub_ice = data.sic.values[inds,:,:]
        
        # if yr != 2023:
        #     sub_ice = data.sic.values[inds,:,:]
        # else:
        #     # Missing December so replace with climatology
        #     # of december
        #     m_inds = np.where((mns == 12) & (yrs >=2010))[0]
        #     december = np.nanmean(data.sic.values[m_inds,:,:], axis = 0)
            
        #     sub_ice = np.concatenate((data.sic.values[inds,:,:],
        #                               december.reshape(1,december.shape[0],december.shape[1])))
        
        
        # Calculate maximum sea ice fraction from the year
        ice[i,:,:]=np.nanmax(sub_ice, axis=0)
        
    ice_ds = xr.Dataset(
                        data_vars=dict(
                            ice=(["x", "y","z"], ice),
                        ),
                        coords=dict(
                            yr=(["x"], yr_range),
                            lat=(["y"], data.latitude.values),
                            lon=(["z"], data.longitude.values),
                        ),
                        attrs=dict(description="Annual minimum ice fraction from HadISST_ice"),
                        )
    return ice_ds
   
    
    
def ReformatSST(datadir, end_year):
    
    datapath = datadir+'raw/'
    yr_range = np.arange(2010, end_year+1).astype('int')
    
    # Load data
    data = xr.open_dataset(datapath+'HadISST_sst.nc')
    data.close()
    
    yrs =  np.array([pd.Timestamp(x.astype('datetime64[s]')).year for x in data.time.values])
    mns =  np.array([pd.Timestamp(x.astype('datetime64[s]')).month for x in data.time.values])
    
    
    sst = np.zeros((yr_range.shape[0],data.latitude.values.shape[0], data.longitude.values.shape[0]))* np.NaN
    
    for i in np.arange(yr_range.shape[0]):
        
        yr = yr_range[i]
        
        inds= np.where(yrs==yr)[0]
        
        sub_sst = data.sst.values[inds,:,:]
        
        # if yr != 2023:
        #     sub_sst = data.sst.values[inds,:,:]
        # else:
        #     # Missing December so replace with climatology
        #     # of december
        #     m_inds = np.where((mns == 12) & (yrs >=2010))[0]
        #     december = np.nanmean(data.sst.values[m_inds,:,:], axis = 0)
            
        #     sub_sst = np.concatenate((data.sst.values[inds,:,:],
        #                               december.reshape(1,december.shape[0],december.shape[1])))
        
        # Calculate mean annual sst
        sst[i,:,:]=np.nanmean(sub_sst, axis=0)
        
    sst_ds = xr.Dataset(
                        data_vars=dict(
                            sst=(["x", "y","z"], sst),
                        ),
                        coords=dict(
                            yr=(["x"], yr_range),
                            lat=(["y"], data.latitude.values),
                            lon=(["z"], data.longitude.values),
                        ),
                        attrs=dict(description="Annual average sst from HadISST_sst"),
                        )
    return sst_ds
    
def ReformatMLD(datadir, end_year):
    
    datapath = datadir+'raw/'
    yr_range = np.arange(2010, end_year+1).astype('int')
    
    # Load data
    data = xr.open_dataset(datapath+'Argo_mixedlayers.nc')
    
    # Get climatological maximum mixed layer depth
    mld = np.nanmax(data.mld_da_mean.values, axis = 2)
    
    mld_ds = xr.Dataset(
                    data_vars=dict(
                        mld=(["x", "y"], mld),
                    ),
                    coords=dict(
                        lat=(["x"], data.lat.values),
                        lon=(["y"], data.lon.values),
                    ),
                    attrs=dict(description="Maximum MLD from density algorithm from argo mld climatology"),
                    )
    
    return mld_ds
    
def ReformatCHL(datadir, end_year):
    
    datapath = datadir+'raw/chl/'
    yr_range = np.arange(2010, end_year+1).astype('int')
    
    # Mean chl
    # N. Hemisphere: mean of April-September
    # S. Hemisphere: mean of December -March
    # Lat division btwn N/S: 10S
    
    chl = np.zeros((yr_range.shape[0],180,360))*np.NaN
    
    for i in np.arange(yr_range.shape[0]):
        yr = yr_range[i]
        
        # Load Northern Hemisphere data
        nh_mns = [str(yr)+'04', str(yr)+'05', str(yr)+'06', str(yr)+'07', str(yr)+'08', str(yr)+'09']
        nh_chl = np.zeros((6, 100, 360))*np.NaN
        
        for j in np.arange(len(nh_mns)):
            m = nh_mns[j]
            # Load data
            
            fname = glob.glob(datapath+'*'+m+'*.nc')
            
            # if len(fname) == 0:
                
            #     # 2023 year...data not available
            #     # mean of 2010-2022
            #     month_vals = np.zeros((13, 100, 360))*np.NaN
                
            #     month_years = np.arange(2010, 2023)
            #     for mi in np.arange(month_years.shape[0]):
            #         m_yr = month_years[mi]
                    
            #         m_list = [str(m_yr)+m[-2:]]
                    
            #         fname_m = glob.glob(Dir+'*'+m_list[0]+'*.nc')[0]
                    
            #         data = xr.open_dataset(fname_m)
            #         # Crop to desired range
            #         lat_inds = np.where(data.lat.values>-10)[0]
            #         subdata = data.chlor_a.values[lat_inds,:]
                    
            #         # Concat array
            #         month_vals[mi,:,:]=subdata
                    
            #         data.close()
                        
                
            #     # Get the monthly mean
            #     nh_chl[j,:,:]=np.nanmean(month_vals, axis = 0)
                
            # else:
                
            fname = fname[0]
            data = xr.open_dataset(fname)

            # Crop to desired range
            lat_inds = np.where(data.lat.values>-10)[0]
            subdata = data.chlor_a.values[lat_inds,:]
            
            # Concat array
            nh_chl[j,:,:]=subdata
            
            data.close()
            
        # Take mean of array
        chl[i,0:100, :]=np.nanmean(nh_chl, axis = 0)
        
        # Load Southern Hemisphere data
        sh_mns = [str(yr-1)+'12', str(yr)+'01', str(yr)+'02', str(yr)+'03']
        sh_chl = np.zeros((6, 80, 360))*np.NaN
        for j in np.arange(len(sh_mns)):
            m = sh_mns[j]
            # Load data
            fname = glob.glob(datapath+'*'+m+'*.nc')[0]
            
            data = xr.open_dataset(fname)

            # Crop to desired range
            lat_inds = np.where(data.lat.values<=-10)[0]
            subdata = data.chlor_a.values[lat_inds,:]
            
            # Concat array
            sh_chl[j,:,:]=subdata
        
            data.close()
            
        chl[i,100:, :]=np.nanmean(sh_chl, axis = 0)
        

    chl_ds = xr.Dataset(
                        data_vars=dict(
                            chl=(["x", "y","z"], chl),
                        ),
                        coords=dict(
                            yr=(["x"], yr_range),
                            lat=(["y"], data.lat.values),
                            lon=(["z"], data.lon.values),
                        ),
                        attrs=dict(description="Northern/Southern Hemisphere summer/spring mean [chlora] from AQUA MODIS data following Fay-McKinley"),
                        )
    return chl_ds

def ReformatData(datadir, end_year = 2023):
    
    ice = ReformatIce(datadir, end_year)
    sst = ReformatSST(datadir, end_year)
    mld = ReformatMLD(datadir, end_year)
    chl = ReformatCHL(datadir, end_year)
    
    return ice, sst, mld, chl
    
    
    
