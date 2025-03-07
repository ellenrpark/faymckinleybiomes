#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:05:20 2024

@author: epark
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import  NearestNDInterpolator
from global_land_mask import globe
import scipy.ndimage
import matplotlib.colors as mcolors
import numpy as np
import os
from matplotlib import gridspec
import matplotlib.ticker as mticker

def GetBiomeColorMap():
    north_pac_color = plt.cm.Blues_r(np.linspace(0, 0.6, 4))
    south_pac_color = plt.cm.Purples(np.linspace(0.3, 1, 3))

    north_atl_color = plt.cm.Reds_r(np.linspace(0, 0.8, 4))
    south_atl_color = plt.cm.Oranges(np.linspace(.5, 1, 2))

    indian_color = plt.cm.YlOrRd(np.linspace(0.3, 0.4, 1))

    south_oc_color = plt.cm.Greens(np.linspace(0.2, 1, 3))

    biome_colors = np.vstack((north_pac_color,
                              south_pac_color,
                              north_atl_color,
                              south_atl_color,
                              indian_color, 
                              south_oc_color))

    cmm = mcolors.LinearSegmentedColormap.from_list('my_colormap',biome_colors,17)
    
    return cmm

def CalculateBiomes(savedir, end_year, ice_data, sst_data, mld_data, chl_data):

    lat_range = ice_data.lat.values
    lon_range = ice_data.lon.values
    
    lon, lat = np.meshgrid(lon_range, lat_range)
    
    # Determine which coordinates are land vs. ocean
    landmask = globe.is_land(lat, lon)
    
    # Initialize biome flag array
    yr_range = np.arange(2010, end_year+1).astype('int')
    biomes_all = np.zeros((yr_range.shape[0]+1, 180, 360))*np.NaN
    change_all = np.zeros((yr_range.shape[0]+1, 180, 360))*np.NaN
    
    # Load original Fay and McKinley biomes to define marginal seas
    # and coastal areas
    fm_biomes = xr.open_dataset('Fay_McKinley_Time_Varying_Biomes.nc')
    
    # Load biome names
    bnames = pd.read_csv('biomes_Fay.csv')
    
    # Get colormap 
    cmm = GetBiomeColorMap()
    
    # Divide oceans into Atlantic, Pacific, Indian, and Southern Oceans
    # NH Pacific lon <= -101 or lon >= 109
    # NH Atlanic -101 < lon < 109
    
    # Make directory to save figures
    if os.path.exists(savedir+'figures/') == False:
        os.mkdir(savedir+'figures/')
    figdir = savedir+'figures/'
    
    # For each year, calculate biome flags
    for i in np.arange(biomes_all.shape[0]):
        
        # For each year, determine biomes
        if i == biomes_all.shape[0]-1:
            # Period average
            print('Calculating mean biomes')
            sea_ice_frac = np.nanmean(ice_data.ice.values,axis = 0)
            sst = np.nanmean(sst_data.sst.values, axis = 0)
            chl = np.nanmean(chl_data.chl.values, axis = 0)
            mld = np.flipud(mld_data.mld.values[:,:])
        else:       
            print('Calculating biomes for ',yr_range[i])
            sea_ice_frac = ice_data.ice.values[i,:,:]
            sst = sst_data.sst.values[i,:,:]
            chl = chl_data.chl.values[i,:,:]
            mld = np.flipud(mld_data.mld.values[:,:])
            
            
       #########################
       ## NORTHERN HEMISPHERE
       #########################
        
        ## ICE BIOMES
        # Split to North Pacific and North Atlantic Ice
        nh_ice_atl = np.where((sea_ice_frac>=0.5)  & \
                              (lat > 0) & (lon>-101) & (lon<109) & \
                                  (np.isnan(biomes_all[i,:,:])==True))
            
        nh_ice_pac = np.where((sea_ice_frac>=0.5) &  (lat > 0) & \
                              ((lon<=-101) | (lon>=109)) & \
                              (np.isnan(biomes_all[i,:,:])==True))
            
        biomes_all[i, nh_ice_atl[0],nh_ice_atl[1]] = 8
        biomes_all[i, nh_ice_pac[0],nh_ice_pac[1]] = 1
        
        ## SUBPOLAR SEASONALLY STRATIFIED SPSS
        
        # North Hemisphere SPSS
        
        # # North Pacific SPSS
        nh_spss_pac = np.where((sst<14) & (chl>=0.25) & \
                               (lat>0) & ((lon<=-101) |  (lon>=109)) & \
                                   (np.isnan(biomes_all[i,:,:])==True))
                                   
        biomes_all[i, nh_spss_pac[0],nh_spss_pac[1]] = 2
        
        # # North Atlantic SPSS
        nh_spss_atl = np.where((sst<14) & (chl>=0.4) & \
                               (lat>0) & (lon>-101) & (lon<109) & \
                                   (np.isnan(biomes_all[i,:,:])==True))
            
        biomes_all[i, nh_spss_atl[0],nh_spss_atl[1]] = 9
        
        ## SUBTROPICAL SEASONALLY STRATIFIED STSS
        
        # North Pacific STSS
        nh_stss_pac = np.where((sst>=11) & (sst<29) & (lat>=25) & \
                               (((chl>=0.16) & (chl<0.4)) | (mld >125)) & \
                                   ((lon<=-101) |  (lon>=109)) & \
                                       (np.isnan(biomes_all[i,:,:])==True))
        biomes_all[i, nh_stss_pac[0],nh_stss_pac[1]] = 3
        
        # North Atlantic STSS
        nh_stss_atl = np.where((sst>=11) & (sst<29) & (lat>=25) & \
                               (((chl>=0.16) & (chl<0.4)) | (mld >125)) & \
                                   (lon>-101) & (lon<109) & \
                                       (np.isnan(biomes_all[i,:,:])==True))
            
        biomes_all[i, nh_stss_atl[0],nh_stss_atl[1]] = 10
        
        ## SUBTROPICAL PERMANENTLY STRATIFIED STPS
    
        # # North Pacific STPS
        nh_stps_pac=np.where((sst>=14) & (sst<29) & \
                              (chl<0.16) & (mld<=125) & \
                                  (((lon<=-83.8) | (lon>103)) & ((lat>0) & (lat<15))) | \
                                      (((lon<=-100) | (lon>103)) & ((lat>15) & (lat<90))) &\
                                  (np.isnan(biomes_all[i,:,:])==True))
            
        inds = np.where((((lon<=-83.8) | (lon>103)) & ((lat>0) & (lat<15))) | \
                         (((lon<=-90) | (lon>103)) & ((lat>15) & (lat<90))))
        biomes_all[i, nh_stps_pac[0],nh_stps_pac[1]] = 4
        
        # # North Atlantic STPS
        nh_stps_atl = np.where((sst>=14) & (sst<29) & \
                              (chl<0.16) & (mld<=125) & \
                              (lat>0) & (lon>-83.8) & (lon<13) & \
                                  (np.isnan(biomes_all[i,:,:])==True))
            
        biomes_all[i, nh_stps_atl[0],nh_stps_atl[1]] = 11
        
        ## EQUATORIAL EQU
        
        # West Pacific EQU
        equ_west_pac = np.where((sst>=29) & (chl<0.25) & \
                                (lat<=15) & (lat>=-15) & \
                                    ((lon>100) | (lon<-100)) & \
                                    (np.isnan(biomes_all[i,:,:])==True))
                                                # lon>100))
        biomes_all[i, equ_west_pac[0],equ_west_pac[1]] = 5
        
        # East Pacific EQU
        
        equ_east_pac = np.where((sst>=19) & (sst<29)  & \
                                (chl>=0.16) & (chl<0.7) & \
                                    (lat<=10) & (lat>=-10) & \
                                        ((lon>160) | (lon<-76.5)) & \
                                        (np.isnan(biomes_all[i,:,:])==True))
                                        
        biomes_all[i, equ_east_pac[0],equ_east_pac[1]] = 6
        
        # # Atlantic EQU
        equ_atl = np.where((sst>=19) & (sst<29) & \
                           (chl>=0.16) & (chl<0.7) & \
                               (lat<=10) & (lat>=-10) & \
                                   (lon>-76.5) & (lon<18) &\
                                       (np.isnan(biomes_all[i,:,:])==True))
            
        biomes_all[i, equ_atl[0],equ_atl[1]] = 12
        
        #########################
        ## SOUTHERN HEMISPHERE
        #########################
        
        ## ICE BIOMES
        sh_ice_inds = np.where((sea_ice_frac>=0.5) & (lat < 0) & \
                               (np.isnan(biomes_all[i,:,:])==True))
            
        biomes_all[i, sh_ice_inds[0],sh_ice_inds[1]] = 17
        
        ## SUBPOLAR SEASONALLY STRATIFIED SPSS
        
        sh_spss_ind = np.where((sst<8) & (lat<=0) & \
                               (np.isnan(biomes_all[i,:,:])==True))
            
        biomes_all[i, sh_spss_ind[0],sh_spss_ind[1]] = 16
        
        ## SUBTROPICAL SEASONALLY STRATIFIED STSS
        
        sh_stss_atl = np.where((sst>=8) & \
                                ((chl>=0.16) | (mld>150)) & \
                                    (lat<-35) & \
                                    (np.isnan(biomes_all[i,:,:])==True))
        biomes_all[i, sh_stss_atl[0],sh_stss_atl[1]] = 15
        
        ## SUBTROPICAL PERMANENTLY STRATIFIED STPS
        
        # Southern Hemisphere
        # South Pacific STPS
    
        sh_stps_pac = np.where((sst>=8) & (chl<0.25) & \
                              ( mld<=150) & (lat<=0) & \
                                  ((lon<=-69.5) | (lon>=150)) & \
                                  (np.isnan(biomes_all[i,:,:])==True))
        biomes_all[i, sh_stps_pac[0],sh_stps_pac[1]] = 7
        
        # # South Atlatic STPS
        sh_stps_atl=np.where((sst>=8) & (chl<0.25) & \
                             (mld<=150) & (lat<=0) & \
                                 (lon>=-69.5) & (lon<=22.7) & \
                                     (np.isnan(biomes_all[i,:,:])==True))
        biomes_all[i, sh_stps_atl[0],sh_stps_atl[1]] = 13
        
        # # Indian Ocean
        # Indian Ocean STPS
        ind_stps=np.where((sst>= 11) & (chl<0.25) & \
                          ((((lon>22.7) & (lon< 117)) & (lat<-10)) | \
                              (((lon>22.7) & (lon< 105)) & (lat>-10))) &\
                              (np.isnan(biomes_all[i,:,:])==True))
            
        biomes_all[i, ind_stps[0],ind_stps[1]] = 14
        
        ####
        ####
        
        ## Smooth data and gap fill
        ## Median filter
        
        unedited = biomes_all[i,:,:].copy()
        
        # Interpolate with nearest neighbor
        inds = np.where(np.isnan(unedited)==False)
        interp = NearestNDInterpolator(list(zip(lon[inds].flatten(), lat[inds].flatten())), unedited[inds].flatten())
        
        gapfilled = interp(lon, lat)
                           
        # Land mask
        gapfilled[np.where(landmask==True)] = np.NaN
        
        smoothed = scipy.ndimage.median_filter(gapfilled, size = 3)
        
        biomes_all[i,:,:]=smoothed
        
        # marginal seas, etc
        bad_inds = np.where(np.isnan(np.flipud(fm_biomes.MeanBiomes.values.T))==True)
        biomes_all[i,bad_inds[0],bad_inds[1]]=np.NaN
        
        # fig = plt.figure(figsize=(10,6))
        # rr = 2; cc = 4
        
        # param_list = [sea_ice_frac, sst, mld, chl]
        # name_list = ['Sea Ice','SST','MLD', 'CHL']
        # for ri in np.arange(cc):
            
        #     ax = fig.add_subplot(rr,cc,int(ri+1))
        #     ax.pcolormesh(lon, lat, param_list[ri])
        #     ax.set_title(name_list[ri])
            
        # param_list = [unedited, gapfilled, smoothed, biomes_all[i,:,:]]
        # name_list = ['Initial','Gap-Filled','Smoothed', 'Final']
        # for ri in np.arange(cc):
            
        #     ax = fig.add_subplot(rr,cc,int(cc+ri+1))
        #     ax.pcolormesh(lon, lat, param_list[ri], cmap = cmm)
        #     ax.set_title(name_list[ri])
            
        # fig.tight_layout()
        # if i == biomes_all.shape[0]-1:
        #     fig.suptitle('Mean Conditions 2010-'+str(int(end_year)))
        #     plt.savefig(figdir+'Mean.jpg', dpi = 300)
        # else:
        #     fig.suptitle(str(yr_range[i]))
        #     plt.savefig(figdir+str(yr_range[i])+'.jpg', dpi = 300)
        
        # plt.close()
        
        # plt.figure(figsize = (6.5, 4))
        # plt.pcolormesh(lon, lat, biomes_all[i,:,:], cmap = cmm)
        # plt.colorbar()
        
        # if i == biomes_all.shape[0]-1:
        #     plt.title('Mean Conditions 2010-'+str(int(end_year)))
        #     fig.tight_layout()
        #     plt.savefig(figdir+'Mean.jpg', dpi = 300)

        # else:
        #     plt.title(str(yr_range[i]))
        #     fig.tight_layout()
        #     plt.savefig(figdir+str(yr_range[i])+'.jpg',dpi = 300)
        # plt.close()
        
        
        if i == biomes_all.shape[0]-1:
            title_str = 'Mean Conditions 2010-'+str(int(end_year))
        else:
            title_str =str(yr_range[i])
            
        fig = plt.figure(figsize = (6.5, 8))
        gs = gridspec.GridSpec(2, 3, figure=fig, 
                               height_ratios=[4, 1], 
                               width_ratios=[1,1,1])
        
        axmap = fig.add_subplot(gs[0,:])
        ax1 = fig.add_subplot(gs[1,0])
        ax2 = fig.add_subplot(gs[1,1])
        ax3 = fig.add_subplot(gs[1,2])
        
        param_list = [unedited, gapfilled, smoothed, biomes_all[i,:,:]]
        name_list = ['Initial','Gap-Filled','Smoothed', 'Final\n'+title_str]
        ax_list = [ax1, ax2,ax3, axmap]
        for ri,ax in enumerate(ax_list):
            
            cax = ax.pcolormesh(lon, lat, param_list[ri], cmap = cmm)
            ax.set_title(name_list[ri])
            
            if ri == 0 or ri == 3:
                ax.set_ylabel('Longitude (ºN)')
                
            if ri == 3:
                cbar = fig.colorbar(cax,ax = ax,
                                    ticks=bnames.loc[:,'BIOME_NUM'].values,
                                    format=mticker.FixedFormatter(bnames.loc[:,'BIOME_ABREV'].values),
                                    location = 'bottom', alpha = 1)
                cbar.ax.tick_params(rotation=90)
                cbar.draw_all()
            
            ax.set_xlabel('Latitude (ºE)')
            
            
        fig.tight_layout()
        if i == biomes_all.shape[0]-1:
            plt.savefig(figdir+'Mean.jpg', dpi = 300)
        else:
            plt.savefig(figdir+str(yr_range[i])+'.jpg', dpi = 300)
        
        plt.close()
        
        # sval = 200
        # lon360 = np.concatenate((lon[:,sval:],lon[:,:sval]), axis=1)
        
        # landmask = np.concatenate((landmask[:,sval:],landmask[:,:sval]), axis=1)
        # smoothed = np.concatenate((smoothed[:,sval:],smoothed[:,:sval]), axis=1)
        # gapfilled = np.concatenate((gapfilled[:,sval:],gapfilled[:,:sval]), axis=1)
        # unedited2 = np.concatenate((unedited[:,sval:],unedited[:,:sval]), axis=1)
        
        # fig = plt.figure(figsize = (15,8))
        # ax1 = fig.add_subplot(1,3,1)
        # ax2 = fig.add_subplot(1,3,2, sharey = ax1, sharex=ax1)
        # ax3 = fig.add_subplot(1,3,3, sharey = ax1, sharex=ax1)
        
        # ax1.pcolormesh(lon, lat, np.where(landmask==True, landmask,np.NaN),cmap='Greys_r')
        # ax2.pcolormesh(lon, lat, np.where(landmask==True, landmask,np.NaN),cmap='Greys_r')
        # ax3.pcolormesh(lon, lat, np.where(landmask==True, landmask,np.NaN),cmap='Greys_r')
        
    
        # inds = np.where((unedited == 4) |(unedited == 5) |(unedited == 6) | (unedited == 7))
        # ax1.scatter(lon[inds], lat[inds], c=cmm(unedited[inds].astype(int)-1), s=20)
    
        
        # inds = np.where((gapfilled == 4) |(gapfilled == 5) |(gapfilled == 6)| (gapfilled == 7))
        # ax2.scatter(lon[inds], lat[inds], c=cmm(gapfilled[inds].astype(int)-1), s=20)
        
    
        # inds = np.where((biomes_all[i,:,:] == 4) |(biomes_all[i,:,:] == 5) |(biomes_all[i,:,:] == 6)| (biomes_all[i,:,:] == 7))
        # ax3.scatter(lon[inds], lat[inds], c=cmm(biomes_all[i,:,:][inds].astype(int)-1), s=20)
        
        # if i == biomes_all.shape[0]-1:
        #     fig.suptitle('MEAN')
        # else:
        #     fig.suptitle((str(yr_range[i])))
        # ax1.set_ylim([-30,30])
        #ax1.set_xlim([-180,-50])
        
    biome_ds = xr.Dataset(
                    data_vars=dict(
                        Biomes=(["x", "y","z"], biomes_all[:-1,:,:]),
                        Biomes_Change=(["x", "y","z"], change_all[:-1,:,:]),
                        MeanBiomes=(["y","z"], biomes_all[-1,:,:]),
                        MeanBiomes_Change=(["y","z"], change_all[-1,:,:]),
                    ),
                    coords=dict(
                        year=(["x"], yr_range),
                        lat=(["y"], lat_range),
                        lon=(["z"], lon_range),
                    ),
                    attrs=dict(description="Fay-McKinley Biomes 2010-"+str(end_year)),
                    )
    
    outfname = savedir + 'biomes_2010_'+str(int(end_year))+'.nc'
    
    biome_ds.to_netcdf(outfname)
    
    return


