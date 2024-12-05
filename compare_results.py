#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 14:59:22 2024

@author: epark
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from makebiome_fxns import GetBiomeColorMap

fay = xr.open_dataset('Fay_McKinley_Time_Varying_Biomes.nc')
biomes = xr.open_dataset('/Volumes/FATDATABABY/faymckinleybiomes/biomes_2010_2023.nc')

fay.close()
biomes.close()

cmm = GetBiomeColorMap()

lon, lat = np.meshgrid(biomes.lon.values, biomes.lat.values)
dif = np.flipud(fay.TimeVaryingBiomes.values[:,:,-1].T)-biomes.Biomes.values[0,:,:]

fig = plt.figure(figsize = (6.5,8))
rr = 2; cc = 1

ax1 = fig.add_subplot(rr,cc,1)
ax1.pcolormesh(lon,lat, np.flipud(fay.TimeVaryingBiomes.values[:,:,-1].T), cmap = cmm)
ax1.set_title('Fay and McKinley')

ax2 = fig.add_subplot(rr,cc,2)
ax2.pcolormesh(lon,lat, biomes.Biomes.values[0,:,:], cmap = cmm)
ax2.set_title('This study')
# ax3 = fig.add_subplot(rr,cc,3)
# dif = np.flipud(fay.TimeVaryingBiomes.values[:,:,-1].T)-biomes.Biomes.values[0,:,:]
# ax3.pcolormesh(lon,lat, dif,
#                cmap = 'bwr', vmin=-1, vmax=1)

change_inds = np.where((np.isnan(dif)==False) & (dif!=0))
change = change_inds[0].shape[0]
total = np.where(np.isnan(biomes.Biomes.values[0,:,:])==False)[0].shape[0]

ts = 'Change Pixels: '+str(change)+'\n'+'% Difference: '+str(np.round(change/total*100,1))
for ax in [ax1,ax2]:
    ax.scatter(lon[change_inds], lat[change_inds],s=1,alpha = 0.3,
               c= 'k', label  =ts)

ax1.legend()
fig.tight_layout()
fig.savefig('figures/compare.jpg',dpi = 300)