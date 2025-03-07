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
from set_paths import outputdir
import matplotlib.ticker as mticker
import pandas as pd
from matplotlib import gridspec

end_year = 2024
fay = xr.open_dataset('Fay_McKinley_Time_Varying_Biomes.nc')
biomes = xr.open_dataset(outputdir+'faymckinleybiomes/biomes_2010_'+str(int(end_year))+'.nc')

fay.close()
biomes.close()

cmm = GetBiomeColorMap()
bnames = pd.read_csv('biomes_Fay.csv')

alpha = 0.1

lon, lat = np.meshgrid(biomes.lon.values, biomes.lat.values)
dif = np.flipud(fay.TimeVaryingBiomes.values[:,:,-1].T)-biomes.Biomes.values[0,:,:]

fig = plt.figure(figsize = (4, 6))
gs = gridspec.GridSpec(2,1, figure=fig, 
                       height_ratios=[3, 4], 
                       )

ax1 = fig.add_subplot(gs[0,0])
ax1.pcolormesh(lon,lat, np.flipud(fay.TimeVaryingBiomes.values[:,:,-1].T), cmap = cmm)
ax1.set_title('Fay and McKinley')

ax2 = fig.add_subplot(gs[1,0])
cax = ax2.pcolormesh(lon,lat, biomes.Biomes.values[0,:,:], cmap = cmm)
ax2.set_title('This work')
# ax3 = fig.add_subplot(rr,cc,3)
# dif = np.flipud(fay.TimeVaryingBiomes.values[:,:,-1].T)-biomes.Biomes.values[0,:,:]
# ax3.pcolormesh(lon,lat, dif,
#                cmap = 'bwr', vmin=-1, vmax=1)

change_inds = np.where((np.isnan(dif)==False) & (dif!=0))
change = change_inds[0].shape[0]
total = np.where(np.isnan(biomes.Biomes.values[0,:,:])==False)[0].shape[0]

ts = 'Change Pixels: '+str(change)+'\n'+'% Difference: '+str(np.round(change/total*100,1))
for ax in [ax1,ax2]:
    ax.scatter(lon[change_inds], lat[change_inds],s=1,alpha = alpha,facecolor = 'None',
               edgecolor = 'k', label  =ts)

ax1.legend()
cbar = fig.colorbar(cax,ax = ax2,
                    ticks=bnames.loc[:,'BIOME_NUM'].values,
                    format=mticker.FixedFormatter(bnames.loc[:,'BIOME_ABREV'].values),
                    location = 'bottom', alpha = 1)
cbar.ax.tick_params(rotation=90)
cbar.draw_all()

fig.tight_layout()
fig.savefig(outputdir+'faymckinleybiomes/figures/compare.jpg',dpi = 300)