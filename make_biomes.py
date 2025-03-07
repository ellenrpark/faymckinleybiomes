#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 11:31:24 2024

@author: epark
"""

from downloaddata_fxns import DownloadData
from reformatdata_fxns import ReformatData
from makebiome_fxns import CalculateBiomes
import os
import shutil
from set_paths import outputdir


end_year = 2024

overwrite_data = False
download_data = False

# Set up directory for data downlaod and output
# Create dir to download data and output final file
if os.path.exists(outputdir+'faymckinleybiomes/') == False:
    
    # Main directory
    os.mkdir(outputdir+'faymckinleybiomes/')
    
    # Where downloaded data will be stored
    os.mkdir(outputdir+'faymckinleybiomes/data/')
    os.mkdir(outputdir+'faymckinleybiomes/data/raw/')
    
elif os.path.exists(outputdir+'faymckinleybiomes/') and (overwrite_data):
    
    # Clear data folders
    shutil.rmtree(outputdir+'faymckinleybiomes/data/raw/')
    os.mkdir(outputdir+'faymckinleybiomes/data/raw/')
      

savedir = outputdir+'faymckinleybiomes/'
datadir = outputdir+'faymckinleybiomes/data/'

# If you want to re-download data
if download_data:
    DownloadData(datadir)
    
# Reformat data
ice, sst, mld, chl = ReformatData(datadir, end_year)

# Calculate biomes
biomes = CalculateBiomes(savedir, end_year, ice, sst, mld, chl)
