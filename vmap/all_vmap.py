#! /usr/bin/env python

"""
Create lists of viable velocity pairs
Pipe output to a text file and excecute

Before running, co-register all DEMs with dem_coreg_all.sh

Create 4-m subsampled versions from 2-m products
#Note: the 4-m subsampling can be done on the fly in vmap, best to do hs_multi on 2-m
#gdal_opt='-co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER'
#parallel "gdalwarp -overwrite $gdal_opt -tr 4 4 -r med {} {.}_4m.tif" ::: *00/dem*/*-DEM_2m_trans.tif

Process all shaded relief maps in current directory
#parallel 'hs_multi.sh {}' ::: *00/dem*/*DEM_2m_trans_4m.tif
parallel 'hs_multi.sh {}' ::: *00/dem*/*align/*DEM.tif

#Create lists of viable pairs
mkdir vmap; cd vmap
all_vmap.py ../*00/dem*/*align/*DEM_hs_multi.tif > all_vmap_pairs.txt
all_vmap.py ../*00/*16b.tif > all_vmap_pairs.txt

#Run vmap
parallel --progress --delay 0.5 -j 16 < all_vmap_pairs.txt

#Combine
post_vmap.sh
"""

import sys
from datetime import datetime, timedelta 
import numpy as np
from pygeotools.lib import timelib

#Set these for min and max time difference
#min_ts_diff = timedelta(days=60)
min_ts_diff = timedelta(days=14)
max_ts_diff = timedelta(days=365*2.5) 

fn_list = np.sort(np.array(sys.argv[1:]))
ts_list = np.array([timelib.fn_getdatetime(fn) for fn in fn_list])

#vmap_arg = ['-filter', '-threads 1', '-tr', '4']
#The bareground masking in HMA can incorrectly remove glacier surfaces
#vmap_arg = ['-filter', '-threads 1', '-tr', '4', '-mask_input', '-dt', 'day']
vmap_arg = ['-filter', '-threads 1', '-tr', '4', '-dt', 'day']
#vmap_arg = ['-filter', '-threads 1', '-tr', '1.0', '-dt', 'day']
#vmap_arg = ['-remove_offsets', '-filter', '-threads 1', '-tr', '4']

for i,fn in enumerate(fn_list):
    ts = ts_list[i]
    ts_c = np.array(ts_list[(i+1):])
    ts_diff = ts_c - ts
    idx = (ts_diff < max_ts_diff) & (ts_diff > min_ts_diff)
    for j in fn_list[(idx.nonzero()[0]+(i+1))]:
        print("vmap.py %s %s %s" % (' '.join(vmap_arg), fn, j))
