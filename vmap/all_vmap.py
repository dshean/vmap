#! /usr/bin/env python

"""
Create lists of viable velocity pairs
Pipe output to a text file and excecute

Process all shaded relief maps in current directory
all_vmap.py *DEM_8m_hs_multi.tif > all_vmap_pairs.txt
parallel -j 16 < all_vmap_pairs.txt
"""

import sys
from datetime import datetime, timedelta 
import numpy as np
from pygeotools.lib import timelib

#Set these for min and max time difference
min_ts_diff = timedelta(days=60)
max_ts_diff = timedelta(days=365*3) 

fn_list = np.sort(np.array(sys.argv[1:]))
ts_list = np.array([timelib.fn_getdatetime(fn) for fn in fn_list])

for i,fn in enumerate(fn_list):
    ts = ts_list[i]
    ts_c = np.array(ts_list[(i+1):])
    ts_diff = ts_c - ts
    idx = (ts_diff < max_ts_diff) & (ts_diff > min_ts_diff)
    for j in fn_list[(idx.nonzero()[0]+(i+1))]:
        print("vmap.py %s %s" % (fn, j))
