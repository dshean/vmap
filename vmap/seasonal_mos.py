#! /usr/bin/env python

"""
Create seasonal velocity composites for specified input relative date ranges
"""

import sys
import subprocess
from pygeotools.lib import timelib
from datetime import datetime, timedelta

fn_list = sys.argv[1:]
dt_list = [timelib.fn_getdatetime_list(fn) for fn in fn_list]

season='winter'
#season='summer'

#Define relative start/end dates for each season
if season == 'winter':
    dt_min = (10, 1)
    dt_max = (6, 30)
else:
    dt_min = (6, 1)
    dt_max = (10, 31)

#Dummy variable for year to compute day of year
yr = 2016
dt_min_doy = timelib.dt2doy(datetime(yr, *dt_min))
dt_max_doy = timelib.dt2doy(datetime(yr, *dt_max))
if dt_min_doy > dt_max_doy:
    max_delta = datetime(yr+1, *dt_max) - datetime(yr, *dt_min)
else:
    max_delta = datetime(yr, *dt_max) - datetime(yr, *dt_min)
min_delta = timedelta(14)

out_fn_list = []
for fn in fn_list:
    dtl = timelib.fn_getdatetime_list(fn)
    #This is for the filenames from wv_cdate.py
    #20160819_2159__20160724_2201-20160914_2158__0051day_vm.tif
    #dt1, dt2 = dtl[1:3]
    #This is for the filenames from disp2v.py
    #20170523_2215_10200100606CAF00_1020010061481F00-DEM_2m_trans_hs_multi__20170808_2200_1020010060548B00_102001006618CD00-DEM_2m_trans_hs_multi_vmap_4m_35px_spm1-F_smooth_vm.tif
    dt1, dt2 = dtl[0:2]
    delta = dt2 - dt1
    #dt1_doy = timelib.dt2doy(dt1)
    #dt2_doy = timelib.dt2doy(dt2)
    #dt_min_doy = timelib.dt2doy(datetime(dt1.year, *dt_max))
    #dt_max_doy = timelib.dt2doy(datetime(dt1.year, *dt_max))
    if delta <= max_delta and delta >= min_delta:
        #print(dt1, dt2)
        #print(timelib.rel_dt_test(dt1, dt_min, dt_max), timelib.rel_dt_test(dt2, dt_min, dt_max))
        if timelib.rel_dt_test(dt1, dt_min, dt_max) and timelib.rel_dt_test(dt2, dt_min, dt_max):
            if dt1.year == dt2.year and dt2 > datetime(dt2.year, *dt_max): 
                continue
            print(fn)
            out_fn_list.append(fn)

print(len(out_fn_list))

if True:
    out_prefix = 'vm_%s_%i' % (season, len(out_fn_list))
    cmd = ['dem_mosaic', '-o', out_prefix]
    cmd.append('--med')
    cmd.extend(out_fn_list)
    subprocess.call(cmd)
