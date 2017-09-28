#! /usr/bin/env python

import sys
from datetime import datetime, timedelta
import numpy as np
from pygeotools.lib import timelib

"""
Create clean filenames with center date for velocity maps generated from WV DEMs
"""

#Rename all
#mkdir vm
#for i in */2*vm.tif; do ln -s ../$i vm/$(~/src/vmap/vmap/wv_cdate.py $i)_vm.tif; done
#parallel "fn=$(~/src/vmap/vmap/wv_cdate.py {}); ln -s ../{} vm/${fn}_vm.tif" ::: */2*vm.tif

import os

outdir="vm"
if not os.path.exists(outdir):
    os.makedirs(outdir)

fn = sys.argv[1]
dt_list = timelib.fn_getdatetime_list(fn)
dt_list = np.sort(dt_list)
dtmin = dt_list.min()
dtmax = dt_list.max()
c_date = timelib.mean_date(dt_list)
ndays = timelib.dt_ptp(dt_list)
ndays = timelib.dt_ptp((dtmin, dtmax))
nyears = ndays/365.25
#s = '%s_%04idays_%s-%s' % (c_date.strftime('%Y%m%d'), ndays, dtmin.strftime('%Y%m%d'), dtmax.strftime('%Y%m%d'))
#s = '%s_%s-%s_%0.2fyr' % (c_date.strftime('%Y%m%d'), dtmin.strftime('%Y%m%d'), dtmax.strftime('%Y%m%d'), nyears)
#s = '%s_%s-%s_%04iday' % (c_date.strftime('%Y%m%d'), dtmin.strftime('%Y%m%d'), dtmax.strftime('%Y%m%d'), ndays)
#Added %H%M here, as there were some inputs acquired on same days
s = '%s__%s-%s__%04iday' % (c_date.strftime('%Y%m%d_%H%M'), dtmin.strftime('%Y%m%d_%H%M'), dtmax.strftime('%Y%m%d_%H%M'), ndays)
print(s)
os.symlink(os.path.join("..", fn), os.path.join(outdir, s+'_vm.tif'))
