#! /usr/bin/env python

import sys
from datetime import datetime, timedelta
import numpy as np
from pygeotools.lib import timelib

"""
Create clean filenames with center date for velocity maps generated from WV DEMs
"""

#Rename all
#for i in */2*vm.tif; do ln -s $i $(~/src/vmap/vmap/wv_cdate.py $i)_vm.tif; done

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
s = '%s_%s-%s_%04iday' % (c_date.strftime('%Y%m%d'), dtmin.strftime('%Y%m%d'), dtmax.strftime('%Y%m%d'), ndays)
print(s)
