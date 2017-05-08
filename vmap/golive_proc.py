#! /usr/bin/env python

import sys
import os
import numpy as np
from osgeo import gdal
from pygeotools.lib import iolib,timelib

#wget -r ftp://sidads.colorado.edu/pub/DATASETS/nsidc0710_landsat8_golive_ice_velocity_v1/p046_r027/*nc
#parallel 'gdal_translate NETCDF:"{}":vv_masked {.}_vv_masked.tif' ::: *nc
#parallel './golive_rename.py {}' ::: L8*vv.tif
#parallel './golive_rename.py {}' ::: L8*vv_masked.tif

#make_stack.py -te '577880 5175800 611800 5205050' 2*vv.tif
#make_stack.py -te '577880 5175800 611800 5205050' 2*vv_masked.tif
#dem_mosaic --median *vv_masked.tif -o vv_masked

in_fn = sys.argv[1]
#L8_046_027_048_2015_350_2016_033_v1_vv_masked.tif
a = in_fn.split('_')
dt1 = timelib.j2dt(int(a[4]), int(a[5]))
dt2 = timelib.j2dt(int(a[6]), int(a[7]))
dtc = timelib.center_date(dt1, dt2)
ndays = (dt2 - dt1).days
out_fn = '%s_%iday_%s-%s_%s' % (dtc.strftime('%Y%m%d'), ndays, dt1.strftime('%Y%m%d'), dt2.strftime('%Y%m%d'), in_fn) 
print(out_fn)
if os.path.exists(out_fn):
    os.remove(out_fn)
os.symlink(in_fn, out_fn)

b_ds = gdal.Open(out_fn)
b = iolib.ds_getma(b_ds)
#Threshold for Rainier
b_thresh = 1.5
b[b>b_thresh] = np.ma.masked
out_fn_masked = os.path.splitext(out_fn)[0]+'_filt.tif'
iolib.writeGTiff(b,out_fn_masked,b_ds)

