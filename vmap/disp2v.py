#! /usr/bin/env python

#David Shean
#dshean@gmail.com

#Convert vmap output disparity map (units of pixels) to velocities (projected coordinates)
#Input disparity should have dates in filename and srs defined in GTiff header

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import timedelta

from pygeotools.lib import iolib
from pygeotools.lib import malib
from pygeotools.lib import geolib
from pygeotools.lib import warplib
from pygeotools.lib import timelib

def make_plot(m, fig_fn, label):
    f, ax = plt.subplots(figsize=(7,7))
    #plt.title('%s to %s' % (t1.strftime('%Y-%m-%d'), t2.strftime('%Y-%m-%d')))
    perc = malib.calcperc(m, (2,98))
    cmap = 'jet'
    imgplot = ax.imshow(m, cmap=cmap)
    imgplot.set_clim(*perc)
    imgplot.axes.get_xaxis().set_visible(False)
    imgplot.axes.get_yaxis().set_visible(False)
    imgplot.axes.patch.set_facecolor('black')
    cb = plt.colorbar(imgplot, orientation='vertical', extend='both', shrink=0.5)
    cb.set_label(label)
    return f, ax

def compute_stride(a):
    return stride

def plotvec(h, v):
    #Probably want to smooth h and v here, or at least filter outliers
    #Interpolation or hole filling would also help with lots of nodata holes 
    #Compute stats within window centered on X and Y
    nvec = np.array([20,20], dtype=np.float)
    stride = np.round(min(h.shape/nvec)).astype('int')
    X,Y = np.meshgrid(np.arange(0,h.shape[1],stride),np.arange(0,h.shape[0],stride))
    h_sub = h[::stride,::stride]
    v_sub = v[::stride,::stride]
    norm = False
    if norm:
        m_sub = np.ma.sqrt(h_sub**2+v_sub**2)
        h_sub /= m_sub
        v_sub /= m_sub
    Q = plt.quiver(X,Y,h_sub,v_sub,pivot='middle',color='white')
    #Add quiver key
    #scale = float(Q.scale)
    #lbl = str(np.round(Q.scale).astype(int)) + r'$\frac{m}{s}$'
    #lbl = '%i m/s' % np.around(scale, decimals=0).astype(int) 
    #qk = plt.quiverkey(Q, 0.05, 0.05, 1, lbl, labelpos='N', fontproperties={'weight':'bold'})

def main():
    if len(sys.argv) != 2:
        sys.exit("Usage is %s F.tif" % os.path.basename(sys.argv[0])) 

    #Output m/t_unit
    t_unit = 'year'
    #t_unit = 'day'

    #Generate plot of velocity magnitude with vectors overlaid
    plot = False

    #Remove offsets
    offset = False 

    #Mask defining static rock surfaces needed to remove horizontal/vertical offsets
    #If None, will attempt to use NLCD or global bare earth
    mask_fn = None

    #Input is 3-band disparity map, extract bands directly
    src_fn = sys.argv[1]
    src_ds = iolib.fn_getds(src_fn)
    #Extract pixel resolution
    h_res, v_res = geolib.get_res(src_ds)
    #Horizontal scale factor
    #If running on disparity_view output (gdal_translate -outsize 5% 5% F.tif F_5.tif)
    #h_res /= 20
    #v_res /= 20
   
    #Load horizontal and vertical disparities
    h = iolib.ds_getma(src_ds, bnum=1) 
    v = iolib.ds_getma(src_ds, bnum=2) 
    #ASP output has northward motion as negative values in band 2
    v *= -1

    t1, t2 = timelib.fn_getdatetime_list(src_fn)
    dt = t2 - t1
    #Default t_factor is in 1/years
    t_factor = timelib.get_t_factor(t1,t2)
    #Input timestamp arrays if inputs are mosaics
    t1_fn = ''
    t2_fn = ''
    #t1_fn, t2_fn = sys.argv[2:4]
    if os.path.exists(t1_fn) and os.path.exists(t2_fn):
        t_factor = timelib.get_t_factor_fn(t1_fn, t2_fn)
    if t_factor is None:
        sys.exit("Unable to determine input timestamps")

    if t_unit == 'day':
        t_factor *= 365.25

    print("Input dates:")
    print(t1)
    print(t2)
    print(dt)
    print(t_factor, t_unit)

    #Scale values for polar stereographic distortion
    srs = geolib.get_ds_srs(src_ds)
    proj_scale_factor = 1.0
    #Want to scale to get correct distances for polar sterographic
    if srs.IsSame(geolib.nps_srs) or srs.IsSame(geolib.sps_srs):
        proj_scale_factor = geolib.scale_ps_ds(src_ds) 

    #Convert disparity values in pixels to m/t_unit
    h_myr = h*h_res*proj_scale_factor/t_factor
    h = None
    v_myr = v*v_res*proj_scale_factor/t_factor
    v = None

    #Velocity Magnitude
    m = np.ma.sqrt(h_myr**2+v_myr**2)
    print("Velocity Magnitude stats")
    malib.print_stats(m)

    offset_str = ''
    #Calibrate over static surfaces
    if offset:
        #By default, use NLCD, otherwise global bare earth
        if mask_fn is None:
            from demcoreg.dem_mask import get_lulc_mask
            #Match disparity ds
            #Can get ds first, check if empty before proceeding
            #if not geolib.ds_IsEmpty(lulc_ds):
            #   print("Unable to compute offsets over static control surfaces")
            print("\nPreparing LULC mask to identify static control surfaces\n")
            mask = get_lulc_mask(src_ds, mask_glaciers=True)
        else:
            mask = iolib.fn_getma(mask_fn)
        print("\nRemoving median x and y offset over static control surfaces")
        h_myr_count = h_myr.count()
        h_myr_static_count = h_myr[mask].count()
        h_myr_med = malib.fast_median(h_myr[mask])
        v_myr_med = malib.fast_median(v_myr[mask])
        h_myr_mad = malib.mad(h_myr[mask])
        v_myr_mad = malib.mad(v_myr[mask])
        print("Static pixel count: %i (%0.1f%%)" % (h_myr_static_count, 100*float(h_myr_static_count)/h_myr_count))
        print("median (+/-NMAD)")
        print("Horizontal: %0.2f (+/-%0.2f) m/%s" % (h_myr_med, h_myr_mad, t_unit))
        print("Vertical: %0.2f (+/-%0.2f) m/%s" % (v_myr_med, v_myr_mad, t_unit))
        h_myr -= h_myr_med
        v_myr -= v_myr_med 
        offset_str = '_corr_h%0.2fma_v%0.2fma' % (h_myr_med, v_myr_med)
        #Velocity Magnitude
        m = np.ma.sqrt(h_myr**2+v_myr**2)
        print("Velocity Magnitude stats after correction")
        malib.print_stats(m)

    if plot:
        fig_fn = os.path.splitext(src_fn)[0]+'.png'
        label='Velocity (m/%s)' % t_unit
        f, ax = make_plot(m, fig_fn, label)
        plotvec(h_myr, v_myr)
        plt.tight_layout()
        plt.savefig(fig_fn, dpi=300, bbox_inches='tight', pad_inches=0, edgecolor='none')

    print("Writing out files") 
    gt = src_ds.GetGeoTransform()
    proj = src_ds.GetProjection()
    dst_fn = os.path.splitext(src_fn)[0]+'_vm%s.tif' % offset_str
    iolib.writeGTiff(m, dst_fn, create=True, gt=gt, proj=proj)
    dst_fn = os.path.splitext(src_fn)[0]+'_vx.tif'
    iolib.writeGTiff(h_myr, dst_fn, create=True, gt=gt, proj=proj)
    dst_fn = os.path.splitext(src_fn)[0]+'_vy.tif'
    iolib.writeGTiff(v_myr, dst_fn, create=True, gt=gt, proj=proj)
    src_ds = None

if __name__ == "__main__":
    main()
