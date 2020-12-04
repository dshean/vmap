#! /usr/bin/env python

#David Shean
#dshean@gmail.com

#Convert vmap output disparity map (units of pixels) to velocities (projected coordinates)
#Input disparity should have dates in filename and srs defined in GTiff header

import os
import sys
import argparse
from datetime import timedelta

import numpy as np
import matplotlib.pyplot as plt

from pygeotools.lib import iolib, malib, geolib, warplib, timelib
import gdal

def make_plot(m, fig_fn, label):
    f, ax = plt.subplots(figsize=(7,7))
    #plt.title('%s to %s' % (t1.strftime('%Y-%m-%d'), t2.strftime('%Y-%m-%d')))
    perc = malib.calcperc(m, (2,98))
    cmap = 'inferno'
    imgplot = ax.imshow(m, cmap=cmap)
    imgplot.set_clim(*perc)
    imgplot.axes.get_xaxis().set_visible(False)
    imgplot.axes.get_yaxis().set_visible(False)
    imgplot.axes.patch.set_facecolor('0.5')
    cb = plt.colorbar(imgplot, orientation='vertical', extend='both', shrink=0.5)
    cb.set_label(label)
    return f, ax

#Streamline plot
def streamline_plt(u,v,bg=None):
    """
    Generate a stremline plot for velocity field
    This needs further development and testing
    """
    dx = 1
    m = np.ma.sqrt(u*u + v*v)
    x,y = np.meshgrid(np.arange(0,u.shape[1],dx),np.arange(0,u.shape[0],dx))
    if bg is None:
        bg = m
    plt.imshow(bg)
    #plt.streamplot(x,y,u,v,color='k', linewidth=5*m/m.max())
    plt.streamplot(x,y,u,v,color='k')
    plt.show()

def compute_stride(a):
    return stride

def plotvec(h, v):
    #Probably want to smooth h and v here, or at least filter outliers
    #Interpolation or hole filling would also help with lots of nodata holes 
    #Compute stats within window centered on X and Y
    nvec = np.array([40,40], dtype=np.float)
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

def build_gauss_overview(fn,side_carr=True):
    """
    build overviews using Gaussian resampling method
    #https://stackoverflow.com/questions/33158526/how-to-correctly-use-gdaladdo-in-a-python-program
    Parameters
    -----------
    fn: str
        filepath
    side_carr: bool
        If true, overviews will be saved in *.ovr files
    """
    if side_carr:
        open_mode = 0 # open dataset in readonly mode
    else:
        open_mode = 1 # open dataset in write, overviews will be attached to dataset itself
    img_ds = gdal.Open(fn,open_mode)
    
    gdal.SetConfigOption('COMPRESS_OVERVIEW=LZW','BIGTIFF_OVERVIEW=YES')
    print("Building Overviews")
    img_ds.BuildOverviews("gauss",[2,4,8])
    del img_ds
    
def read_overviews(ds,bnum=1,ovr_level=1):
    """
    Function to read specied overview level from a gdal dataset
    Parameters
    -----------
    ds: gdal dataset
        dataset whose overview needs to be read
    bnum: int
        band number for which overview needs to be read
    ovr_level: int
        overview level to be read
    Returns
    -----------
    ma_ovr: np.ma array
        masked array of the read overview
    """
    # read the corresponding band
    b = ds.GetRasterBand(bnum)
    b_ndv = geolib.get_ndv_b(b)
    # read overview into array
    ovr_array = b.GetOverviews(ovr_level)
    # mask for no-data values
    ma_ovr = np.ma.masked_values(ovr_array,b_ndv)
    return ma_ovr

    
def getparser():
    parser = argparse.ArgumentParser(description="Convert ASP disparity map to velocity map(s)")
    parser.add_argument('-dt', type=str, choices=['yr','day'], default='yr', help='Time increment (default: %(default)s)')
    parser.add_argument('-remove_offsets', action='store_true', help='Remove median offset from stable control surfaces, requires demcoreg or input mask')
    parser.add_argument('-mask_fn', type=str, default=None, help='Provide existing raster or vector mask for offset correction. If None, will automatically determine valid control surfaces (excluding glaciers, vegetation, etc.). This should work for any supported GDAL/OGR formats. For raster input, expects integer pixel values of 1 for valid control surfaces, 0 for surfaces to ignore during offset determination.')
    parser.add_argument('-plot', action='store_true', help='Generate plot of velocity magnitude with vectors overlaid')
    final_res_choices = [1,2,4,8]
    parser.add_argument('-final_res_factor',type='int',default=1,choices=final_res_choices,
                        help='Factor (Overview level) of input resultion (-tr) at which final velocity will be posted (default: %(default)s)')
    parser.add_argument('disp_fn', type=str, help='Disparity map filename (e.g., *-RD.tif, *-F.tif)')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    t_unit = args.dt
    plot = args.plot 
    remove_offsets = args.remove_offsets 
    mask_fn = args.mask_fn
    if mask_fn is not None:
        remove_offsets = True

    #Input is 3-band disparity map, extract bands directly
    src_fn = args.disp_fn 
    if not iolib.fn_check(src_fn):
        sys.exit("Unable to locate input file: %s" % src_fn)

    src_ds = iolib.fn_getds(src_fn)
    if src_ds.RasterCount != 3:
        sys.exit("Input file must be ASP disparity map (3 bands: x, y, mask)")
    #Extract pixel resolution
    h_res, v_res = geolib.get_res(src_ds)

    #Horizontal scale factor
    #If running on disparity_view output (gdal_translate -outsize 5% 5% F.tif F_5.tif)
    #h_res /= 20
    #v_res /= 20

    # adding if clauses for overview reads
    # current options include 2,4,8 subsampling, corresponding to gaussian overview levels of 1,2,3
    # would need to scale resolution accoringly
    res_factor = args.final_res_factor

    if res_factor>1:
        h_res = h_res*res_factor
        v_res = v_res*res_factor
        # Build gaussian overviews
        build_gauss_overview(src_fn)
        # determine required_overview level
        possible_ovr_levels = np.array([1,2,3])
        possible_res_levels = np.power(2,possible_ovr_levels)
        # add 1 to account for numpy indexing
        over_idx = np.where(possible_res_levels==res_factor)[0]+1
        # finally read the two arrays at the given overview levels
        h = read_overviews(src_ds,1,over_idx)
        v = read_overviews(src_ds,2,over_idx)
        # string for future filename saving
        res_factor_string=f'_{res_factor}x_'
 
    else:
        # read the arrays at full resolution of correlation
        #Load horizontal and vertical disparities
        h = iolib.ds_getma(src_ds, bnum=1) 
        v = iolib.ds_getma(src_ds, bnum=2) 
        res_factor_string = ''

    #ASP output has northward motion as negative values in band 2
    v *= -1

    t1, t2 = timelib.fn_getdatetime_list(src_fn)
    dt = t2 - t1
    #Default t_factor is in 1/years
    t_factor = timelib.get_t_factor(t1,t2)

    #Input timestamp arrays if inputs are mosaics
    if False:
        t1_fn = ''
        t2_fn = ''
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

    #Remove x and y offsets over control surfaces
    offset_str = ''
    if remove_offsets:
        if mask_fn is None:
            from demcoreg.dem_mask import get_mask
            print("\nUsing demcoreg to prepare mask of stable control surfaces\n")
            #mask = get_lulc_mask(src_ds, mask_glaciers=True, filter='rock+ice+water')
            # just have glaciers for now
            #TODO:
            # Extend to LULC or barground based on raster extent (CONUS or not)
            mask = get_mask(src_ds,mask_list=['glaciers'],dem_fn=src_fn)
        else:
            print("\nWarping input raster mask")
            #This can be from previous dem_mask.py run (e.g. *rockmask.tif)
            mask_ds = warplib.memwarp_multi_fn([mask_fn,], res=src_ds, extent=src_ds, t_srs=src_ds)[0]
            mask = iolib.ds_getma(mask_ds)
            #The default from ds_getma is a masked array, so need to isolate boolean mask
            #Assume input is 0 for masked, 1 for unmasked (valid control surface)
            mask = mask.filled().astype('bool')
            #This should work, as the *rockmask.py is 1 for unmasked, 0 for masked, with ndv=0
            #mask = np.ma.getmaskarray(mask)
            #Vector mask - untested
            if os.path.splitext(mask_fn)[1] == 'shp':
                mask = geolib.shp2array(mask_fn, src_ds)

        print("\nRemoving median x and y offset over static control surfaces")
        h_myr_count = h_myr.count()
        h_myr_static_count = h_myr[mask].count()
        h_myr_med = malib.fast_median(h_myr[mask])
        v_myr_med = malib.fast_median(v_myr[mask])
        h_myr_mad = malib.mad(h_myr[mask])
        v_myr_mad = malib.mad(v_myr[mask])
        print("Static pixel count: %i (%0.1f%%)" % (h_myr_static_count, 100*float(h_myr_static_count)/h_myr_count))
        print("median (+/-NMAD)")
        print("x velocity offset: %0.2f (+/-%0.2f) m/%s" % (h_myr_med, h_myr_mad, t_unit))
        print("y velocity offset: %0.2f (+/-%0.2f) m/%s" % (v_myr_med, v_myr_mad, t_unit))
        h_myr -= h_myr_med
        v_myr -= v_myr_med 
        offset_str = '_offsetcorr_h%0.2f_v%0.2f' % (h_myr_med, v_myr_med)
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
    # the geotranform would need to be adjusted for overview
    gt = src_ds.GetGeoTransform()
    proj = src_ds.GetProjection()
    dst_fn = os.path.splitext(src_fn)[0]+'_vm%s.tif' % offset_str
    iolib.writeGTiff(m, dst_fn, create=True, gt=gt, proj=proj)
    dst_fn = os.path.splitext(src_fn)[0]+'_vx%s.tif' % offset_str
    iolib.writeGTiff(h_myr, dst_fn, create=True, gt=gt, proj=proj)
    dst_fn = os.path.splitext(src_fn)[0]+'_vy%s.tif' % offset_str
    iolib.writeGTiff(v_myr, dst_fn, create=True, gt=gt, proj=proj)
    src_ds = None

if __name__ == "__main__":
    main()
