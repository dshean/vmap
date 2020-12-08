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
from demcoreg import dem_mask
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

def ovr_resample(src_fn,res_factor,src_ds=None):
    """
    build overviews using Gaussian resampling method and returned subsampled masked arrays and dataset at a given scale
    #https://stackoverflow.com/questions/33158526/how-to-correctly-use-gdaladdo-in-a-python-program
    # glued together from several pygeotools functions

    Parameters
    -----------
    src_ds: gdal.Dataset
        path to gdal dataset
    src_fn: str
        path to ASP 3-band disparity map
    res_factor: int
        final scale at which subsampling is desired
    Returns
    -----------
    h_sub: np.ma.array
        subsampled horizontal disparities
    v_sub: np.ma.array
        subsampled vertical disparities
    src_ds: gdal.Dataset
        dataset at sampled at specified scale
    """
    # need to create padded dataset to handle dataset with odd dimensions
    # prevents sub-pixel shifts due to resampling
    if src_ds is None:
        src_ds = iolib.fn_getds(src_fn)
    ns = src_ds.RasterXSize
    nl = src_ds.RasterYSize
    # now pad with max overview scale
    p_olist = np.power(2,np.arange(1,10))
    olist = list(p_olist[p_olist<=res_factor])
    ns_pad = int((np.floor(ns/res_factor)+1)*res_factor)
    nl_pad = int((np.floor(nl/res_factor)+1)*res_factor)
    dtype = src_ds.GetRasterBand(1).DataType
    ns_pad_diff = ns_pad - ns
    nl_pad_diff = nl_pad - nl
    xres,yres = geolib.get_res(src_ds)
    # need to maintain geotransform so pad extent properly
    xshift = ns_pad_diff*xres # this will be added to xmax
    yshift = nl_pad_diff*yres # this will be subtracted from ymin
    init_extent = geolib.ds_extent(src_ds)
    final_extent = [init_extent[0],init_extent[1]-yshift,init_extent[2]+xshift,init_extent[3]]
    projwin_extent = [final_extent[0],final_extent[3],final_extent[2],final_extent[1]]
    # get padded in-memory vrt with gdal_translate
    # due to a bug, need to write vrt to disk using the src_fn file path and then re-read the vrt to get overviews
    # see here: # http://osgeo-org.1560.x6.nabble.com/gdal-dev-Python-bindings-BuildOverviews-not-supported-for-VRT-dataset-td5429453.html
    # Latest gdal should have this fixed
    tgt_vrt = os.path.splitext(src_fn)[0]+'_gdal.vrt'
    mem_vrt = gdal.Translate(tgt_vrt,src_fn,width=ns_pad,height=nl_pad,format='VRT',projWin = projwin_extent)
    # close the mem_vrt dataset and then re-read
    mem_vrt = None
    mem_vrt = iolib.fn_getds(tgt_vrt)
    # now we build the overviews
    gdal.SetConfigOption('COMPRESS_OVERVIEW','LZW')
    gdal.SetConfigOption('BIGTIFF_OVERVIEW','YES')
    print("Building Overviews")
    mem_vrt.BuildOverviews("gauss",olist)
    # now read the bands and the corresponding dataset
    h_sub,src_ds = iolib.ds_getma_sub(mem_vrt,bnum=1,scale=res_factor,return_ds=True)
    v_sub = iolib.ds_getma_sub(mem_vrt,bnum=2,scale=res_factor)
    mem_vrt = None
    return h_sub,v_sub,src_ds
    
 
def getparser():
    parser = argparse.ArgumentParser(description="Convert ASP disparity map to velocity map(s)")
    parser.add_argument('-dt', type=str, choices=['yr','day'], default='yr', help='Time increment (default: %(default)s)')
    parser.add_argument('-remove_offsets', action='store_true', help='Remove median offset from stable control surfaces, requires demcoreg or input mask')
    parser.add_argument('-mask_fn', type=str, default=None, help='Provide existing raster or vector mask for offset correction. If None, will automatically determine valid control surfaces (excluding glaciers, vegetation, etc.). This should work for any supported GDAL/OGR formats. For raster input, expects integer pixel values of 1 for valid control surfaces, 0 for surfaces to ignore during offset determination.')
    parser.add_argument('-plot', action='store_true', help='Generate plot of velocity magnitude with vectors overlaid')
    final_res_choices = [1,2,4,8]
    parser.add_argument('-final_res_factor',type=int,default=1,choices=final_res_choices,
                        help='Factor (Overview level) of input resultion (-tr) at which final velocity will be posted (default: %(default)s)')
    parser.add_argument('disp_fn', type=str, help='Disparity map filename (e.g., *-RD.tif, *-F.tif)')
    parser.add_argument('-mask_list', nargs='*', type=str, default=['glaciers'], choices=dem_mask.mask_choices, \
            help='Define masks to use to limit stable control surfaces for -remove_offsets option\
                 (**NOTE**: If provided, mask_list option should be specified as last input)')
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
        # the disparity map is in units of pixel, calculated at native resolution
        # when converting to displacement, the unit there should remain the same
        # so we should not be applying a res_factor to resolution
        #h_res = h_res*res_factor
        #v_res = v_res*res_factor
        # Build  gaussian overviews and read arrays,datset at specified res_factor
        # src_ds is swapped at this stage
        h,v,src_ds = ovr_resample(src_fn,res_factor,src_ds)
        # string for future filename saving
        res_factor_string=f'_{res_factor}x'

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
            print("\nUsing demcoreg to prepare mask of stable control surfaces\n")
            #mask = get_lulc_mask(src_ds, mask_glaciers=True, filter='rock+ice+water')
            # just have glaciers for now
            #TODO:
            # Probably best to predefine mask over static low sloped area and pass via mask_fn option
            mask_list = args.mask_list
            mask = dem_mask.get_mask(src_ds,mask_list=['glaciers'],dem_fn=src_fn)
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
    # src_ds is now updated if low res overview
    gt = src_ds.GetGeoTransform()
    proj = src_ds.GetProjection()
    dst_fn = f'{os.path.splitext(src_fn)[0]}_vm{res_factor_string}{offset_str}.tif'
    iolib.writeGTiff(m, dst_fn, create=True, gt=gt, proj=proj)
    dst_fn = f'{os.path.splitext(src_fn)[0]}_vx{res_factor_string}{offset_str}.tif'
    iolib.writeGTiff(h_myr, dst_fn, create=True, gt=gt, proj=proj)
    dst_fn = f'{os.path.splitext(src_fn)[0]}_vy{res_factor_string}{offset_str}.tif'
    iolib.writeGTiff(v_myr, dst_fn, create=True, gt=gt, proj=proj)
    src_ds = None

if __name__ == "__main__":
    main()
