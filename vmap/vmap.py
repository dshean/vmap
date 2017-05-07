#! /usr/bin/env python

#David Shean
#dshean@gmail.com

#This script uses ASP correlator to produce disparity maps from two inputs
#Input data should be orthorectified/mapped in the same projected coordinate system
#Run disp2v.py to convert to surface velocities

import sys
import os
import subprocess
from multiprocessing import cpu_count
from datetime import datetime, timedelta
from distutils.spawn import find_executable

from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt

from pygeotools.lib import warplib
from pygeotools.lib import geolib
from pygeotools.lib import iolib
from pygeotools.lib.malib import calcperc 
from pygeotools.lib.timelib import get_t_factor_fn

#Generate and execute stereo commands
def run_cmd(bin, args, **kw):
    #Note, need to add full executable
    binpath = find_executable(bin)
    call = [binpath]
    call.extend(args)
    print(' '.join(call))
    try:
        code = subprocess.call(call, shell=False)
    except OSError as e:
        raise Exception('%s: %s' % (binpath, e))
    if code != 0:
        raise Exception('Stereo step ' + kw['msg'] + ' failed')

def get_stereo_opt(maxnthreads=28, kernel=(21,21), nlevels=5, spr=1, timeout=360, erode=0):
    stereo_opt = []

    #This is irrelevant
    stereo_opt.extend(['-t', 'pinhole'])
    
    #Set number of threads/cores to use
    #Use all CPUs if possible
    nthreads = cpu_count() 
    if nthreads > maxnthreads:
        nthreads = maxnthreads
    stereo_opt.extend(['--threads', str(nthreads)])

    #This assumes input images are already mapped 
    stereo_opt.extend(['--alignment-method', 'None'])

    #This should be explored further
    stereo_opt.append('--individually-normalize')

    #Set correlator kernel size
    #Smaller kernels will offer more detail but are prone to more noise
    #No triangulation filters here

    #Because surface is changing between image interval, using larger kernel can help
    #Good choices are 11, 21 (ASP default), 35, 71 
    stereo_opt.extend(['--corr-kernel', str(kernel[0]), str(kernel[1])])
   
    #Numer of gaussian pyramids to use
    #Can look at texture in GDAL overviews to make a decision
    #If you can see plenty of texture at 1/32 resolution, go with 5 
    #For featureless areas, limiting to 2 can help, or even 0
    stereo_opt.extend(['--corr-max-levels', str(nlevels)])
  
    if timeout > 0:
        stereo_opt.extend(['--corr-timeout', str(timeout)])

    #Define the search area
    #Useful if you know your orhotorectification is good to, say 100 pixels in any direction
    #stereo_opt.extend(['--corr-search', '-100', '-100', '100', '100'])

    #Sub-pixel refinement
    #1)Parabolic, 2)Bayes, 3)AffineAdaptive
    #See ASP doc or Shean et al, ISPRS, (2016)
    #1 is fast but lower quality, 2 is slow but highest quality, 
    #3 is a good compromise for speed and quality
    stereo_opt.extend(['--subpixel-mode', str(spr)])

    #Sub-pixel kernel size
    #ASP default is 35
    stereo_opt.extend(['--subpixel-kernel', str(kernel[0]), str(kernel[1])])

    #Note: stereo_fltr throws out a lot of good data when noisy
    #Want to play with the following options
    #--rm-half-kernel 5 5
    #--rm_min_matches 60
    #--rm-threshold 3

    if erode > 0:
        stereo_opt.extend(['--erode-max-size', str(erode)])

    return stereo_opt

def make_ln(outdir, outprefix, ext):
    #Create symbolic links with appropriate names 
    ln_fn = os.path.join(outdir, outdir+ext)
    if os.path.lexists(ln_fn):
        os.remove(ln_fn)
    os.symlink(os.path.split(outprefix)[1]+ext, ln_fn)
    return ln_fn

#Create a dummy camera model file, so we don't have to find it
#This may no longer be necessary, but was a clever hack 
#Same as StereoPipeline/data/K10/black_left.tsai
def dummy_tsai():
    import tempfile
    fn = None
    with tempfile.NamedTemporaryFile(suffix='.tsai',delete=False) as temp:
        temp.write("""\
fu = 620.857971191406
fv = 622.298034667969
cu = 526.128784179688
cv = 374.273742675781
u_direction = 1 0 0
v_direction = 0 1 0
w_direction = 0 0 1
C = 0.0751 -0.065 -0.902
R = 0.000794376506820371 -0.371711462117335 0.928347972420126 0.997550683961614 -0.0646367856775087 -0.0267342264708707 0.0699428473375328 0.92609539188351 0.370749677327501
k1 = -0.331998556852341
k2 = 0.125557452440262
p1 = -0.000432605884270742
p2 = 0.00110327918082476""")
        temp.flush()
        temp.close()
        fn = temp.name
    return fn 

def gen_d_sub(d_sub_fn, dx, dy, pad_perc=0.1, ndv=-9999):
    nl = dx.shape[0]
    ns = dx.shape[1]
    #Use GDT_Byte or GDT_Int16 to save space?
    dtype = gdal.GDT_Int32
    opt = iolib.gdal_opt
    d_sub_ds = iolib.gtif_drv.Create(d_sub_fn, ns, nl, 3, dtype, opt)
    d_sub_ds.GetRasterBand(1).WriteArray(np.rint(dx.filled(ndv)).astype(np.int32))
    d_sub_ds.GetRasterBand(2).WriteArray(np.rint(dy.filled(ndv)).astype(np.int32))
    d_sub_ds.GetRasterBand(3).WriteArray((~dx.mask).astype(np.int32))
    for n in range(1, d_sub_ds.RasterCount+1):
        band = d_sub_ds.GetRasterBand(n)
        band.SetNoDataValue(float(ndv))
    d_sub_ds = None

    #Now write D_sub_spread.tif - defines spread around D_sub values
    d_sub_ds = iolib.fn_getds(d_sub_fn)
    d_sub_spread_fn = os.path.splitext(d_sub_fn)[0]+'_spread.tif'
    d_sub_spread_ds = gtif_drv.CreateCopy(d_sub_ds)
    dx_spread = np.ma.abs(dx * pad_perc)
    dy_spread = np.ma.abs(dy * pad_perc)
    d_sub_spread_ds.GetRasterBand(1).WriteArray(np.rint(dx_spread.filled(ndv)).astype(np.int32))
    d_sub_spread_ds.GetRasterBand(2).WriteArray(np.rint(dy_spread.filled(ndv)).astype(np.int32))
    d_sub_spread_ds.GetRasterBand(3).WriteArray((~dx_spread.mask).astype(np.int32))
    for n in range(1, d_sub_spread_ds.RasterCount+1):
        band = d_sub_spread_ds.GetRasterBand(n)
        band.SetNoDataValue(float(ndv))
    d_sub_spread_ds = None
    #Copy proj/gt to D_sub and D_sub_spread?

#Return ndarray with h, v, m
def get_vel(fn, fill=True):
    ds = gdal.Open(fn)
    if fill:
        import dem_downsample_fill
        ds = dem_downsample_fill.gdalfill_ds(ds)
    u_b = ds.GetRasterBand(1)
    v_b = ds.GetRasterBand(2)
    u = iolib.b_getma(u_b)
    v = iolib.b_getma(v_b)
    m = np.ma.sqrt(u*u + v*v)
    return u, v, m

#Streamline plot
def streamline_plt(u,v,bg=None):
    dx = 1
    m = np.ma.sqrt(u*u + v*v)
    x,y = np.meshgrid(np.arange(0,u.shape[1],dx),np.arange(0,u.shape[0],dx))
    if bg is None:
        bg = m
    plt.imshow(bg)
    #plt.streamplot(x,y,u,v,color='k', linewidth=5*m/m.max())
    plt.streamplot(x,y,u,v,color='k')
    plt.show()

#Quiver plot
#Currently in disp2v.py

#Inputs can be DEM, DRG, image, etc.
#Only 2 input datsets allowed for this - want to stay modular
def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: %s img1.tif img2.tif" % os.path.basename(sys.argv[0]))

    print('\n%s' % datetime.now())
    print('%s UTC\n' % datetime.utcnow())

    #Should accept these through argparse

    #Integer correlator seeding
    seedmode = 'default' #'sparse_disp', 'velocity

    #Maximum thread count
    maxnthreads = 28

    #Correlator tile timeout
    #With proper seeding, correlation should be very fast
    #timeout = 360 
    timeout = 1200 
    
    #Correlation kernel size
    kernel = (35, 35)
    #kernel = (21, 21)

    #Sub-pixel refinement
    spr = 1

    #Erode disparity map islands smaller than this (area px), set to 0 to turn off
    erode = 1024

    #Set this to smooth the output F.tif with Gaussian filter
    smoothF = False 

    #Set this to mask input to remove vegetation
    #Currently only supports sites where NLCD is available
    #TODO: Need to update to global bare earth
    mask_input = False

    #User-input low-res velocity maps for seeding
    #TODO: Add functions that fetch best available velocities for Ant/GrIS or user-defined low-res velocities
    vx_fn = '' 
    vy_fn = '' 

    #Open input files
    fn1 = sys.argv[1]
    fn2 = sys.argv[2]

    #outdir = '%s__%s_vmap' % (os.path.splitext(os.path.split(fn1)[1])[0], os.path.splitext(os.path.split(fn2)[1])[0])
    outdir = '%s__%s_vmap_%ipx_spm%i' % (os.path.splitext(os.path.split(fn1)[1])[0], os.path.splitext(os.path.split(fn2)[1])[0], kernel[0], spr)
    #Note, issues with boost filename length here, just use vmap prefix
    outprefix = '%s/vmap' % (outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ds1_clip_fn = os.path.join(outdir, os.path.splitext(os.path.basename(fn1))[0]+'_warp.tif')
    ds2_clip_fn = os.path.join(outdir, os.path.splitext(os.path.basename(fn2))[0]+'_warp.tif')

    if not os.path.exists(ds1_clip_fn) or not os.path.exists(ds2_clip_fn):
        #This should write out files to new subdir
        res = 'min'
        #Resample images to something easier to work with
        #res = 2.0
        ds1_clip, ds2_clip = warplib.diskwarp_multi_fn([fn1, fn2], extent='intersection', res=res, r='cubic', outdir=outdir)
        #However, if inputs have identical extent/res/proj, then link to original files
        if not os.path.exists(ds1_clip_fn):
            os.symlink(os.path.abspath(fn1), ds1_clip_fn)
        if not os.path.exists(ds2_clip_fn):
            os.symlink(os.path.abspath(fn2), ds2_clip_fn)

    #Mask support - limit correlation only to rock/ice surfaces, no water/veg
    #This masks input images - guarantee we won't waste time correlating over vegetation
    #TODO: Add support to load arbitrary raster or shp mask
    #TODO: Add support to load global bare earth, not just NLCD
    if mask_input:
        ds1_masked_fn = os.path.splitext(ds1_clip_fn)[0]+'_masked.tif'
        ds2_masked_fn = os.path.splitext(ds2_clip_fn)[0]+'_masked.tif'

        if not os.path.exists(ds1_masked_fn) or not os.path.exists(ds2_masked_fn):
            #Load NLCD or bareground mask
            from demcoreg.dem_mask import get_nlcd, mask_nlcd, get_bareground, mask_bareground

            nlcd_fn = get_nlcd()
            ds1_clip = iolib.fn_getds(ds1_clip_fn)
            #Note: use nearest here to avoid interpolated values
            nlcd_fn_warp = os.path.join(outdir, os.path.split(os.path.splitext(nlcd_fn)[0]+'_warp.tif')[-1])
            if not os.path.exists(nlcd_fn_warp):
                nlcd_ds = warplib.diskwarp_multi_fn([nlcd_fn,], extent=ds1_clip, res=ds1_clip, t_srs=ds1_clip, r='near', outdir=outdir)[0]
            else:
                nlcd_ds = gdal.Open(nlcd_fn_warp)
           
            #Need to check if NLCD or bareground is more appropriate

            nlcd_mask_fn = os.path.join(outdir, 'nlcd_validmask.tif')
            #if not os.path.exists(nlcd_mask_fn):
            validmask = mask_nlcd(nlcd_ds, valid='not_forest', mask_glaciers=False)
            iolib.writeGTiff(validmask, nlcd_mask_fn, nlcd_ds) 

            #bg_fn = get_bareground()
            #bg_ds = warplib.diskwarp_multi_fn([bg_fn,], extent=ds1_clip, res=ds1_clip, t_srs=ds1_clip, r='near', outdir=outdir)[0]
            #validmask = mask_bareground(bg_ds, mask_glaciers=False)

            #Now apply to original images 
            #validmask = validmask.astype(int)
            for fn in (ds1_clip_fn, ds2_clip_fn):
                ds = iolib.fn_getds(fn)
                a = iolib.ds_getma(ds)
                a = np.ma.array(a, mask=~(validmask))
                if a.count() > 0:
                    out_fn = os.path.splitext(fn)[0]+'_masked.tif'
                    iolib.writeGTiff(a,out_fn,ds)
                    a = None
                else:
                    sys.exit("No unmasked pixels over bare earth")

        ds1_clip_fn = ds1_masked_fn
        ds2_clip_fn = ds2_masked_fn

    #Load warped versions on disk
    ds1_clip = iolib.fn_getds(ds1_clip_fn)
    ds2_clip = iolib.fn_getds(ds1_clip_fn)

    #Should have extra kwargs option here
    stereo_opt = get_stereo_opt(maxnthreads=maxnthreads, kernel=kernel, timeout=timeout, erode=erode, spr=spr)
    
    #Stereo arguments
    stereo_args = [ds1_clip_fn, ds2_clip_fn, outprefix]

    #Run stereo_pprc
    if not os.path.exists(outprefix+'-R_sub.tif'):
        run_cmd('stereo_pprc', stereo_opt+stereo_args, msg='0: Preprocessing')
        #Copy proj info to outputs, this should happen automatically now?
        for ext in ('L', 'R', 'L_sub', 'R_sub', 'lMask', 'rMask', 'lMask_sub', 'rMask_sub'):
            geolib.copyproj(ds1_clip_fn, '%s-%s.tif' % (outprefix,ext))

    #Prepare seeding for stereo_corr
    #TODO: these are untested after refactoring
    if not os.path.exists(outprefix+'_D_sub.tif'):
        #Don't need to do anything for default seed-mode 1
        if seedmode == 'sparse_disp':
            #Sparse correlation of full-res images
            stereo_opt.extend(['--corr-seed-mode', '3'])
            sparse_disp_opt = []
            sparse_disp_opt.extend(['--Debug', '--coarse', '512', '--fine', '256', '--no_epipolar_fltr']) 
            sparse_disp_opt.extend(['-P', str(nthreads)])
            sparse_disp_args = [outprefix+'-L.tif', outprefix+'-R.tif', outprefix]
            run_cmd('sparse_disp', sparse_disp_opt+sparse_disp_args, msg='0.5: D_sub generation')
        elif seedmode == 'velocity':
            if os.path.exists(vx_fn) and os.path.exists(vy_fn):
                ds1_res = geolib.get_res(ds1_clip, square=True)[0]

                #Compute L_sub res - use this for output dimensions
                L_sub_fn = outprefix+'-L_sub.tif' 
                L_sub_ds = gdal.Open(L_sub_fn)
                L_sub_x_scale = float(ds1_clip.RasterXSize) / L_sub_ds.RasterXSize
                L_sub_y_scale = float(ds1_clip.RasterYSize) / L_sub_ds.RasterYSize
                L_sub_scale = np.max([L_sub_x_scale, L_sub_y_scale])
                L_sub_res = ds1_res * L_sub_scale

                #Since we are likely upsampling here, use cubicspline
                vx_ds_clip, vy_ds_clip = warplib.memwarp_multi_fn([vx_fn, vy_fn], extent=ds1_clip, res=L_sub_res, r='cubicspline')
                vx = iolib.ds_getma(vx_ds_clip)
                vy = iolib.ds_getma(vy_ds_clip)

                #Determine time interval between inputs
                #Use to scaling of known low-res velocities
                t_factor = timelib.get_t_factor_fn(ds1_clip_fn, ds2_clip_fn, ds=vx_ds_clip)

                if t_factor is not None:
                    #Compute expected offset in scaled pixels 
                    dx = (vx*t_factor)/L_sub_res
                    #Note: Joughin and Rignot's values are positive y up!
                    #ASP is positive y down, so need to multiply these values by -1
                    dy = -(vy*t_factor)/L_sub_res

                    #Should smooth/fill dx and dy

                    #If absolute search window is only 30x30
                    #Don't seed, just use fixed search window 
                    #search_window_area_thresh = 900
                    search_window_area_thresh = 0 
                    search_window = np.array([dx.min(), dy.min(), dx.max(), dy.max()])
                    dx_p = calcperc(dx, perc=(0.5, 99.5))
                    dy_p = calcperc(dy, perc=(0.5, 99.5))
                    search_window = np.array([dx_p[0], dy_p[0], dx_p[1], dy_p[1]])
                    search_window_area = (search_window[2]-search_window[0]) * (search_window[3]-search_window[1])
                    if search_window_area < search_window_area_thresh:
                        stereo_opt.extend(['--corr-seed-mode', '0'])
                        stereo_opt.append('--corr-search')
                        stereo_opt.extend([str(x) for x in search_window])
                        #pad_perc=0.1
                        #stereo_opt.extend(['--corr-sub-seed-percent', str(pad_perc)]
                    #Otherwise, generate a D_sub map from low-res velocity
                    else:
                        stereo_opt.extend(['--corr-seed-mode', '3'])
                        #This is relative to the D_sub scaled disparities
                        d_sub_fn = L_sub_fn.split('-L_sub')[0]+'-D_sub.tif' 
                        gen_d_sub(d_sub_fn, dx, dy)

    #If the above didn't generate a D_sub.tif for seeding, run stereo_corr to generate Low-res D_sub.tif
    if not os.path.exists(outprefix+'-D_sub.tif'):
        newopt = ['--compute-low-res-disparity-only',]
        run_cmd('stereo_corr', newopt+stereo_opt+stereo_args, msg='1.1: Low-res Correlation')
    #Copy projection info to D_sub
    geolib.copyproj(outprefix+'-L_sub.tif', outprefix+'-D_sub.tif')
      
    #Mask D_sub to limit correlation over bare earth surfaces
    #This _should_ be a better approach to masking, but stereo_corr doesn't honor D_sub
    #Now mask input images before stereo_pprc
    #Left this in here for reference, or if this changes in ASP
    if False:
        D_sub_ds = gdal.Open(outprefix+'-D_sub.tif', gdal.GA_Update)

        #Mask support - limit correlation only to rock/ice surfaces, no water/veg
        from demcoreg.dem_mask import get_nlcd, mask_nlcd
        nlcd_fn = get_nlcd()
        nlcd_ds = warplib.diskwarp_multi_fn([nlcd_fn,], extent=D_sub_ds, res=D_sub_ds, t_srs=D_sub_ds, r='near', outdir=outdir)[0]
        #validmask = mask_nlcd(nlcd_ds, valid='rock+ice')
        validmask = mask_nlcd(nlcd_ds, valid='not_forest', mask_glaciers=False)
        nlcd_mask_fn = os.path.join(outdir, 'nlcd_validmask.tif')
        iolib.writeGTiff(validmask, nlcd_mask_fn, nlcd_ds) 

        #Now apply to D_sub (band 3 is valid mask)
        #validmask = validmask.astype(int)
        for b in (1,2,3):
            dsub = iolib.ds_getma(D_sub_ds, b)
            dsub = np.ma.array(dsub, mask=~(validmask))
            D_sub_ds.GetRasterBand(b).WriteArray(dsub.filled())
        D_sub_ds = None

    #OK, finally run stereo_corr full-res integer correlation with appropriate seeding
    if not os.path.exists(outprefix+'-D.tif'):
        run_cmd('stereo_corr', stereo_opt+stereo_args, msg='1: Correlation')
        geolib.copyproj(ds1_clip_fn, outprefix+'-D.tif')

    #Run stereo_rfne
    if not os.path.exists(outprefix+'-RD.tif'):
        run_cmd('stereo_rfne', stereo_opt+stereo_args, msg='2: Refinement')
        geolib.copyproj(ds1_clip_fn, outprefix+'-RD.tif')

    d_fn = make_ln(outdir, outprefix, '-RD.tif')

    #Run stereo_fltr
    if not os.path.exists(outprefix+'-F.tif'):
        run_cmd('stereo_fltr', stereo_opt+stereo_args, msg='3: Filtering')
        geolib.copyproj(ds1_clip_fn, outprefix+'-F.tif')

    d_fn2 = make_ln(outdir, outprefix, '-F.tif')

    if smoothF and not os.path.exists(outprefix+'-F_smooth.tif'):
        print('Smoothing F.tif')
        from pygeotools.lib import filtlib 
        #Fill holes and smooth F
        F_fill_fn = outprefix+'-F_smooth.tif'
        F_ds = gdal.Open(outprefix+'-F.tif', gdal.GA_ReadOnly)
        #import dem_downsample_fill
        #F_fill_ds = dem_downsample_fill.gdalfill_ds(F_fill_ds)
        print('Creating F_smooth.tif')
        F_fill_ds = iolib.gtif_drv.CreateCopy(F_fill_fn, F_ds, 0, options=iolib.gdal_opt)
        F_ds = None
        for n in (1, 2):
            print('Smoothing band %i' % n)
            b = F_fill_ds.GetRasterBand(n)
            b_fill_bma = iolib.b_getma(b)
            #b_fill_bma = iolib.b_getma(dem_downsample_fill.gdalfill(b))
            #Filter extreme values (careful, could lose areas of valid data with fastest v)
            #b_fill_bma = filtlib.perc_fltr(b_fill_bma, perc=(0.01, 99.99))
            #These filters remove extreme values and fill data gaps
            #b_fill_bma = filtlib.median_fltr_skimage(b_fill_bma, radius=7, erode=0)
            #b_fill_bma = filtlib.median_fltr(b_fill_bma, fsize=7, origmask=True)
            #Gaussian filter
            b_fill_bma = filtlib.gauss_fltr_astropy(b_fill_bma, size=9)
            b.WriteArray(b_fill_bma)
        F_fill_ds = None
        d_fn = make_ln(outdir, outprefix, '-F_smooth.tif')

    #Delete intermediate files
    ds1_clip = None
    ds2_clip = None

    print('\n%s' % datetime.now())
    print('%s UTC\n' % datetime.utcnow())
    
    #Generate output velocity products and figure
    cmd = ['disp2v.py', d_fn]
    print("Converting disparities to velocities")
    print(cmd)
    subprocess.call(cmd)

if __name__ == "__main__":
    main()
