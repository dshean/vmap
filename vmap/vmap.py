#! /usr/bin/env python

#David Shean
#dshean@gmail.com

#This script uses ASP correlator to produce disparity maps from two inputs
#Input data should be orthorectified/mapped in the same projected coordinate system
#Run disp2v.py to convert to surface velocities

import sys
import os
import argparse
import subprocess
from datetime import datetime, timedelta
from distutils.spawn import find_executable

from osgeo import gdal
import numpy as np

from pygeotools.lib import warplib, geolib, iolib
from pygeotools.lib.malib import calcperc 
from pygeotools.lib.timelib import get_t_factor_fn

#Generate and execute stereo commands
def run_cmd(bin, args, **kw):
    #Note, need to add full executable
    binpath = find_executable(bin)
    if binpath is None:
        msg = ("Unable to find executable %s\n" 
        "Install ASP and ensure it is in your PATH env variable\n" 
        "https://ti.arc.nasa.gov/tech/asr/intelligent-robotics/ngt/stereo/" % bin)
        sys.exit(msg)
    call = [binpath,]
    call.extend(args)
    print(' '.join(call))
    try:
        code = subprocess.call(call, shell=False)
    except OSError as e:
        raise Exception('%s: %s' % (binpath, e))
    if code != 0:
        raise Exception('Stereo step ' + kw['msg'] + ' failed')

def get_stereo_opt(threads=28, kernel=(35,35), nlevels=5, spr=1, timeout=360, erode=0, align='None'):
    stereo_opt = []
    #This is irrelevant
    stereo_opt.extend(['-t', 'pinhole'])
    #Set number of threads/cores to use (will use all CPUs if possible)
    stereo_opt.extend(['--threads', str(threads)])
    #This assumes input images are already mapped 
    stereo_opt.extend(['--alignment-method', align])
    #This will attempt to remove most of the offset between two images, for relative offsets
    #stereo_opt.extend(['--alignment-method', 'Homography'])
    #stereo_opt.append('--ip-debug-images')
    #This should be explored further
    stereo_opt.append('--individually-normalize')
    #Integer correlator kernel size
    stereo_opt.extend(['--corr-kernel', str(kernel[0]), str(kernel[1])])
    stereo_opt.extend(['--corr-max-levels', str(nlevels)])
    if timeout > 0:
        stereo_opt.extend(['--corr-timeout', str(timeout)])
    #Define the search area
    #Useful if you know your orhotorectification is good to, say 100 pixels in any direction
    #stereo_opt.extend(['--corr-search', '-100', '-100', '100', '100'])
    stereo_opt.extend(['--subpixel-mode', str(spr)])
    #If using Semi-global matching (spr 0):
    if spr > 3:
        #Use SGM
        stereo_opt.extend(['--stereo-algorithm', 'asp_sgm'])
        #Use MGM
        #stereo_opt.extend(['--stereo-algorithm', 'asp_mgm'])
        #bro nodes had 128 GB of RAM, 28 threads, ~4.7 GB/thread
        stereo_opt.extend(['--corr-tile-size', '3600'])
        stereo_opt.extend(['--xcorr-threshold', '-1'])
        stereo_opt.extend(['--median-filter-size', '5'])
        stereo_opt.extend(['--texture-smooth-size', '11'])
    else:
        #Sub-pixel kernel size (ASP default is 35)
        #Set to same as integer correlator kernel
        stereo_opt.extend(['--subpixel-kernel', str(kernel[0]), str(kernel[1])])
        stereo_opt.extend(['--stereo-algorithm','asp_bm'])
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
    d_sub_spread_ds = iolib.gtif_drv.CreateCopy(d_sub_spread_fn, d_sub_ds, 0)
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

def getparser():
    parser = argparse.ArgumentParser(description="Generate velocity map via feature-tracking")
    parser.add_argument('-outdir', default=None, help='Output directory')
    parser.add_argument('-threads', type=int, default=iolib.cpu_count(), help='Number of threads to use(default: %(default)s)')
    parser.add_argument('-tr', default='min', help='Output resolution (default: %(default)s)')
    #Set correlator kernel size
    parser.add_argument('-kernel', type=int, default=35, help='Correlator kernel size. Smaller kernels offer more detail but are prone to more noise. Odd integers required (~9-51 px recommended). (default: %(default)s)')
    align_choices = ['AffineEpipolar', 'Homography', 'Epipolar', 'None']
    parser.add_argument('-align', default='None', choices=align_choices, help='Alignment method to warp second image to match first image, if not already orthorectified. Provides flexibility for L1B inputs')
    #Integer correlator seeding
    #D_sub is low-resolution correlation (default), which works well for most situations
    #sparse_disp will use sparse seeding from full-res chips, useful for ice sheets with limited low-frequency texture
    #existing_velocity will accept existing vx and vy rasters.  Useful for limiting search range and limiting blunders.  Measures products are useful for ice sheets.
    seedmode_choices = ['D_sub', 'sparse_disp', 'existing_velocity']
    parser.add_argument('-seedmode', type=str, choices=seedmode_choices, default='D_sub', help='Seeding option (default: %(default)s)')
    parser.add_argument('-vx_fn', type=str, default=None, help='Seed E-W velocity map filename')
    parser.add_argument('-vy_fn', type=str, default=None, help='Seed N-S velocity map filename')
    
    #Sub-pixel refinement
    #0) None, 1) Parabolic, 2) Bayes, 3) AffineAdaptive
    #See ASP doc or Shean et al, ISPRS, (2016)
    #1 is fast but lower quality
    #2 is slow but highest quality, 
    #3 is a good compromise for speed and quality
    refinement_choices = list(range(12))
    parser.add_argument('-refinement', type=int, default=1, help='Sub-pixel refinement type (see ASP doc): 0) None, 1) Parabolic, 2) Bayes, 3) AffineAdaptive 4) LK, 5) Bayes w/gamma, 6) SGM Linear, 7) SGM Poly4, 8) SGM Cos, 9) SGM Parabola, 10) SGM None, 11) SGM Blend (default: %(default)s)')
    #Numer of gaussian pyramids to use
    #Can look at texture in GDAL overviews to make a decision
    #If you can see plenty of texture at 1/32 resolution, go with 5 
    #For featureless areas, limiting to 2 can help, or even 0
    parser.add_argument('-pyramid-levels', type=int, default=5, help='Number of pyramid levels for correlation (default: %(default)s)')
    #This helps get rid of bogus "islands" in the disparity maps
    parser.add_argument('-erode', type=int, default=1024, help='Erode isolated blobs smaller than this many pixels. Set to 0 to disable (default: %(default)s)')
    parser.add_argument('-filter', action='store_true', help='Filter the output F.tif, smoothing with Gaussian filter')
    #This masks input images to improve performance.  Useful for forested areas.
    parser.add_argument('-mask_input', action='store_true', help='Mask any vegetation/water in input images. Requires demcoreg')
    parser.add_argument('-remove_offsets', action='store_true', help='Remove median horizontal and vertical offsets over stable control surfaces')
    parser.add_argument('-dt', type=str, choices=['yr','day','none'], default='yr', help='Time increment (default: %(default)s)')

    #Inputs can be images, DEMs, shaded relief maps
    #Personal experience suggests multi-directional hillshades with identical illumination work well
    #Only 2 input datsets allowed for this - want to stay modular
    parser.add_argument('fn1', type=str, help='Raster filename 1')
    parser.add_argument('fn2', type=str, help='Raster filename 2')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()
    if args.seedmode == 'existing_velocity':
        if args.vx_fn is None or args.vy_fn is None:
            parser.error('"-seedmode existing_velocity" requires "-vx_fn" and "-vy_fn"')

    print('\n%s' % datetime.now())
    print('%s UTC\n' % datetime.utcnow())

    align = args.align
    seedmode = args.seedmode
    spr = args.refinement
    erode = args.erode
    #Correlator tile timeout
    #With proper seeding, correlation should be very fast
    #timeout = 360 
    timeout = 1200 
    threads = args.threads

    kernel = (args.kernel, args.kernel)
    #SGM correlator
    if spr > 3:
        #kernel = (7,7)
        kernel = (11,11)
        erode = 0

    #Smooth the output F.tif 
    smoothF = args.filter 

    res = args.tr
    #Resample input to something easier to work with
    #res = 4.0

    #Open input files
    fn1 = args.fn1
    fn2 = args.fn2 

    if not iolib.fn_check(fn1) or not iolib.fn_check(fn2):
        sys.exit("Unable to locate input files")

    if args.outdir is not None:
        outdir = args.outdir
    else:
        outdir = '%s__%s_vmap_%sm_%ipx_spm%i' % (os.path.splitext(os.path.split(fn1)[1])[0], \
                os.path.splitext(os.path.split(fn2)[1])[0], res, kernel[0], spr)

    #Note, can encounter filename length issues in boost, just use vmap prefix
    outprefix = '%s/vmap' % (outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #Check to see if inputs have geolocation and projection information
    ds1 = iolib.fn_getds(fn1)
    ds2 = iolib.fn_getds(fn2)

    if geolib.srs_check(ds1) and geolib.srs_check(ds2):
        ds1_clip_fn = os.path.join(outdir, os.path.splitext(os.path.basename(fn1))[0]+'_warp.tif')
        ds2_clip_fn = os.path.join(outdir, os.path.splitext(os.path.basename(fn2))[0]+'_warp.tif')

        if not os.path.exists(ds1_clip_fn) or not os.path.exists(ds2_clip_fn):
            #This should write out files to new subdir
            ds1_clip, ds2_clip = warplib.diskwarp_multi_fn([fn1, fn2], extent='intersection', res=res, r='average', outdir=outdir)
            ds1_clip = None
            ds2_clip = None
            #However, if inputs have identical extent/res/proj, then link to original files
            if not os.path.exists(ds1_clip_fn):
                os.symlink(os.path.abspath(fn1), ds1_clip_fn)
            if not os.path.exists(ds2_clip_fn):
                os.symlink(os.path.abspath(fn2), ds2_clip_fn)
            align = 'None'

        #Mask support - limit correlation only to rock/ice surfaces, no water/veg
        #This masks input images - guarantee we won't waste time correlating over vegetation
        #TODO: Add support to load arbitrary raster or shp mask
        if args.mask_input:
            ds1_masked_fn = os.path.splitext(ds1_clip_fn)[0]+'_masked.tif'
            ds2_masked_fn = os.path.splitext(ds2_clip_fn)[0]+'_masked.tif'

            if not os.path.exists(ds1_masked_fn) or not os.path.exists(ds2_masked_fn):
                #Load NLCD or bareground mask
                from demcoreg.dem_mask import get_lulc_mask

                ds1_clip = iolib.fn_getds(ds1_clip_fn)
                lulc_mask_fn = os.path.join(outdir, 'lulc_mask.tif')
                #if not os.path.exists(nlcd_mask_fn):
                lulc_mask = get_lulc_mask(ds1_clip, mask_glaciers=False, filter='not_forest')
                iolib.writeGTiff(lulc_mask, lulc_mask_fn, ds1_clip) 
                ds1_clip = None

                #Now apply to original images 
                #This could be problematic for huge inputs, see apply_mask.py
                #lulc_mask = lulc_mask.astype(int)
                for fn in (ds1_clip_fn, ds2_clip_fn):
                    ds = iolib.fn_getds(fn)
                    a = iolib.ds_getma(ds)
                    a = np.ma.array(a, mask=~(lulc_mask))
                    if a.count() > 0:
                        out_fn = os.path.splitext(fn)[0]+'_masked.tif'
                        iolib.writeGTiff(a,out_fn,ds)
                        a = None
                    else:
                        sys.exit("No unmasked pixels over bare earth")
            ds1_clip_fn = ds1_masked_fn
            ds2_clip_fn = ds2_masked_fn
    else:
        ds1_clip_fn = fn1
        ds2_clip_fn = fn2
        #Now let user specify alignment methods as option - don't hardcode
        #align = 'Homography'
        #align = 'AffineEpipolar'
    ds1 = None
    ds2 = None

    #Should have extra kwargs option here
    stereo_opt = get_stereo_opt(threads=threads, kernel=kernel, timeout=timeout, \
            erode=erode, spr=spr, align=align)
    
    #Stereo arguments
    #Latest version of ASP should accept tif without camera models
    #stereo_args = [ds1_clip_fn, ds2_clip_fn, outprefix]
    #Nope - still need to provide dummy camera models, and they must be unique files
    #Use the dummy.tsai file bundled in the vmap repo
    dummy_tsai = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'dummy.tsai')
    dummy_tsai2 = os.path.splitext(dummy_tsai)[0]+'2.tsai'
    if not os.path.exists(dummy_tsai2):
        dummy_tsai2 = os.symlink(dummy_tsai, os.path.splitext(dummy_tsai)[0]+'2.tsai')
    stereo_args = [ds1_clip_fn, ds2_clip_fn, dummy_tsai, dummy_tsai2, outprefix]

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
            sparse_disp_opt.extend(['-P', str(threads)])
            sparse_disp_args = [outprefix+'-L.tif', outprefix+'-R.tif', outprefix]
            run_cmd('sparse_disp', sparse_disp_opt+sparse_disp_args, msg='0.5: D_sub generation')
        elif seedmode == 'existing_velocity':
            #User-input low-res velocity maps for seeding
            #TODO: Add functions that fetch best available velocities for Ant/GrIS or user-defined low-res velocities
            #Automatically query GoLive velocities here
            vx_fn = args.vx_fn 
            vy_fn = args.vy_fn 
            #Check for existence

            #HMA seeding
            #vdir = '/nobackup/deshean/rpcdem/hma/velocity_jpl_amaury_2013-2015'
            #vx_fn = os.path.join(vdir, 'PKH_WRS2_B8_2013_2015_snr5_n1_r170_res12.x_vel.TIF')
            #vy_fn = os.path.join(vdir, 'PKH_WRS2_B8_2013_2015_snr5_n1_r170_res12.y_vel.TIF')

            if os.path.exists(vx_fn) and os.path.exists(vy_fn):
                ds1_clip = iolib.fn_getds(ds1_clip_fn)
                ds1_res = geolib.get_res(ds1_clip, square=True)[0]

                #Compute L_sub res - use this for output dimensions
                L_sub_fn = outprefix+'-L_sub.tif' 
                L_sub_ds = gdal.Open(L_sub_fn)
                L_sub_x_scale = float(ds1_clip.RasterXSize) / L_sub_ds.RasterXSize
                L_sub_y_scale = float(ds1_clip.RasterYSize) / L_sub_ds.RasterYSize
                L_sub_scale = np.max([L_sub_x_scale, L_sub_y_scale])
                L_sub_res = ds1_res * L_sub_scale

                #Since we are likely upsampling here, use cubicspline
                vx_ds_clip, vy_ds_clip = warplib.memwarp_multi_fn([vx_fn, vy_fn], extent=ds1_clip, \
                        t_srs=ds1_clip, res=L_sub_res, r='cubicspline')

                ds1_clip = None

                #Get vx and vy arrays
                vx = iolib.ds_getma(vx_ds_clip)
                vy = iolib.ds_getma(vy_ds_clip)

                #Determine time interval between inputs
                #Use to scaling of known low-res velocities
                t_factor = get_t_factor_fn(ds1_clip_fn, ds2_clip_fn, ds=vx_ds_clip)

                if t_factor is not None:
                    #Compute expected offset in scaled pixels 
                    dx = (vx*t_factor)/L_sub_res
                    dy = (vy*t_factor)/L_sub_res
                    #Note: Joughin and Rignot's values are positive y up!
                    #ASP is positive y down, so need to multiply these values by -1
                    #dy = -(vy*t_factor)/L_sub_res

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
    #This _should_ be a better approach than masking input images, but stereo_corr doesn't honor D_sub
    #Still need to mask input images before stereo_pprc
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
    if spr > 0:
        if not os.path.exists(outprefix+'-RD.tif'):
            run_cmd('stereo_rfne', stereo_opt+stereo_args, msg='2: Refinement')
            geolib.copyproj(ds1_clip_fn, outprefix+'-RD.tif')
        d_fn = make_ln(outdir, outprefix, '-RD.tif')
    else:
        ln_fn = outprefix+'-RD.tif'
        if os.path.lexists(ln_fn):
            os.remove(ln_fn)
        os.symlink(os.path.split(outprefix)[1]+'-D.tif', ln_fn)

    #Run stereo_fltr
    if not os.path.exists(outprefix+'-F.tif'):
        run_cmd('stereo_fltr', stereo_opt+stereo_args, msg='3: Filtering')
        geolib.copyproj(ds1_clip_fn, outprefix+'-F.tif')

    d_fn = make_ln(outdir, outprefix, '-F.tif')

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

    print('\n%s' % datetime.now())
    print('%s UTC\n' % datetime.utcnow())

    #If time interval is specified, convert pixel displacements to rates
    if args.dt != 'none':
        #Check if vm.tif already exists
        #Should probably just overwrite by default
        #if os.path.exists(os.path.splitext(d_fn)[0]+'_vm.tif'):
        #    print("\nFound existing velocity magnitude map!\n"
        #else:
        #Generate output velocity products and figure
        #Requires that vmap repo is in PATH
        cmd = ['disp2v.py', d_fn]
        #Note: this will attempt to automatically determine control surfaces
        #disp2v.py will accept arbitrary mask, could pass through here
        if args.remove_offsets:
            cmd.append('-remove_offsets')
        cmd.extend(['-dt', args.dt])
        print("Converting disparities to velocities")
        print(cmd)
        subprocess.call(cmd)

if __name__ == "__main__":
    main()
