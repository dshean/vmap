#! /bin/bash

infile=$1

#Percent size of original, default is 5%
rs=5
#Resize to 5% unless we're dealing with subsampled D map
if [[ "$infile" == *D_sub* ]] ; then
    rs=100
fi

#Default ASP disparity ndv is set to the minimum value

#Can grab this from gdalinfo -stats, but slow.  
#Might be best to load into masked array, using the appropriate LMask and RMask files
#Then compute stats to find minimum value (but what if disparity range extends beyond specified min)
ndv=0

#Run parallel gdal_translate
#Use default nearest neighbor for fast downsampling
echo "Resampling b1 and b2 at ${rs}% of original size"
parallel --bar "if [ ! -e ${infile%.*}_b{}_${rs}.tif ] ; then gdal_translate -b {} -a_nodata $ndv -outsize ${rs}% ${rs}% $infile ${infile%.*}_b{}_${rs}.tif; fi" ::: 1 2 

#View the result
imview.py -link ${infile%.*}_b1_${rs}.tif ${infile%.*}_b2_${rs}.tif 
