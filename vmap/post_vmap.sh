#! /bin/bash

#Generate standardized stacks/products from a directory of vmap output

mkdir rename 

for j in vm vx vy
do
    echo $j
    mkdir rename/${j}
    #This isn't working, fn doesn't get assigned
    #parallel "fn=$(~/src/vmap/vmap/wv_cdate.py {}); ln -s ../{} rename/${j}/${fn}_${j}.tif" ::: 2*/2*${j}.tif
    for i in 2*/2*${j}.tif
    do 
        ln -sv ../../$i rename/${j}/$(~/src/vmap/vmap/wv_cdate.py $i)_${j}.tif
    done
done

cd rename

#Create median, stddev and count products
parallel "dem_mosaic --{2} -o {1} {1}/*{1}.tif" ::: vm vx vy ::: median stddev count

#Clip to RGI poylgons
parallel "clip_raster_by_shp.py -extent raster {}-tile-0-median.tif rgi" ::: vm vx vy

#Clip to rock surfaces
parallel "dem_mask.py {}-tile-0-median.tif" ::: vm vx vy

#Compute stats for rock surfaces
for i in *ref.tif ; do robust_stats.py $i > ${i%.*}_stats.txt ; done

#tar -hczvf khumbu_vel.tar.gz rename
