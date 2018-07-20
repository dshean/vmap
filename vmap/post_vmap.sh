#! /bin/bash

#Generate standardized stacks/products from a directory of vmap output

if [ ! -d rename ] ; then
    mkdir rename 
    for j in vm vx vy
    do
        echo $j
        mkdir rename/${j}
        for i in 2*/2*${j}.tif
        do 
            dst_fn=rename/${j}/$(~/src/vmap/vmap/wv_cdate.py $i)_${j}.tif
            #For some sites with multiple DEMs acquired on same pass, can have conflicts for timestamps
            #Hack to append _1 until we have unique filenames
            while [ -e $dst_fn ]
            do
                dst_fn=$(echo $dst_fn | sed "s/_${j}/_1_${j}/")
            done
            ln -sv ../../$i $dst_fn
        done
    done
fi

cd rename

#Create median, stddev and count products
parallel --verbose "dem_mosaic --{2} -o {1} {1}/*{1}.tif" ::: vm vx vy ::: median stddev count
#The median run can use a lot of memory (50 GB for 700 inputs at Rainier)
#Can lead to memory errors on Pleiades nodes, so split processing
#parallel --verbose -j 2 "dem_mosaic --threads 10 --median -o {1} {1}/*{1}.tif" ::: vm vx vy 
#parallel --verbose "dem_mosaic --{2} -o {1} {1}/*{1}.tif" ::: vm vx vy ::: stddev count

#Clip to RGI poylgons
parallel "clip_raster_by_shp.py -extent raster {}-tile-0-median.tif rgi" ::: vm vx vy

#Clip to rock surfaces
parallel "dem_mask.py {}-tile-0-median.tif" ::: vm vx vy

#Compute stats for rock surfaces
for i in *ref.tif ; do robust_stats.py $i > ${i%.*}_stats.txt ; done

make_stack.py --med -outdir stack vm/*vm.tif

#tar -hczvf khumbu_vel.tar.gz rename
