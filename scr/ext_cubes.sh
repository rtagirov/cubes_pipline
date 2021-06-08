#!/bin/bash

scr=/mnt/SSD/sim/cubes_pipline/scr
inp=/mnt/SSD/sim/cubes_pipline/inp
exe=/mnt/SSD/sim/cubes_pipline/src/cpline.exe

files=$(ls cube_small.*.nc)

for f in $files; do

    IFS='.' read -ra arr <<< $f

    num=${arr[1]}

    rm -rf $num

    mkdir $num

    cd $num

#    rm -f ${num}.nc abund.inp dims.inp

    ln -s ../cube_small.${num}.nc ./${num}.nc

    ln -s $inp/abund_cubes_veronika.dat ./abund.inp
    ln -s $inp/dims.inp                 ./

#    rm -f ext_cube.py det_ext.sh aux_func.sh

    ln -s $scr/det_ext.sh  ./
    ln -s $scr/aux_func.sh ./
    ln -s $scr/ext_cube.py ./

    cp $exe ./

    cp $inp/INPUT_ext       ./INPUT
    cp $inp/ross_table.dat  ./
    cp $inp/kappa_table.dat ./

    echo $num > snapshot.inp

    ./det_ext.sh $num nc &

    cd ..

done
