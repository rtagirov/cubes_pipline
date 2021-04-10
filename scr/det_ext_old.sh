#!/bin/bash

source aux_func.sh

declare -i Nz Nze i

Nz=150

for i in `seq 24 30`; do

    echo Extending the cube...

    echo

    python ext_cube.py $1 $2 $i

    exit_if_error

    echo $i >> max_tau200.out

    Nze=Nz+i

    sed_exp="1s/.*/$Nze/"

    sed -i $sed_exp dims.inp

    exit_if_error

    echo "max(tau200) calculation"

    ./cpline.exe

    exit_if_error

    echo

done
