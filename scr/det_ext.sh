#!/bin/bash

source aux_func.sh

declare -i Nz Nze i

Nz=150

i=20

rm -f extension.log

while true; do

    echo "\n"Extending the cube..."\n" >> extension.log

    python ext_cube.py $1 $2 $i &

    wait

    exit_if_error

    echo $i >> max_tau200.out

    Nze=Nz+i

    sed_exp="1s/.*/$Nze/"

    sed -i $sed_exp dims.inp

    exit_if_error

    echo "max(tau200) calculation" >> extension.log

    ./cpline.exe >> extension.log 2>&1

    wait

    exit_if_error

    maxtau=$(tail -1 max_tau200.out)

    if (($(echo "$maxtau < 0.2" | bc -l))); then

        mv *.bin ../

        break

    fi

    i=i+1

done
