#!/bin/bash

# make clean

export OMP_NUM_THREADS=16

Ngrid=250

datalabel=cambridge2018CPI

for model in cambridgeSIthetaTVPlambdaTVP cambridgeSIthetaTVPlambdaCONST cambridgeSIthetaCONSTlambdaTVP cambridgeSIthetaCONSTlambdaCONST
do
    echo $model
    for Nparticles in 100000 10000 1000
    do
	echo $Nparticles
	make THIS=$model runstderr FCmode=prod Nparticles=$Nparticles Ngrid=$Ngrid doInflationNoise=1 DATALABEL=$datalabel
	make THIS=$model runstderr FCmode=prod Nparticles=$Nparticles Ngrid=$Ngrid doInflationNoise=0 DATALABEL=$datalabel
    done
done
