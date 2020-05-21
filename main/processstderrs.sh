#!/bin/bash

make clean

Ngrid=250

for model in cambridgeSIthetaTVPlambdaTVP cambridgeSIthetaTVPlambdaCONST cambridgeSIthetaCONSTlambdaTVP cambridgeSIthetaCONSTlambdaCONST
do
    echo $model
    for Nparticles in 10000 #100000 10000 1000
    do
	echo $Nparticles
	make THIS=$model runstderr FCmode=prod Nparticles=$Nparticles Ngrid=$Ngrid doInflationNoise=1
	make THIS=$model runstderr FCmode=prod Nparticles=$Nparticles Ngrid=$Ngrid doInflationNoise=0
    done
done
