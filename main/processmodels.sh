#!/bin/bash

make clean

export OMP_NUM_THREADS=16
Nparticles=100000
Nsmoother=1000
smootherNparticles=10000
datalabel=cambridge2018GDPD

for model in cambridgeSIthetaTVPlambdaTVP cambridgeSIthetaTVPlambdaCONST cambridgeSIthetaCONSTlambdaTVP cambridgeSIthetaCONSTlambdaCONST
do
    echo $model
    make THIS=$model run FCmode=prod doInflationNoise=1 Nparticles=$Nparticles doSmoother=1 Nsmoother=$Nsmoother DATALABEL=$datalabel smootherNparticles=$smootherNparticles
    make THIS=$model run FCmode=prod doInflationNoise=0 Nparticles=$Nparticles doSmoother=1 Nsmoother=$Nsmoother DATALABEL=$datalabel smootherNparticles=$smootherNparticles
done
