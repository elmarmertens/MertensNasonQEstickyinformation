#!/bin/bash

# make cleanall
export OMP_NUM_THREADS=32

Ngrid=640
Nparticles=10000
doInflationNoise=1
T=200

make -f mdd.makefile run FCmode=prod Nparticles=$Nparticles MDDgrid=10 Ngrid=$Ngrid doInflationNoise=$doInflationNoise THIS=cambridgesimMDD T=$T
