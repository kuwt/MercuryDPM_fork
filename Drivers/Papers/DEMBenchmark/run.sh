#!/bin/bash
# script to run all simulations for the OS paper
nohup ./OSSilo -useSmallOrifice -useM1 > SiloSmallM1.out &
nohup ./OSSilo -useM1 > SiloLargeM1.out &
nohup ./OSSilo -useSmallOrifice > SiloSmallM2.out &
nohup ./OSSilo > SiloLargeM2.out &
nohup ./OSPenetration -size 25K > Penetration25K.out &
nohup ./OSPenetration -size 50K > Penetration50K.out &
nohup ./OSPenetration -size 100K > Penetration100K.out &
nohup ./OSDrum > Drum.out &
