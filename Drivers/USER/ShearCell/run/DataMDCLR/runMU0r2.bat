#!/bin/tcsh
# batch file to prepare and run a shear simulation

# cut a quarter shear-cell out of c3d.old
# and write it to c3d.new
# Input paramters: inner radius, outer radius
#                  maximal height
#                  scaling factor

foreach Zmax (30)

set R1=10
set R2=75

# Rscaled = 0.11 *R/R2
set R1scaled=0.01466667
set R2scaled=0.11000000
set Rcscaled=0.08500000
# the third radius in the first line 
# is the bottom gap radc=0.085/0.11 *R2
# scalefactor=R1scaled/R1

# c3dmake << EOF
# $R1 $R2 57.955
# $Zmax
# 0.001466667
# 2 2 1
# EOF

# prepare simulation startup configuration
\cp c3d_MU0r.restart_30   c3d.ini
\cp drmax_MU0r.restart_30 drmax.ini
echo 14 > c3d_rot.ini

# prepare wall file
echo  1 -1            >  walls.ini
echo  3 -1            >> walls.ini
echo  5 10  0  0.001  >> walls.ini
echo  2 10  0  0.001  >> walls.ini
echo  4 10  0  0.001  >> walls.ini
echo  6 10  0  0.001  >> walls.ini
echo  7  210  $R1scaled  1  0.00  >> walls.ini
echo  8  210  $R2scaled -1  0.01  >> walls.ini
echo  9  233  $Rcscaled  1  0.01  >> walls.ini
echo 10  233  $Rcscaled -1  0.00  >> walls.ini

# prepare parameter file
echo 11 1 -9.81       >  par.ini
echo 250.00 5e-6      >> par.ini
echo 0.01  0.1        >> par.ini
echo 2000 1e2 2e-3 0  >> par.ini
echo 0 0              >> par.ini

# prepare restart at #250
echo 1250 > fstat.ini

# prepare linked cells
echo 0                >  lcell.ini
echo 35 35 20         >> lcell.ini

# prepare the walls for xballs
echo 100 $R1scaled $R1scaled $R1scaled   0  0  0 >  walls.conf
echo 100 $R2scaled $R2scaled $R2scaled   0  0  1 >> walls.conf 
echo 100 $Rcscaled $Rcscaled $Rcscaled   0  0  0 >> walls.conf

echo -1 -10  1            > specnum.ini
echo 0 34518 1           >> specnum.ini
echo 34519 34633 -1 -7   >> specnum.ini
echo 34634 35361 -1 -8   >> specnum.ini
echo 35362 36026 -1 -9   >> specnum.ini
echo 36027 37011 -1 -10  >> specnum.ini

echo 1 > species.ini
echo 2000 0.4  1e2 1.1e2 0  1.2e1 0 0  0.01 0.01 0 0  2e-3 0.5e-3 0 0  0 0 0.05 0 >> species.ini

# perform program
time nice +19 ./mdCLR.icc

# save the output data
\cp c3d.restart c3d_MU0r2.restart
\cp drmax.restart drmax_MU0r2.restart
\cp ene         ene_MU0r2
\cp c3d         c3d_MU0r2

end
