#!/bin/tcsh

# line 1: boundary mode (0|1|11=quarter periodic cylinder)
# line 2: output mode (0/1=without/with stiffness tensor)
# line 3: c3d/data input file-pattern
# line 4: fstat input file-pattern
# line 5: 4-digit number range
# line 6: time-window: tmin, tmax
# line 7: x bins: number, min, max (if min>max use walls from c3d/data)
# line 8: y bins: number, min, max (if min>max use walls from c3d/data)
# line 9: z bins: number, min, max (if min>max use walls from c3d/data)
# line 10: smoothing mode (-1, 0, ...), ishape mode (+2, +/-3, +/-1### +/-3###)
#                                       with ###/100 as stretching factor
# line 11: particle material density

echo 14 > c3d_rot.ini

#gunzip -9 fstat.15*.gz
#gunzip -9 c3d.15*.gz
#gunzip -9 fstat.16*.gz
#gunzip -9 c3d.16*.gz
#gunzip -9 fstat.17*.gz
#gunzip -9 c3d.17*.gz
#gunzip -9 fstat.18*.gz
#gunzip -9 c3d.18*.gz
#gunzip -9 fstat.19*.gz
#gunzip -9 c3d.19*.gz

#gunzip -9 fstat.2*.gz
#gunzip -9 c3d.2*.gz

foreach IS (-3200)     # use 3200, 3400.
foreach NX (0)        # 24.
foreach NY (0)   
foreach NZ (0)   
foreach NS (100)



nice +19 time ~/tav10_new/tav10 << EOF
1
1
./mu0-p3000.data.
./mu0-p3000.fstat.
 0000 -0548 1  #the second number to be the last snapshot 
 0 25000
 $NX 0.2 0.1
 $NY 0.2 0.1
 $NZ 0.2 0.1
 $NS $IS
 2000
EOF

\mv tav10.data tav10new.data

#gzip -9 fstat.16*
#gzip -9 c3d.16*
#gzip -9 fstat.17*
#gzip -9 c3d.17*
#gzip -9 fstat.18*
#gzip -9 c3d.18*
#gzip -9 fstat.19*
#gzip -9 c3d.19*

#gzip -9 fstat.2*
#gzip -9 c3d.2*

end
end
end
end
end

