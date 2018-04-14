#!/bin/bash
#
# Regenerate gnuplot file. Useful after modifications to template gnuplot file.

set -e

if [ $# != 1 ]; then
    printf "Usage: $0 Exp\n"
    exit 1
fi

length=`cat $1/Exp.config | grep length | awk '{print $2;}'`
length=`bc <<< "$length*20"`
width=`cat $1/Exp.config | grep width | awk '{print $2;}'`
width=`bc <<< "$width*1.2"`
height=`cat $1/Exp.config | grep height | awk '{print $2;}'`
angle=`cat $1/Exp.config | grep alpha | awk '{print $2;}'`
depth=`bc -l <<< "$length*s($angle)/c($angle)"`
depth=`bc <<< "$depth/1+1"`
timeMax=`cat $1/Exp.config | grep timeMax | awk '{print $2;}'`
timeStep=`cat $1/Exp.config | grep timeStep | awk '{print $2;}'`
saveCount=`cat $1/Exp.config | grep saveCount | awk '{print $2;}'`
final=`bc <<< "$timeMax/$timeStep/$saveCount"`
radius=`cat $1/Exp.config | grep radius_large | awk '{print $2;}'`
sed "s/\$Exp/$1/g;s/\$length/$length/g;s/\$width/$width/g;s/\$height/$height/g;s/\$depth/$depth/g;s/\$final/$final/g;s/\$radius/$radius/g" visualize.gnuplot > $1/visualize.gnuplot
chmod +x $1/visualize.gnuplot
