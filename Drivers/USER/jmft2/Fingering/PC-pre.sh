#!/bin/bash
#
# Generate config file and gnuplot file for a simulation.

set -e

if [ $# == 1 ]; then
    template="test.config"
elif [ $# == 2 ]; then
    if [[ $2 == *.config && -f $2 ]]; then
        template=$2
    elif [[ -f $2".config" ]]; then
        template=$2".config"
    elif [[ -f $2/Exp.config ]]; then
        template=$2/Exp.config
    else
        printf "Cannot find the config file $2!\n"
        exit 1
    fi
else
    printf "Usage: $0 Exp [template]\n"
    exit 1
fi

if [ -d $1 ]; then
    printf "$1 already exists!\n"
    exit 1
fi

mkdir $1
cp $template $1/Exp.config
# vim $1/Exp.config
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
sed "s/\$Exp/$1/g;s/\$length/$length/g;s/\$width/$width/g;s/\$height/$height/g;s/\$depth/$depth/g;s/\$final/$final/g;s/\$radius/$radius/g" PC-visualize.gnuplot > $1/visualize.gnuplot
chmod +x $1/visualize.gnuplot
# vim $1/visualize.gnuplot

read -p "Start running $1? [y/N]: " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    ./PC-run.sh $1 &
fi
