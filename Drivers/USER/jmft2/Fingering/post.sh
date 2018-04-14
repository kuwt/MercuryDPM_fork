#!/bin/bash
#
# Generate visualization and summary page.

set +e

if [ $# != 1 ]; then
    printf "Usage: $0 Exp\n"
    exit 1
fi

./HHT.py $1
./profile.py $1
./segregation.py $1
mkdir ~/public_html/Fingering-results/$1
#chmod 755 ~/public_html/Fingering-results/$1
cat $1/Exp.config > ~/public_html/Fingering-results/$1/Exp.config
$1/visualize.gnuplot
for name in side bottom top height speed vfield; do avconv -nostdin -y -loglevel panic -i $1/$name%03d.png -c:v libx264 -pix_fmt yuv444p -q 0 ~/public_html/Fingering-results/$1/$name.mp4; done
read N[{1..3}] <<< `head -n 3 $1/Exp.txt | awk '{print $1;}'`
sed "s/\$Exp/$1/g;s/\$NSmall/${N[1]}/g;s/\$NLarge/${N[2]}/g;s/\$NBase/${N[3]}/g" Summary.html > ~/public_html/Fingering-results/$1/Summary.html
#chmod 644 ~/public_html/Fingering-results/$1/*

echo "http://damtp.cam.ac.uk/user/$USER/Fingering-results/$1/Summary.html" | mail -s "$1 finished!" $USER
