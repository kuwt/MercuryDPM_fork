#!/usr/bin/gnuplot

set term pngcairo
set xrange [0:$length]
set yrange [-$width:$width]
set zrange [-$depth:$height]
set cbrange [0:1]
set palette model RGB defined (0 'red', 1 'blue')
unset colorbox
do for [i=0:$final] {
  set output sprintf('$Exp/side%03d.png',i)
  splot '$Exp/Exp.data.'.i u 1:3:2:($14+1)/5:14 every ::2 w p pal pt 7 ps var tit sprintf('%i',i)
}

set view equal xyz
set view 180,0,2
do for [i=0:$final] {
  set output sprintf('$Exp/bottom%03d.png',i)
  splot '$Exp/Exp.data.'.i u 1:3:2:($14+1)/5:14 every ::2 w p pal pt 7 ps var tit sprintf('%i',i)
}

set view 0,0,2
do for [i=0:$final] {
  set output sprintf('$Exp/top%03d.png',i)
  splot '$Exp/Exp.data.'.i u 1:3:2:($14+1)/5:14 every ::2 w p pal pt 7 ps var tit sprintf('%i',i)
}

unset view
set pm3d map interp 0,0
set colorbox
set cbrange [0:0.3]
do for [i=0:$final] {
  set output sprintf('$Exp/height%03d.png',i)
  splot '$Exp/height.data.'.i matrix u ($1*$radius*4):($2*$radius*4-$width):3 tit sprintf('%i',i)
  set output sprintf('$Exp/speed%03d.png',i)
  splot '$Exp/speed.data.'.i matrix u ($1*$radius*4):($2*$radius*4-$width):3 tit sprintf('%i',i)
}

set autoscale
unset key
unset colorbox
set cbrange [0:1]
set output '~/public_html/Fingering-results/$Exp/front.png'
plot '$Exp/Exp.data.$final' u 3:1:($14+1)/5:14 every ::2 w p pal pt 7 ps var, '$Exp/front.data' u 1:2 w l lc rgb "black" lw 3

set output '~/public_html/Fingering-results/$Exp/ene.png'
plot '$Exp/Exp.ene' u 1:4 w l

set output '~/public_html/Fingering-results/$Exp/txt.png'
plot '$Exp/Exp.txt' u 1:3 w l

set output '~/public_html/Fingering-results/$Exp/EMD.png'
N=`head -n 1 $Exp/EMD.data | wc -w`
plot [-$width:$width] for [i=2:N] '$Exp/EMD.data' u 1:i w l

set output '~/public_html/Fingering-results/$Exp/Hilbert.png'
N=`head -n 1 $Exp/Hilbert.data | wc -w`
plot [-$width:$width] for [i=2:N] '$Exp/Hilbert.data' u 1:i w l

set size ratio -1
set colorbox
set cbrange [0:*]
do for [i=0:$final] {
  set output sprintf('$Exp/vfield%03d.png',i)
  set multiplot layout 2,1 tit sprintf('%i',i)
  stats '$Exp/segregation.data.'.i nooutput
  plot [floor(STATS_max_x-0.5):floor(STATS_max_x+1.5)][-0.02:0.2] '$Exp/segregation.data.'.i u 1:2:3:4:5 w vec lc pal
  plot [floor(STATS_max_x-0.5):floor(STATS_max_x+1.5)][-0.02:0.2] '$Exp/segregation.data.'.i u 1:2:6:7:8 w vec lc pal
  unset multiplot
}
