../../snapshot.exe ../../../FlowRulePaper/run/full_runs/restart_flowrule/H10A22L0M0.5B0.5
../../snapshot.exe ../../../FlowRulePaper/run/full_runs/restart_flowrule/H10A26L2M0.5B0.5
../../snapshot.exe ../static2d/static2d

convert 000000.H10A22L0M0.5B0.5.pdf L0.gif
convert -rotate 22 L0.eps L0.eps
convert 000000.H10A26L2M0.5B0.5.pdf L2.gif
convert -rotate 26 L2.gif L2.gif
convert 000000.static2d.pdf 2D.gif

for file in $(ls ??.gif)
do
	convert -trim $file $file
done
