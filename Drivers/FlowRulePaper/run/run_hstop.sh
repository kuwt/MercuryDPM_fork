direc=~/DRIVERS/FlowRulePaper/run
cd $direc
mkdir hstop
ssh node02 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 9 &> hstop/hstop_09.out' &
sleep 20
ssh node02 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 10 &> hstop/hstop_10.out' &
sleep 20
ssh node02 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 11 &> hstop/hstop_11.out' &
sleep 20
ssh node02 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 12 0.5 &> hstop/hstop_12.out' &


