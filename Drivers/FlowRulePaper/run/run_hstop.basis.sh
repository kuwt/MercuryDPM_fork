direc=~/DRIVERS/FlowRulePaper/run
cd $direc
mkdir hstop
ssh node03 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 0.0 0.5 &> hstop/hstop1.out' &
sleep 20
ssh node03 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 0.5 0.5 &> hstop/hstop2.out' &
sleep 20
ssh node03 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 0.666666666666666 0.5 &> hstop/hstop3.out' &
sleep 20
ssh node03 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 0.833333333333333 0.5 &> hstop/hstop4.out' &
sleep 20
ssh node03 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 1.0 0.5 &> hstop/hstop5.out' &
sleep 20
ssh node03 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 1.5 0.5 &> hstop/hstop6.out' &
sleep 20
ssh node06 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 2.0 0.5 &> hstop/hstop7.out' &
sleep 20
ssh node06 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 1.0 0.0 &> hstop/hstop8.out' &
sleep 20
ssh node06 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 1.0 1.0 &> hstop/hstop10.out' &
sleep 20
ssh node06 'cd '$direc'; nohup nice -19 time sc/quick_run hstop 1.0 1e20 &> hstop/hstop11.out' &
