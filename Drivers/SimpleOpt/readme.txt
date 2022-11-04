These driver files demonstrate a simple simulation-guided optimization with MercuryDPM simulation and 
numpy/scipy optimizer accessed via SimpleOpt tool

Simulation-guided optimization iteratively runs the driver file and varies parameters of optimization
to ensure extremum of the chosen functional.

File opt.cpp is an example of such driver. Optimization parameters are passed as command line arguments. 
Note that if you need to use command line arguments in solve() routine, you should remove optimization 
parameters from the arguments passed to solve(argc, argv).

The functional to be extremized is defined in main funciton of opt.cpp, an every optimization step it is 
computed and and saved to the file functional.txt in the executable home directory.

The file opt_main.cpp calls the optimization tool with extra command line parameters: MercuryDPM source dir
and mercuryDPM build dir - this allows SimpleOpt tool to access the executable and the functional.

The tool uses one of the optimization algorithms from scipy.optimize library to perform a simulation-guided
optimization.    
