# How to calibrate material parameters

This guide shows you how to calibrate bulk material parameters in MercuryDPM using GrainLearning. I will first explain the idea, followed by step-by-step instructions on how to run the calibration scripts.

## The idea

Assume you want to simulate a particulate process using MercuryDPM for a specific material A. You have calibration data, such as:
 - the static angle of repose from a heap test, 
 - the dynamic angle of repost from a rotating drum, and/or
 - a set of incipient shear stresses and corresponding normal stresses from a shear tester

To simulate this in MercuryDPM, you need to (a) decide on a particle and contact model and (b) find the parameters of these particle and contact models such that the simulated material behaves like the real material. This is not trivial, since the above calibration tests result are bulk measures, i.e. have no 1:1-correspondence to the contact model parameters. 

To calibrate, we will use sensible assumptions and direct measurements to whittle down the number of free particle and contact model parameters, then use GrainLearning to caliibrate the free parameters from teh calibration data.

## A test example

Before we run a realistic example, consider the following test. Say we have a particle and contact model with one input parameter, `param0`, and a single calibration test that returns a scalar output value `0.5`. This calibration test is implemented in the executable `CalibrationTest -fit identity1`; it simply returns the input parameter as the output parameter. To test whether GrainLearning can find the value of `param0` that produces the output value `0.5`, do the following:

- Make a file named `calibrationTest_identity1.txt` and place it in the Drivers/Calibration/py folder. For this example, the file already exists, please take a look. The file name contains the name of the output directory (`Test_identity1`), prepended `calibrate` and appended `.txt`. It contains three sections, each started by a header line:
   - Section 1: arguments that should be passed to the executable via command line interface (in this case, `-fit identity1`)
   - Section 2: arguments that should be chosen by GrainLearning (in this case, `-param0 val`, with *val* a value chosen in the specified range of 0 to 1 )
   - Section 3: characterisation data (in this case the value 0.5) and the associated executable (in this case, `CalibrationTest`)
- GrainLearning is executed by opening a terminal and running the script `./calibrate.py --material Test_identity1` in the `py/` directory. It will iterate three times, with 10 iterations per iteration. 
- The simulation output is stored in the output directory `Test_identity11`. 
  - The characterisation data will be stored in subdirectory `Exp/`
  - The 0-th iteration data will be stored in subdirectory `Sim0/`, the 1-st iteration in `Sim1/`, etc
- The console output should be as follows:
    ````
    Running in serial
    Material: Test_identity1
    Parameters to be identified: param0
    Parameter ranges: [0.0, 1.0]
    Executables: CalibrationTest
    Data values: [0.5]
    Weights: [1.0]
    Command line parameters: -speciesType Test -fit identity1 -param0 %s
    Build directory for executables: /Users/weinhartt/Code/Calibration/cmake-build-debug/Drivers/Calibration
    Covariance: 1
    Goal for effective sample size: 0.500000
    Number of Gaussian components: 2
    Number of samples: 10
    Prior weight: 0.500000
    
    Creating smcTable_Test_identity1_0.txt.
    Output will be written to Test_identity1/Sim0/.
    Preparing simulations. After simulations have finished, rerun this script to continue.
    Running make in the build directory /Users/weinhartt/Code/Calibration/cmake-build-debug/Drivers/Calibration
    Runs started. Rerun this script for analysis.
    Starting analysis of iteration 0
    Skipping DEM simulations, read in pre-existing simulation data now
    The optimal parameters for Test_identity1 in iteration 0 are:
    Best params for Test_identity1 in iteration 0 (error 0.000000):
    param0: 0.5
    data: CalibrationTest: 0.5 (exp: 0.5)
    Means and Covariances:
    param0 = mean 0.501352 with covariance 0.222814
    
    Creating smcTable_Test_identity1_1.txt.
    Output will be written to Test_identity1/Sim1/.
    Preparing simulations. After simulations have finished, rerun this script to continue.
    Running make in the build directory /Users/weinhartt/Code/Calibration/cmake-build-debug/Drivers/Calibration
    Runs started. Rerun this script for analysis.
    Starting analysis of iteration 1
    Skipping DEM simulations, read in pre-existing simulation data now
    The optimal parameters for Test_identity1 in iteration 1 are:
    Best params for Test_identity1 in iteration 1 (error 0.031246):
    param0: 0.515623
    data: CalibrationTest: 0.515623 (exp: 0.5)
    Means and Covariances:
    param0 = mean 0.505875 with covariance 0.093983
    
    Creating smcTable_Test_identity1_2.txt.
    Output will be written to Test_identity1/Sim2/.
    Preparing simulations. After simulations have finished, rerun this script to continue.
    Running make in the build directory /Users/weinhartt/Code/Calibration/cmake-build-debug/Drivers/Calibration
    Runs started. Rerun this script for analysis.
    Starting analysis of iteration 2
    Skipping DEM simulations, read in pre-existing simulation data now
    The optimal parameters for Test_identity1 in iteration 2 are:
    Best params for Test_identity1 in iteration 2 (error 0.011758):
    param0: 0.505879
    data: CalibrationTest: 0.505879 (exp: 0.5)
    Means and Covariances:
    param0 = mean 0.509932 with covariance 0.032157
    Finished Calibration
    ````
  - As you can see, it slowly converges to the value `0.5` for `param0`.

## A second test example

For a second test example, see `calibrateTest_calibration44.txt`.

## A realistic example

Let's say we have the following calibration data for Material A:
 - A static angle of repose 
 - A dynamic angle of repose

We chose to work with spherical particles and with the LinearViscoelasticFrictionReversibleAdhesiveSpecies contact model. Thus, your particle and contact model parameters are:
 - particle size distribution
 - particle density
 - normal stiffness and dissipation coefficient
 - sliding, rolling and torsion stiffness and dissipation coefficient
 - sliding, rolling, and torsion friction coefficient
 - adhesion stiffness and maximum adhesive force

Before we calibrate, we reduce the amount of free parameters using the following assumptions:
 - Set the particle density and particle size distribution using experimentally determined values
 - Set the normal, sliding, rolling and torsion stiffness and dissipation coefficient by setting an appropriate collision time and a to-be-determined restitution coefficient (usually, we use the same values for normal, sliding, rolling and torsion)
 - Assume torsion friction is negligible (zero)
 - Assume adhesion stiffness equals the normal stiffness, and set the maximum adhesive force via the (minimum) bond number, i.e. f_a^max = Bo m_min g

Now the remaining free contact model parameters are:
 - Restitution coefficient
 - Sliding friction
 - Rolling friction
 - Bond number

Now we have four parameters to calibrate, and two calibration data values. Next, I will describe the step-by-step procedure to run the calibration suite:
 - Make a file named calibrationMaterialA.txt and place it in the Drivers/Calibration/py folder. For this example, the file already exists, please take a look. The file name contains the name of the material to be calibrated (no spaces or underscores), prepended "calibrate" and appended ".txt". It contains three sections, each started by a header line: 
   - Section 1: fixed material parameters
   - Section 2: variable parameters
   - Section 3: experimental data
 - Open a command line and enter "./calibrate.py --material MaterialA". This will start the calibration procedure. However, it takes many days to complete, so it's better to run this in parallel on a cluster.
 - Install MercuryDPM on the msm3 cluster, run "./calibrate.py --material MaterialA -node 01 -cores 40" and follow the instructions (Todo: the cluster now has slurm, so we need to update this section). This will run the calibration on forty parallel cores on node01.
 - The calibration script will create output in the folder MaterialA. 
   - The characterisation data will be stored in subdirectory Exp/ 
   - The 0-th iteration data will be stored in subdirectory Sim0/, the 1-st iteration in Sim1/, etc
 - Check that the output created by the calibrate.py script is sensible. It should look like this:
   ````
   Running on 40 core(s) on node(s) 01
   Material: MaterialA
   Parameters to be identified: restitutionCoefficient, slidingFriction, rollingFriction, relativeCohesionStiffness
   Parameter ranges: [0.2, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 2.0]
   Executables: ShearCellBu
   Data values: [770.0, 532.0, 327.0]
   Weights: [1.0, 1.0, 1.0]
   Command line parameters: -speciesType MaterialA -restitutionCoefficient %s -slidingFriction %s -rollingFriction %s -relativeCohesionStiffness %s -normalStress 1027 827 527 226
   Build directory for executables: Drivers/Calibration
   Covariance: 1
   Goal for effective sample size: 0.500000
   Number of Gaussian components: 2
   Number of samples: 40
   Prior weight: 0.500000
   
   Creating smcTable_MaterialA_0.txt.
   Output will be written to MaterialA/Sim0/.
   Preparing simulations. After simulations have finished, rerun this script to continue.
   Running make in the build directory /storage1/usr/people/weinhartt/Code/Lab/build/Drivers/USER/MercuryLab/Buttercup/Calibration
   Number of commands 40, per script 2, number of scripts: 20
   Adjusting script for msm3
   1) cd '/storage1/usr/people/weinhartt/Code/Lab/Drivers/USER/MercuryLab/Buttercup/Calibration/py/MaterialA/Sim0/'
   2) Adjust node numbers run.sh.
   3) Execute 'source ./run.sh'
      Runs started. Rerun this script for analysis.
   ````
 - Repeat the step above until all simulations have finished (this might take a day or two). 
 - Now run "./calibrate.py --material MaterialA -node 01 -cores 40" again to start the next iteration.
 
