/*
    *** COMPRESSION TEST ***
 Particles are placed inside a cylindrical casing and let settle under gravity.
 When the system is at rest a piston is loaded and starts compressing the particles.
 The values of pressure towards the piston, the compression length and ther quantities are printed at each time step.
 
 *** IMPORTANT ***
 The interaction law is purely cohesive with no adhesion.
 If needed uncomment the proper parts and substitute
 LinearPlasticViscoelasticFrictionSpecies -> LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies
 
 LAST UPDATE: 07/12/17
*/

#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include "CompressionPiston.h"
#include <math.h>
#include <fstream>
#include "MultiStepCompressionTestDataDrivenNewMethod.h"

int main(int argc, char *argv[])
{
    // geometry size setting
    double cylinderRadius = 0.0125;
    double cylinderHeight = 0.019;
    
    // particle size setting
    double cylinderToParticleRadiusRatio = 40.0;
    double particleRadius = cylinderRadius/cylinderToParticleRadiusRatio;
    
    // fixed hysteretic interaction parameters
    double kWall = 3000;
    double k2max = 2000;
    double k1_0 = 0.1*k2max;
    double kC_0 = 0.1*k1_0;
    double phi = 0.05;
    
    // height and pressure arrays for the data calibration procedure
    const int nDataPoints = 16;
    double heightArray[nDataPoints] = {0.0185328, 0.0184762, 0.0184496, 0.0184252, 0.0184058, 0.0183886, 0.0183917, 0.0184044, 0.0184047, 0.0183935, 0.0183803, 0.0183798, 0.018384, 0.0183902, 0.0183958, 0.0184039};
    double pressureArray[nDataPoints] = {4790.42, 5763.35, 6749.81, 7771.34, 8751.82, 9717.71, 7929.76, 6001.37, 5993.35, 8013.84, 9934.54, 9056.17, 8030.74, 7008.09, 5953.65, 4975.92};
    
    // calibration routine parameters
    double minimumDispersionTolerance = 0.10;
    int maximumNumberRoutineTrials = 50;
    int runNumber = 0;
    double initialK1Array[maximumNumberRoutineTrials];
    double initialKcArray[maximumNumberRoutineTrials];
       
    // sets the initial values of the stiffnesses to 0 and the dispersion to 1
    for (int i=0; i<maximumNumberRoutineTrials; i++)
    {
        initialK1Array[i] = 0.0;
        initialKcArray[i] = 0.0;
    }
    
    // the initial value of the array is set to the initial specified value
    initialK1Array[0] = k1_0;
    initialKcArray[0] = kC_0;
    
    // bools for determining which calibration branch is run
    bool *calibrationDone;
    calibrationDone = new bool[1];
    calibrationDone[0] = {false};
    
    // name setting variables
    std::ostringstream name;
    name.str("");
    name.clear();

    do
    {
        // creation of the model
        CompressionTest_parameterCalibrationRoutine cTestRoutine;
        
        // sets the bool restartedFile true or false accorging to the argc value
        if (argc > 1) {cTestRoutine.restartedFile = true;}
        else {cTestRoutine.restartedFile = false;}
        
        // reads the .restart file from command line (if any)
        cTestRoutine.readArguments(argc, argv);
        
//        cTestRoutine.readRestartFile(argv[2]);
        
        // name setting
        if (cTestRoutine.restartedFile)
        {
            std::cout.unsetf(std::ios::floatfield);
            name << "CompressionTest_NEWROUTINE_CARRIER_Nrun_" << runNumber << "_rRatio_" << std::fixed << std::setprecision(0) <<  cylinderToParticleRadiusRatio << "_k1_" << initialK1Array[runNumber] << "_k2max_" << k2max << "_kc_" << initialKcArray[runNumber] << "_phi_" << std::fixed << std::setprecision(2) << phi;
            cTestRoutine.setName(name.str());
            
            std::cout << "This file has been restarted.\nThe new name of the associated output of this simulation will be:\n" << name.str() << std::endl;
        }
        else
        {
            std::cout << "\n*** WARNING: NO INITIAL RESTART FILE SPECIFIED. TERMINATING.\n";
            exit(0);
        }
        
        // sets simulation parameters
        cTestRoutine.setTimeStep(3.0e-6);
        cTestRoutine.setTimeMax(100.0);
        cTestRoutine.setGravity(Vec3D(0., 0., -9.81));
        cTestRoutine.setSystemDimensions(3);
        cTestRoutine.setVerbose(true);
        
        // global simulation parameters
        cTestRoutine.setCdatOutputTimeInterval(0.005);
        cTestRoutine.setEnergyRatioTolerance(1.0e-4);
        cTestRoutine.setSaveCount(0.01/cTestRoutine.getTimeStep());
        
        // sets the particle intrinsic properties
        cTestRoutine.setParticleDensity(1452.7, 1452.7);
        cTestRoutine.setParticleRadiusAndDispersity(particleRadius, 0.10, particleRadius, 0.10);
        
        // sets the small-to-big total mass ratio
        cTestRoutine.setTotalMassRatio(0.01);
        
        // sets the stiffnesses
        cTestRoutine.setWallStiffness(kWall);
        cTestRoutine.setParticlePlasticProperties(initialK1Array[runNumber], k2max, initialKcArray[runNumber], phi);
        //    cTestRoutine.setParticlesAdhesionProperties(0.0, 0.0);
        
        // sets the particle-wall friction coefficients
        cTestRoutine.setParticleWallSlidingFrictionCoeff(0.3, 0.3);
        cTestRoutine.setParticleWallRollingFrictionCoeff(0.01, 0.01);
        cTestRoutine.setParticleWallTorsionFrictionCoeff(0.0, 0.0);
        
        // sets the particle-particle friction coefficients
        cTestRoutine.setParticleParticleSlidingFrictionCoeff(0.5, 0.5, 0.5);
        cTestRoutine.setParticleParticleRollingFrictionCoeff(0.3, 0.3, 0.3);
        cTestRoutine.setParticleParticleTorsionFrictionCoeff(0.0, 0.0, 0.0);
        
        // sets the particle restitution coefficients
        cTestRoutine.setParticleWallRestitutionCoeff(0.7, 0.7);
        cTestRoutine.setParticleParticleRestitutionCoeff(0.5, 0.5, 0.5);
        
        // sets the casing dimensions
        cTestRoutine.setCasingDimensions(cylinderRadius, cylinderHeight);
        
        // sets the bulk packing fraction of the powder prior to compression
        cTestRoutine.setbulkPackingFractionPreCompression(0.60);
	

        // NEW STUFF HERE -  -  -  -  -  -  -  -  -
        

        // sets the bool for teh calibration done or not
        cTestRoutine.setCalibrationStatus(calibrationDone);
	
        // sets teh data points 
        cTestRoutine.setCalibrationDataPoints(nDataPoints, heightArray, pressureArray, 0.01);

        // sets the convergence threshold value
        cTestRoutine.setPressureConvergengeThresholdValue(minimumDispersionTolerance);
        
        // passes the run number (to avoid updating the pointers to the initialArray[runNumber + 1]) and the max
        // *** WARNING! THIS MUST BE ALLOCATED BEFORE setInitialStiffnessPointer
        cTestRoutine.setMaxRunNumber(maximumNumberRoutineTrials);
        cTestRoutine.setRunNumber(runNumber);

        // passes the pointer to the initial k arrays
        cTestRoutine.setInitialStiffnessPointer(initialK1Array, initialKcArray);
        
        // sets additional built-in arguments for the xballs visualization
        cTestRoutine.setXBallsAdditionalArguments("-h 800 -p 10 -o 200 -3dturn 3");

        cTestRoutine.solve();
        
        // clearing the streams and incrementing the run number
        name.str("");
        name.clear();
        runNumber++;
    }
    while (runNumber < maximumNumberRoutineTrials && !calibrationDone);
    
    return 0;
}







