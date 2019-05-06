/*
    *** COMPRESSION TEST ***
 Particles are placed inside a cylindrical casing and let settle under gravity.
 When the system is at rest a piston is loaded and starts compressing the particles.
 The values of pressure towards the piston, the compression length and ther quantities are printed at each time step.
 
 *** IMPORTANT ***
 The interaction law is purely cohesive with no adhesion.
 If needed uncomment the proper parts and substitute
 LinearPlasticViscoelasticFrictionSpecies -> LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies
 
 LAST UPDATE: 20/10/17
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
#include "MultiStepCompressionTest.h"

// computes the mean spring stiffness
double computeMean(int n, double *arrayPointer)
{
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += arrayPointer[i];
    
    return mean /= n;
}

// computes the mean absolute deviation
double computeMeanAbsoluteDeviation(int n, double *arrayPointer)
{
    double deviation = 0.0;
    double mean = computeMean(n, arrayPointer);
    
    for (int i = 0; i < n; i++) deviation += fabs(arrayPointer[i] - mean);
    
    return deviation /= n;
}

int main(int argc, char *argv[])
{
    // geometry size setting
    double cylinderRadius = 0.0125;
    double cylinderHeight = 0.019;
    
    // particle size setting
    double cylinderToParticleRadiusRatio = 40.0;
    double particleRadius = cylinderRadius/cylinderToParticleRadiusRatio;
    
    // fixed hysteretic interaction parameters
    double kWall = 3000.0;
    double k2max = 2000.0;
    double k1 = 400.0;
    double kc = 40.0;
    double phi = 0.05;
    
    // height and pressure arrays for the data calibration procedure
    //int nDataPoints = 16;
    //double heightArray[nDataPoints] = {0.0185328, 0.0184762, 0.0184496, 0.0184252, 0.0184058, 0.0183886, 0.0183917, 0.0184044, 0.0184047, 0.0183935, 0.0183803, 0.0183798, 0.018384, 0.0183902, 0.0183958, 0.0184039};
    //double pressureArray[nDataPoints] = {4790.42, 5763.35, 6749.81, 7771.34, 8751.82, 9717.71, 7929.76, 6001.37, 5993.35, 8013.84, 9934.54, 9056.17, 8030.74, 7008.09, 5953.65, 4975.92};
    //const int nDataPoints = 10;
    //double heightArray[nDataPoints] = {0.0185328, 0.0184762, 0.0184496, 0.0184252, 0.0184058, 0.0183886, 0.0183917, 0.0184044, 0.0184047, 0.0183935};
    //double pressureArray[nDataPoints] = {4790.42, 5763.35, 6749.81, 7771.34, 8751.82, 9717.71, 7929.76, 6001.37, 5993.35, 8013.84};
    const int nDataPoints = 6;
    double heightArray[nDataPoints] = {0.0185328, 0.0184762, 0.0184496, 0.0184252, 0.0184058, 0.0183886};
    double pressureArray[nDataPoints] = {4790.42, 5763.35, 6749.81, 7771.34, 8751.82, 9717.71};
    
    // calibration routine parameters
    double minimumDispersionTolerance = 0.05;
    int maximumNumberRoutineTrials = 10;
    int runNumber = 0;
    double *initialKArray, *disArray, *calibratedKArray, *meanCalibratedKArray;
    double *kArray;

    // friction parameters
	double muSBigBig = 0.16;
	double muRBigBig = 0.05;
    
    // initialization of the calibration routine arrays
    initialKArray = new double[maximumNumberRoutineTrials];
    disArray = new double[maximumNumberRoutineTrials];
    meanCalibratedKArray = new double[maximumNumberRoutineTrials];
    calibratedKArray = new double[nDataPoints];
    
    kArray = new double[nDataPoints];
    
    for (int i=0; i<maximumNumberRoutineTrials; i++)
    {
        initialKArray[i] = 0.0;
        disArray[i] = 1.0;
    }
    
    initialKArray[0] = k1;
    
    // name setting variables
    std::ostringstream name;
//     std::ostringstream kdataFilename;
    
    name.str("");
    name.clear();
//     kdataFilename.str("");
//     kdataFilename.clear();

    do
    {
        // setting the compression test elastic constant
        k1 = initialKArray[runNumber];
        
        // creation of the model
        CompressionTest_parameterCalibrationRoutine cTestRoutine;
        
        // sets the bool restartedFile true or false accorging to the argc value
        if (argc > 1) {cTestRoutine.restartedFile = true;}
        else {cTestRoutine.restartedFile = false;}
        
        // reads the .restart file from command line (if any)
        cTestRoutine.readArguments(argc, argv);
                
        // sets simulation parameters
        cTestRoutine.setTimeStep(1.0e-6);
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
        cTestRoutine.setParticlePlasticProperties(k1, k2max, kc, phi);
        //    cTestRoutine.setParticlesAdhesionProperties(0.0, 0.0);
        
        // sets the flag for automatic k2 scaling during data-driven calibration (default is false)
        cTestRoutine.setK2automaticScaling(false);
        
        // sets the particle-wall friction coefficients
        cTestRoutine.setParticleWallSlidingFrictionCoeff(0.3, 0.3);
        cTestRoutine.setParticleWallRollingFrictionCoeff(0.01, 0.01);
        cTestRoutine.setParticleWallTorsionFrictionCoeff(0.0, 0.0);
        
        // sets the particle-particle friction coefficients
        cTestRoutine.setParticleParticleSlidingFrictionCoeff(muSBigBig, 0.0, 0.0);
        cTestRoutine.setParticleParticleRollingFrictionCoeff(muRBigBig, 0.0, 0.0);
        cTestRoutine.setParticleParticleTorsionFrictionCoeff(0.0, 0.0, 0.0);
        
        // sets the particle restitution coefficients
        cTestRoutine.setParticleWallRestitutionCoeff(0.7, 0.7);
        cTestRoutine.setParticleParticleRestitutionCoeff(0.5, 0.5, 0.5);
        
        // sets the casing dimensions
        cTestRoutine.setCasingDimensions(cylinderRadius, cylinderHeight);
        
        // sets the bulk packing fraction of the powder prior to compression
        cTestRoutine.setbulkPackingFractionPreCompression(0.60);
        
        cTestRoutine.setDataPoints(nDataPoints, heightArray, pressureArray, 0.01);
        
        // sets additional built-in arguments for the xballs visualization
        cTestRoutine.setXBallsAdditionalArguments("-h 800 -p 10 -o 200 -3dturn 3");
        
        // passes the pointer to the calibrated k array for the computation of the mean and the deviation
        cTestRoutine.kArraySetter(nDataPoints, kArray);
        
        // name setting
        if (cTestRoutine.restartedFile)
        {
			if (nDataPoints == 6)
			{
				std::cout.unsetf(std::ios::floatfield);
				name << "CompressionTest_dataDriven_CARRIER_LOADING_runNumber_" << runNumber << "_rRatio_" << std::fixed << std::setprecision(0) << cylinderToParticleRadiusRatio
				<< "_k1_" << k1 << "_k2max_" << k2max << "_kc_" << kc << "_phi_" << std::fixed << std::setprecision(2) << phi
				<< "_muS_" << muSBigBig << "_muR_" << muRBigBig;
				cTestRoutine.setName(name.str());
			}
			else if (nDataPoints == 10)
			{
				std::cout.unsetf(std::ios::floatfield);
				name << "CompressionTest_dataDriven_CARRIER_RELOADING_runNumber_" << runNumber << "_rRatio_" << std::fixed << std::setprecision(0) << cylinderToParticleRadiusRatio
				<< "_k1_" << k1 << "_k2max_" << k2max << "_kc_" << kc << "_phi_" << std::fixed << std::setprecision(2) << phi
				<< "_muS_" << muSBigBig << "_muR_" << muRBigBig;
				cTestRoutine.setName(name.str());
			}
			else if (nDataPoints == 16)
			{
				std::cout.unsetf(std::ios::floatfield);
				name << "CompressionTest_dataDriven_CARRIER_FULL_runNumber_" << runNumber << "_rRatio_" << std::fixed << std::setprecision(0) << cylinderToParticleRadiusRatio
				<< "_k1_" << k1 << "_k2max_" << k2max << "_kc_" << kc << "_phi_" << std::fixed << std::setprecision(2) << phi
				<< "_muS_" << muSBigBig << "_muR_" << muRBigBig;
				cTestRoutine.setName(name.str());
			}
            
            std::cout << "This file has been restarted.\nThe new name of the associated output of this simulation will be:\n" << name.str() << std::endl;
        }
        else
        {
            std::cout.unsetf(std::ios::floatfield);
            name << "CompressionTest_SETTLING_CARRIER_rRatio_" << std::fixed << std::setprecision(0) << cylinderToParticleRadiusRatio
				<< "_k1_" << k1 << "_k2max_" << k2max << "_kc_" << kc << "_phi_" << std::fixed << std::setprecision(2) << phi
				<< "_muS_" << muSBigBig << "_muR_" << muRBigBig;
//             kdataFilename  << "CompressionTest_parameterEstimation_RESTARTFILE_CARRIER_rRatio_" << cylinderToParticleRadiusRatio << "_k1_" << k1 << "_k2max_" << k2max << "_kc_" << kc << "_phi_" << std::fixed << std::setprecision(2) << phi << ".kdata";
            cTestRoutine.setName(name.str());
            
            std::cout << "The new name of the associated output of this simulation will be:\n" << name.str() << std::endl;
        }
        
        cTestRoutine.solve();
        
        // computing the mean value and the absolute deviation of the calibrated stiffnesses
        meanCalibratedKArray[runNumber] = computeMean(nDataPoints, kArray);
        disArray[runNumber] = computeMeanAbsoluteDeviation(nDataPoints, kArray);
        
        std::cout << "k0 = " << initialKArray[runNumber] << "\n";
        std::cout << "mean(k) = " << meanCalibratedKArray[runNumber] << "\n";
        std::cout << "|disp(k)| = " << disArray[runNumber] << "\n";
        std::cout << "|disp(k)|/mean(k) = " << disArray[runNumber]/meanCalibratedKArray[runNumber] << "\n";
        std::cout << "|mean(k) - k0|: " << fabs(meanCalibratedKArray[runNumber] - initialKArray[runNumber]) << "\n";
        std::cout << "|mean(k) - k0|/k0: " << fabs(meanCalibratedKArray[runNumber] - initialKArray[runNumber])/initialKArray[runNumber] << "\n";
        std::cout << "\n";
        
        // update the initial k value
        if (runNumber < maximumNumberRoutineTrials - 1)
        {
            initialKArray[runNumber + 1] = 0.5*(meanCalibratedKArray[runNumber] + initialKArray[runNumber]);
            std::cout << "New starting k0: " << initialKArray[runNumber + 1] << std::endl;
        }
        
        // clearing the streams and incrementing the run number
        name.str("");
        name.clear();
//         kdataFilename.str("");
//         kdataFilename.clear();
        runNumber++;
    }
    while (runNumber < maximumNumberRoutineTrials && !(disArray[runNumber - 1]/meanCalibratedKArray[runNumber - 1] < minimumDispersionTolerance && fabs(meanCalibratedKArray[runNumber - 1] - initialKArray[runNumber - 1])/initialKArray[runNumber - 1] < minimumDispersionTolerance));
    
    if (disArray[runNumber - 1]/meanCalibratedKArray[runNumber - 1] > minimumDispersionTolerance || fabs(meanCalibratedKArray[runNumber - 1] - initialKArray[runNumber - 1])/initialKArray[runNumber - 1] > minimumDispersionTolerance)
    {
        std::cout << "\nFAILED TO CONVERGE TO AN ACCEPTABLE DISPERSION\n";
    }
    else
    {
        std::cout << "\nOptimal elastic spring stiffness found: " << meanCalibratedKArray[runNumber - 1] << ", with a dispersion of " << disArray[runNumber - 1] << " after " << runNumber << " cycles.\n";
    }
    
    return 0;
}







