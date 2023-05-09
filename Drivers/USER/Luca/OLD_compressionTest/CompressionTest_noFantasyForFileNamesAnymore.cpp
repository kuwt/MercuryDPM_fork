/*
    *** COMPRESSION TEST ***
 Particles are placed inside a cylindrical casing and let settle under gravity.
 When the system is at rest a piston is loaded and starts compressing the particles.
 The values of pressure towards the piston, the compression length and ther quantities are printed at each time step.
 
 *** IMPORTANT ***
 The interaction law is purely cohesive with no adhesion.
 If needed uncomment the proper parts and substitute
 LinearPlasticViscoelasticFrictionSpecies -> LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies
 
 LAST UPDATE: 18/10/17
*/

#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <math.h>
#include <fstream>
#include "MultiStepCompressionTest.h"


/*
 1 - run data analysis for loading branch
 2 - find the ideal k1
 3 - run height-driven simulation with the k1 found
 4 - if (OK) {loadingBranchDone = true;} else {k1_0 = k1;} 
 5 - run data analysis for reloading branch
 6 - find the ideal k1
 7 - run height-driven simulation with the k1 found
 8 - if (OK) {loadingBranchDone = true;} else {}
 
 
 */

class Pippo : public Mercury3D
{
private:
    
    void setupInitialConditions() override
    {
        
    }
    
    void actionsOnRestart() override
    {
        
    }
    
    void actionsAfterTimeStep() override
    {
        if (!loadingBranchDone && !reloadingBranchDone && !unloadingBranchDone)
        {
            
        }
    }
    
    void actionAfterSolve()
    {
        delete[] heightDataAnalysisPoints;
        delete[] pressureDataAnalysisPoints;
    }
    
public:
    
    // reads the data points and initializes the height-pressure array for the loading branch
    void setLoadingDataPoints(int n, double *heightArrayPointer, double *pressureArrayPointer, double duration)
    {
        loadingCompressionInstances = n;
        loadingCompressionInstanceDuration = duration;
        heightStepsLoading = new double[n];
        pressureStepsLoading = new double[n];
        calibratedLoadingK1 = new double[n];
        calibratedLoadingK2 = new double[n];
        calibratedLoadingKC = new double[n];
        calibratedLoadingPhi = new double[n];
        
        for (int i = 0; i < n; i++)
        {
            heightStepsLoading[i] = heightArrayPointer[i];
            pressureStepsLoading[i] = pressureArrayPointer[i];
        }
        
        if (verbose)
        {
            std::cout << "\tLoading height-pressure steps:\n";
            for (int i = 0 ; i < n; i++) std::cout << "\t( " << heightStepsLoading[i] << " , " << pressureStepsLoading[i] << " )\n";
            std::cout << std::endl;
        }
    }
    
    // reads the data points and initializes the height-pressure array for the reloading branch
    void setReloadingDataPoints(int n, double *heightArrayPointer, double *pressureArrayPointer, double duration)
    {
        reloadingCompressionInstances = n;
        reloadingCompressionInstanceDuration = duration;
        heightStepsReloading = new double[n];
        pressureStepsReloading = new double[n];
        calibratedReloadingK1 = new double[n];
        calibratedReloadingK2 = new double[n];
        calibratedReloadingKC = new double[n];
        calibratedReloadingPhi = new double[n];
        
        for (int i = 0; i < n; i++)
        {
            heightStepsReloading[i] = heightArrayPointer[i];
            pressureStepsReloading[i] = pressureArrayPointer[i];
        }
        
        if (verbose)
        {
            std::cout << "\tLoading height-pressure steps:\n";
            for (int i = 0 ; i < n; i++) std::cout << "\t( " << heightStepsReloading[i] << " , " << pressureStepsReloading[i] << " )\n";
            std::cout << std::endl;
        }
    }
    
    // reads the data points and initializes the height-pressure array for the unloading branch
    void setUnloadingDataPoints(int n, double *heightArrayPointer, double *pressureArrayPointer, double duration)
    {
        unloadingCompressionInstances = n;
        unloadingCompressionInstanceDuration = duration;
        heightStepsUnloading = new double[n];
        pressureStepsUnloading = new double[n];
        calibratedUnloadingK1 = new double[n];
        calibratedUnloadingK2 = new double[n];
        calibratedUnloadingKC = new double[n];
        calibratedUnloadingPhi = new double[n];
        
        for (int i = 0; i < n; i++)
        {
            heightStepsUnloading[i] = heightArrayPointer[i];
            pressureStepsUnloading[i] = pressureArrayPointer[i];
        }
        
        if (verbose)
        {
            std::cout << "\tLoading height-pressure steps:\n";
            for (int i = 0 ; i < n; i++) std::cout << "\t( " << heightStepsUnloading[i] << " , " << pressureStepsUnloading[i] << " )\n";
            std::cout << std::endl;
        }
    }
    
    
    
    void setDataAnalysisPoints()
    {
        // if in the loading branch allocates only the loading data points
        if (!loadingBranchDone && !reloadingBranchDone && !unloadingBranchDone)
        {
            nDataAnalysisPoints = nLoadingCompressionInstances;
            heightDataAnalysisPoints = new double[nDataAnalysisPoints];
            pressureDataAnalysisPoints = new double[nDataAnalysisPoints];
            
            for (int i = 0; i < nLoadingCompressionInstances; i++)
            {
                heightDataAnalysisPoints[i] = heightStepsLoading[i];
                pressureDataAnalysisPoints[i] = pressureStepsLoading[i];
            }
        }
        
        // if in the reloading branch allocates the loading and teh reloading data points
        if (loadingBranchDone && !reloadingBranchDone && !unloadingBranchDone)
        {
            nDataAnalysisPoints = nLoadingCompressionInstances + nReloadingCompressionInstances;
            heightDataAnalysisPoints = new double[nDataAnalysisPoints];
            pressureDataAnalysisPoints = new double[nDataAnalysisPoints];
            
            for (int i = 0; i < nLoadingCompressionInstances; i++)
            {
                heightDataAnalysisPoints[i] = heightStepsLoading[i];
                pressureDataAnalysisPoints[i] = pressureStepsLoading[i];
            }
            for (int i = 0; i < nReloadingCompressionInstances; i++)
            {
                heightDataAnalysisPoints[i + nLoadingCompressionInstances] = heightStepsReloading[i];
                pressureDataAnalysisPoints[i + nLoadingCompressionInstances] = pressureStepsReloading[i];
            }
        }
        
        // if in the unloading branch allocates all the data points
        if (loadingBranchDone && reloadingBranchDone && !unloadingBranchDone)
        {
            nDataAnalysisPoints = nLoadingCompressionInstances + nReloadingCompressionInstances + nUnloadingCompressionInstances;
            heightDataAnalysisPoints = new double[nDataAnalysisPoints];
            pressureDataAnalysisPoints = new double[nDataAnalysisPoints];
            
            for (int i = 0; i < nLoadingCompressionInstances; i++)
            {
                heightDataAnalysisPoints[i] = heightStepsLoading[i];
                pressureDataAnalysisPoints[i] = pressureStepsLoading[i];
            }
            for (int i = 0; i < nReloadingCompressionInstances; i++)
            {
                heightDataAnalysisPoints[i + nLoadingCompressionInstances] = heightStepsReloading[i];
                pressureDataAnalysisPoints[i + nLoadingCompressionInstances] = pressureStepsReloading[i];
            }
            for (int i = 0; i < nUnloadingCompressionInstances; i++)
            {
                heightDataAnalysisPoints[i + nLoadingCompressionInstances + nReloadingCompressionInstances] = heightStepsReloading[i];
                pressureDataAnalysisPoints[i + nLoadingCompressionInstances + nReloadingCompressionInstances] = pressureStepsReloading[i];
            }
        }
    }
    
    // calibrates the loading stiffness
    void calibrateLoadingStiffness()
    {
        // halves and reverses the piston velocity if needed
        if ((pistonVelocity < 0.0 && pistonPointer -> getHeight() < targetHeight) || (pistonVelocity > 0.0 && pistonPointer -> getHeight() > targetHeight))
        {
            pistonVelocity *= 0.5;
            pistonPointer -> setVelocity(pistonVelocity);
        }
        
        // moving the piston when the target height is not met
        if (!heightTargetMet && fabs((pistonPointer -> getHeight() - targetHeight)/targetHeight) > 1.0e-4)
        {
            pistonPointer -> movePiston(getTimeStep());
        }
        
        // when the target height is met updates the bool
        if (!heightTargetMet && fabs((pistonPointer -> getHeight() - targetHeight)/targetHeight) < 1.0e-4)
        {
            std::cout << "Height target reached. Start tuning the spring stiffness...\n";
            heightTargetMet = true;
        }
        
        // if the height is met and the pressure is not met tunes k1
        if (heightTargetMet && !pressureTargetMet)
        {
            if (fabs((pistonPressure - targetPressure)/targetPressure) > 1.0e-4)
            {
                adjustK1();
            }
            else
            {
                std::cout << "Pressure target reached, now chillin' for " << compressionInstanceDuration << " seconds...\n";
                timePlaceholder = getTime();
                pressureTargetMet = true;
            }
        }
        
        // when both height and pressure are met waits and updates the target values
        if (heightTargetMet && pressureTargetMet && getTime() > timePlaceholder + compressionInstanceDuration)
        {
            calibratedK1[compressionStep] = particleStiffnessBig;
            calibratedK2[compressionStep] = particleMaxUnloadingStiffness;
            calibratedKC[compressionStep] = particleCohesiveStiffness;
            calibratedPhi[compressionStep] = particlePlasticityDepth;
            compressionStep++;
            
            if (compressionStep < nDataAnalysisPoints)
            {
                targetHeight = loadingHeightSteps[compressionStep];
                targetPressure = loadingPressureSteps[compressionStep];
                std::cout << "New target height: " << targetHeight << "\n";
                std::cout << "New target pressure: " << targetPressure << "\n";
                resetPistonVelocity();
                heightTargetMet = false;
                pressureTargetMet = false;
            }
            else
            {
                std::cout << "Loading cycles finished...\n";
                setTimeMax(getTime() + compressionInstanceDuration);
                stage++;
            }
        }
    }
    
    void calibrateCohesiveStiffness()
    {
        // halves and reverses the piston velocity if needed
        if ((pistonVelocity < 0.0 && pistonPointer -> getHeight() < targetHeight) || (pistonVelocity > 0.0 && pistonPointer -> getHeight() > targetHeight))
        {
            pistonVelocity *= 0.5;
            pistonPointer -> setVelocity(pistonVelocity);
        }
        
        // moving the piston when the target height is not met
        if (!heightTargetMet && fabs((pistonPointer -> getHeight() - targetHeight)/targetHeight) > 1.0e-4)
        {
            pistonPointer -> movePiston(getTimeStep());
        }
        
        // when the target height is met updates the bool
        if (!heightTargetMet && fabs((pistonPointer -> getHeight() - targetHeight)/targetHeight) < 1.0e-4)
        {
            std::cout << "Height target reached. Start tuning the spring stiffness...\n";
            heightTargetMet = true;
        }
        
        // if the height is met and the pressure is not met tunes k1
        if (heightTargetMet && !pressureTargetMet)
        {
            if (fabs((pistonPressure - targetPressure)/targetPressure) > 1.0e-4)
            {
                adjustKC();
            }
            else
            {
                std::cout << "Pressure target reached, now chillin' for " << compressionInstanceDuration << " seconds...\n";
                timePlaceholder = getTime();
                pressureTargetMet = true;
            }
        }
        
        // when both height and pressure are met waits and updates the target values
        if (heightTargetMet && pressureTargetMet && getTime() > timePlaceholder + compressionInstanceDuration)
        {
            calibratedK1[compressionStep] = particleStiffnessBig;
            calibratedK2[compressionStep] = particleMaxUnloadingStiffness;
            calibratedKC[compressionStep] = particleCohesiveStiffness;
            calibratedPhi[compressionStep] = particlePlasticityDepth;
            compressionStep++;
            
            if (compressionStep < unloadingCompressionInstances)
            {
                targetHeight = unloadingHeightSteps[compressionStep];
                targetPressure = unloadingPressureSteps[compressionStep];
                std::cout << "New target height: " << targetHeight << "\n";
                std::cout << "New target pressure: " << targetPressure << "\n";
                resetPistonVelocity();
                heightTargetMet = false;
                pressureTargetMet = false;
            }
            else
            {
                std::cout << "Unloading cycles finished...\n";
                setTimeMax(getTime() + compressionInstanceDuration);
                stage++;
            }
        }
        
    }
    
    
    int nDataAnalysisPoints;
    double *heightDataAnalysisPoints;
    double *pressureDataAnalysisPoints;
    
    int nLoadingCompressionInstances;
    double loadingCompressionInstanceDuration;
    int nReloadingCompressionInstances;
    double reloadingCompressionInstanceDuration;
    int nUnloadingCompressionInstances;
    double unloadingCompressionInstanceDuration;
    
    bool loadingBranchDone = false;
    bool reloadingBranchDone = false;
    bool unloadingBranchDone = false;
};














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
    double k1_0;
    double k2max = 3000;
    double kc_0 = 0.0;
    double phi = 0.03;
    
    k1_0 = 0.05*k2max;
    
    // height and pressure arrays for the data calibration procedure
    // LOADING
    const int nLoadingDataPoints = 6;
    double heightArrayLoading[nLoadingDataPoints] = {0.0185328, 0.0184762, 0.0184496, 0.0184252, 0.0184058, 0.0183886};
    double pressureArrayLoading[nLoadingDataPoints] = {4790.42, 5763.35, 6749.81, 7771.34, 8751.82, 9717.71};
    // RELOADING
    const int nReloadingDataPoints = 4;
    double heightArrayReloading[nReloadingDataPoints] = {0.0183917, 0.0184044, 0.0184047, 0.0183935};
    double pressureArrayReloading[nReloadingDataPoints] = {7929.76, 6001.37, 5993.35, 8013.84};
    // UNLOADING
    const int nUnloadingDataPoints = 6;
    double heightArrayUnloading[nUnloadingDataPoints] = {0.0183803, 0.0183798, 0.018384, 0.0183902, 0.0183958, 0.0184039};
    double pressureArrayUnloading[nUnloadingDataPoints] = {9934.54, 9056.17, 8030.74, 7008.09, 5953.65, 4975.92};
    
    const int nDataPoints = nLoadingDataPoints + nReloadingDataPoints + nUnloadingDataPoints;
    
    // calibration routine parameters
    double minimumDispersionTolerance = 0.05;
    int maximumNumberRoutineTrials = 10;
    int runNumber = 0;
    double *initialK1Array, *dispersionK1Array, *calibratedK1Array, *meanCalibratedK1Array;
    
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
    
    initialKArray[0] = k1_1;
    
    // name setting variables
    std::ostringstream name;
    std::ostringstream kdataFilename;
    
    name.str("");
    name.clear();
    kdataFilename.str("");
    kdataFilename.clear();
    
    while (runNumber - 1 < maximumNumberRoutineTrials && !(disArray[runNumber - 1]/meanCalibratedKArray[runNumber - 1] < minimumDispersionTolerance && fabs(meanCalibratedKArray[runNumber - 1] - initialKArray[runNumber - 1])/initialKArray[runNumber - 1] < minimumDispersionTolerance))
    {
        // setting the compression test elastic constant
        k1 = initialKArray[runNumber];
        
        // creation of the model
        CompressionTest_parameterCalibrationRoutine cTestRoutine;
        
        // sets the bool restartedFile true or false accorging to the argc value
        if (argc > 1) {cTestRoutine.restartedFile = true;}
        else {cTestRoutine.restartedFile = false;}
        
        // sets simulation parameters
        cTestRoutine.setTimeStep(5.0e-6);
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
        cTestRoutine.setWallStiffness(k2max);
        cTestRoutine.setParticlePlasticProperties(k1, k2max, kc, phi);
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
        
        // sets teh data points
        cTestRoutine.setLoadingDataPoints(nLoadingDataPoints, heightArrayLoading, pressureArrayLoading, 0.05);
        cTestRoutine.setReloadingDataPoints(nReloadingDataPoints, heightArrayReloading, pressureArrayReloading, 0.05);
        cTestRoutine.setUnloadingDataPoints(nUnloadingDataPoints, heightArrayUnloading, pressureArrayUnloading, 0.05);
        
        // sets additional built-in arguments for the xballs visualization
        cTestRoutine.setXBallsAdditionalArguments("-h 800 -p 10 -o 200 -3dturn 3");
        
        // passes the pointer to the calibrated k array for the computation of the mean and the deviation
        cTestRoutine.kArraySetter(nDataPoints, calibratedKArray);
        
        // reads the .restart file from command line (if any)
        cTestRoutine.readArguments(argc, argv);
        
        // name setting
        if (cTestRoutine.restartedFile)
        {
            std::cout.unsetf(std::ios::floatfield);
            name << "CompressionTest_dataDriven_EXPERIMENTAL_CARRIER_rRatio_" << cylinderToParticleRadiusRatio << "_k1_" << k1 << "_k2max_" << k2max << "_kc_" << kc << "_phi_" << std::fixed << std::setprecision(2) << phi;
            kdataFilename << "CompressionTest_dataDriven_EXPERIMENTAL_CARRIER_rRatio_" << cylinderToParticleRadiusRatio << "_k1_" << k1 << "_k2max_" << k2max << "_kc_" << kc << "_phi_" << std::fixed << std::setprecision(2) << phi << ".kdata";
            cTestRoutine.setName(name.str());
            
            std::cout << "This file has been restarted.\nThe new name of the associated output of this simulation will be:\n" << name.str() << std::endl;
        }
        else
        {
            std::cout.unsetf(std::ios::floatfield);
            name << "CompressionTest_parameterEstimation_RESTARTFILE_CARRIER_rRatio_" << cylinderToParticleRadiusRatio << "_k1_" << k1 << "_k2max_" << k2max << "_kc_" << kc << "_phi_" << std::fixed << std::setprecision(2) << phi;
            kdataFilename  << "CompressionTest_parameterEstimation_RESTARTFILE_CARRIER_rRatio_" << cylinderToParticleRadiusRatio << "_k1_" << k1 << "_k2max_" << k2max << "_kc_" << kc << "_phi_" << std::fixed << std::setprecision(2) << phi << ".kdata";
            cTestRoutine.setName(name.str());
            
            std::cout << "The new name of the associated output of this simulation will be:\n" << name.str() << std::endl;
        }
        
        cTestRoutine.solve();
        
        // computing the mean value and the absolute deviation of the calibrated stiffnesses
        meanCalibratedKArray[runNumber] = computeMean(nLoadingDataPoints, calibratedKArray);
        disArray[runNumber] = computeMeanAbsoluteDeviation(nLoadingDataPoints, calibratedKArray);
        
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
        kdataFilename.str("");
        kdataFilename.clear();
        runNumber++;
    }
    
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







