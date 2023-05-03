/*
    *** COMPRESSION TEST ***
 Particles are placed inside a cylindrical casing and let settle under gravity.
 When the system is at rest a piston is loaded and starts compressing the particles.
 The values of pressure towards the piston, the compression length and ther quantities are printed at each time step.
 
 *** IMPORTANT ***
 The interaction law is purely cohesive with no adhesion.
 If needed uncomment the proper parts and substitute
 LinearPlasticViscoelasticFrictionSpecies -> LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies
 
 LAST UPDATE: 7/07/17
*/

#ifndef MultiStepCompressionTestDataDriven_H
#define MultiStepCompressionTestDataDriven_H

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
 
class CompressionTest_parameterCalibrationRoutine : public Mercury3D
{
private:
    
    void setupInitialConditions() override
    {
        stage = 1;
        setsInserted = 0;
        compressionStep = 0;
        
        // sets particle masses
        std::cout << "Setting particles masses...\n";
        setParticleMassAndVolume();
        std::cout << "DONE\n";
        
        // sets the species
        std::cout << "Setting up the species...\n";
        setSpecies();
        std::cout << "DONE\n";
        
        // computes the number of particle sets needed
        std::cout << "Computing the number of particle sets and the total number of particles needed...\n";
        computeNumberOfSetsAndParticles();
        std::cout << "DONE\n";
        
        // computes the filling region dimensions
        std::cout << "Computing the initialization region dimensions...\n";
        computeInitializationRegionDimensions();
        std::cout << "DONE\n";

        // computes the total volume of particles
        std::cout << "Computing the total volume of particles...\n";
        computeParticleTotalVolume();
        std::cout << "DONE\n";
        
        // computes the particle level and the system max height
        std::cout << "Computing the height of the settled particle bed and the total height of the system...\n";
        computeHeights();
        std::cout << "DONE\n";
        
        // setting the simulation domain
        std::cout << "Setting the simulation domain...\n";
        setBoundaries();
        std::cout << "DONE\n";

        // creating the geometrical components of the casing
        std::cout << "Creating the geometry...\n";
        makeGeometry();
        std::cout << "DONE\n";
        
        // first round of particle insertion
        std::cout << "Inserting the particles...\n";
        makeParticleSet();
        
        stage++;
    }
    
    void actionsOnRestart() override
    {
        stage = 3;
        setsInserted = nSets;
        compressionStep = 0;
        
        // sets particle masses
        std::cout << "Setting particles masses...\n";
        setParticleMassAndVolume();
        std::cout << "DONE\n";
        
        // sets the species
        std::cout << "Setting up the species...\n";
        setSpecies();
        std::cout << "DONE\n";
        
        // computes the number of particle sets needed
        std::cout << "Computing the number of particle sets and the total number of particles needed...\n";
        computeNumberOfSetsAndParticles();
        std::cout << "DONE\n";
    }
    
    void actionsAfterTimeStep() override
    {
        // particle insertion loop and trimming
        if (stage == 2 && !restartedFile)
        {
            // other rounds of particle insertion
            if (setsInserted < nSets && getTime() > 0.05*setsInserted)
            {
                makeParticleSet();
            }
            
            // end of insertion and settling
            if (setsInserted == nSets && getKineticEnergy()/getElasticEnergy() < energyRatioTolerance)
            {
                std::cout << "Particle insertion terminated.\nTrimming the particles...\n";
                particleTrimmer();
                std::cout << "DONE\n";
                stage++;
            }
        }
        
        // resetting times, begin of the compression
        if (stage == 3)
        {            
            std::cout << "Creating the piston...\n";
            makePiston();
            std::cout << "DONE\nCreating output data file...\n";
            makeOutputFile();
            
            // resetting time and tMax
            std::cout << "DONE\nResetting time...\n";
            setTime(0.0);
            std::cout << "DONE\nNew tMax = " << getTimeMax() << "\n";
            
            // sets the first compression step
            std::cout << "Starting the compression stage...\n";
            targetHeight = heightSteps[compressionStep];
            targetPressure = pressureSteps[compressionStep];
            std::cout << "Target height: " << targetHeight << "\n";
            std::cout << "Target pressure: " << targetPressure << "\n";
            resetPistonVelocity();
            heightTargetMet = false;
            pressureTargetMet = false;
            
            timePlaceholder = getTime();
            stage++;
        }
        
        // computes teh pressure and writes the data in the output file
        if (stage >= 4)
        {           
            makeDataAnalysis();
            
            if (fmod(getTime(),cdatOutputTimeInterval) < getTimeStep())
            {
                writeDataToOutptFile();
            }
        }

        // goes on until all the compression and relaxation instances are ultimated
        if (stage == 4)
        {            
            calibrateLoadingBranch();
        }
        
        // writes teh stiffness capibration data to file if data driven
        if (stage == 5)
        {
            writeCalibrationData();
            exportCalibrationData();
        }
    }
    
    void actionAfterSolve()
    {
        // closes the stream to the output file
        outputFile.close();
        
        // deletes the dynamically allocated arrays
        delete[] heightSteps;
        delete[] pressureSteps;
        delete[] calibratedK1;
        delete[] calibratedK2;
        delete[] calibratedKC;
        delete[] calibratedPhi;
    }
    
public:
    
// ----- FUNCTIONS CALLED IN THE MAIN -----
    // sets verbosity on or off
    void setVerbose(bool v)
    {
        verbose = v;
    }
    
    // set the particle density (OVERLOADED)
    void setParticleDensity(double rhoB)
    {
        particleDensityBig = rhoB;
        particleDensitySmall = rhoB;
    }
    void setParticleDensity(double rhoB, double rhoS)
    {
        particleDensityBig = rhoB;
        particleDensitySmall = rhoS;
    }
    
    // set the mean and the dispersity of the particles (OVERLOADED)
    void setParticleRadiusAndDispersity(double mu, double sigma)
    {
        meanRadiusBig = mu;
        meanRadiusSmall = mu;
        dispersityBig = sigma;
        dispersitySmall = sigma;
        
        biModal = false;
    }
    void setParticleRadiusAndDispersity(double muBig, double sigmaBig, double muSmall, double sigmaSmall)
    {
        meanRadiusBig = muBig;
        meanRadiusSmall = muSmall;
        dispersityBig = sigmaBig;
        dispersitySmall = sigmaSmall;
        
        biModal = true;
    }
    
    // set the small-to-big total mass ratio
    void setTotalMassRatio(double ratio)
    {
        smallToBigMassRatio = ratio;
    }
    
    // set wall stiffness
    void setWallStiffness(double kW)
    {
        wallStiffness = kW;
    }
    
    // set particle stiffness (OVERLOADED)
    void setParticleStiffness(double k)
    {
        particleStiffnessBig = k;
        particleStiffnessSmall = k;
    }
    void setParticleStiffness(double kBig, double kSmall)
    {
        particleStiffnessBig = kBig;
        particleStiffnessSmall = kSmall;
    }
    
    // set particle-wall sliding friction coefficients (OVERLOADED)
    void setParticleWallSlidingFrictionCoeff(double bigWallMu)
    {
        bigWallSlidingFrictionCoeff = bigWallMu;
        smallWallSlidingFrictionCoeff = bigWallMu;
    }
    void setParticleWallSlidingFrictionCoeff(double bigWallMu, double smallWallMu)
    {
        bigWallSlidingFrictionCoeff = bigWallMu;
        smallWallSlidingFrictionCoeff = smallWallMu;
    }
    
    // set particle-wall rolling friction coefficients (OVERLOADED)
    void setParticleWallRollingFrictionCoeff(double bigWallMu)
    {
        bigWallRollingFrictionCoeff = bigWallMu;
        smallWallRollingFrictionCoeff = bigWallMu;
    }
    void setParticleWallRollingFrictionCoeff(double bigWallMu, double smallWallMu)
    {
        bigWallRollingFrictionCoeff = bigWallMu;
        smallWallRollingFrictionCoeff = smallWallMu;
    }
    
    // set particle-wall torsion friction coefficients (OVERLOADED)
    void setParticleWallTorsionFrictionCoeff(double bigWallMu)
    {
        bigWallTorsionFrictionCoeff = bigWallMu;
        smallWallTorsionFrictionCoeff = bigWallMu;
    }
    void setParticleWallTorsionFrictionCoeff(double bigWallMu, double smallWallMu)
    {
        bigWallTorsionFrictionCoeff = bigWallMu;
        smallWallTorsionFrictionCoeff = smallWallMu;
    }
    
    // set particle-particle sliding friction coefficients (OVERLOADED)
    void setParticleParticleSlidingFrictionCoeff(double bigBigMu)
    {
        bigBigSlidingFrictionCoeff = bigBigMu;
        smallSmallSlidingFrictionCoeff = bigBigMu;
        bigSmallSlidingFrictionCoeff = bigBigMu;
    }
    void setParticleParticleSlidingFrictionCoeff(double bigBigMu, double smallSmallMu, double bigSmallMu)
    {
        bigBigSlidingFrictionCoeff = bigBigMu;
        smallSmallSlidingFrictionCoeff = smallSmallMu;
        bigSmallSlidingFrictionCoeff = bigSmallMu;
    }
    
    // set particle-particle rolling friction coefficients (OVERLOADED)
    void setParticleParticleRollingFrictionCoeff(double bigBigMu)
    {
        bigBigRollingFrictionCoeff = bigBigMu;
        smallSmallRollingFrictionCoeff = bigBigMu;
        bigSmallRollingFrictionCoeff = bigBigMu;
    }
    void setParticleParticleRollingFrictionCoeff(double bigBigMu, double smallSmallMu, double bigSmallMu)
    {
        bigBigRollingFrictionCoeff = bigBigMu;
        smallSmallRollingFrictionCoeff = smallSmallMu;
        bigSmallRollingFrictionCoeff = bigSmallMu;
    }
    
    // set particle-particle torsion friction coefficients (OVERLOADED)
    void setParticleParticleTorsionFrictionCoeff(double bigBigMu)
    {
        bigBigTorsionFrictionCoeff = bigBigMu;
        smallSmallTorsionFrictionCoeff = bigBigMu;
        bigSmallTorsionFrictionCoeff = bigBigMu;
    }
    void setParticleParticleTorsionFrictionCoeff(double bigBigMu, double smallSmallMu, double bigSmallMu)
    {
        bigBigTorsionFrictionCoeff = bigBigMu;
        smallSmallTorsionFrictionCoeff = smallSmallMu;
        bigSmallTorsionFrictionCoeff = bigSmallMu;
    }
    
    // set particle-wall restitution coefficients (OVERLOADED)
    void setParticleWallRestitutionCoeff(double bigWallE)
    {
        bigWallRestitutionCoeff = bigWallE;
        smallWallRestitutionCoeff = bigWallE;
    }
    void setParticleWallRestitutionCoeff(double bigWallE, double smallWallE)
    {
        bigWallRestitutionCoeff = bigWallE;
        smallWallRestitutionCoeff = smallWallE;
    }
    
    // set particle-particle restitution coefficients (OVERLOADED)
    void setParticleParticleRestitutionCoeff(double bigBigE)
    {
        bigBigRestitutionCoeff = bigBigE;
        smallSmallRestitutionCoeff = bigBigE;
        bigSmallRestitutionCoeff = bigBigE;
    }
    void setParticleParticleRestitutionCoeff(double bigBigE, double smallSmallE, double bigSmallE)
    {
        bigBigRestitutionCoeff = bigBigE;
        smallSmallRestitutionCoeff = smallSmallE;
        bigSmallRestitutionCoeff = bigSmallE;
    }
    
    // set particle plastic properties
    void setParticlePlasticProperties(double k1, double k2max, double kC, double phiC)
    {
        particleStiffnessBig = k1;
        particleStiffnessSmall = k1;
        
        particleMaxUnloadingStiffness = k2max;
        particleCohesiveStiffness = kC;
        particlePlasticityDepth = phiC;
    }
    
//    // set particle adhesion properties
//    void setParticlesAdhesionProperties(double kA, double fAmax)
//    {
//        particleAdhesionStiffness = kA;
//        particleMaxAdhesiveForce = fAmax;
//    }
    
    // set compression cylinder dimensions
    void setCasingDimensions(double radius, double height)
    {
        casingRadius = radius;
        casingHeight = height;
        casingVolume = constants::pi*pow(casingRadius,2.0)*casingHeight;
    }
    
    // set bulk density of the mixture prior to compression
    void setbulkPackingFractionPreCompression(double bulkPF)
    {
        bulkPackingFractionPreCompression = bulkPF;
    }
        
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
       
    // sets the time output between two cdat data prints
    void setCdatOutputTimeInterval(double dt)
    {
        cdatOutputTimeInterval = dt;
    }
    
    // sets the kinetic energy/elastic energy ratio for when the system is static
    void setEnergyRatioTolerance(double eRatio)
    {
        energyRatioTolerance = eRatio;
    }
    
    // allocates the calibration data to a pointer for further usage
    void kArraySetter(int n, double *arrayPointer)
    {
        calibrationStiffnessArrayLength = n;
        pointerToCalibrationRoutineStiffnessArray = arrayPointer;
    }
    
 
// ----- FUNCTIONS CALLED IN THE CLASS -----
    // set the particle masses
    void setParticleMassAndVolume()
    {
        particleVolumeBig = 4.0*constants::pi*pow(meanRadiusBig,3.)/3.;
        particleVolumeSmall = 4.0*constants::pi*pow(meanRadiusSmall,3.)/3.;
        
        particleMassBig = particleDensityBig*particleVolumeBig;
        particleMassSmall = particleDensitySmall*particleVolumeSmall;
        
        if (verbose)
        {
            if (biModal)
            {
                std::cout << "\tParticle mass BIG: " << particleMassBig << "\n";
                std::cout << "\tParticle mass SMALL: " << particleMassSmall << "\n";
            }
            else std::cout << "\tParticle mass: " << particleMassBig << "\n";
        }
    }
    
    // set the particle species
    void setSpecies()
    {
        speciesHandler.clear();
        
        // BIG-BIG
//        specieBig = new LinearViscoelasticFrictionSpecies;
        specieBig = new LinearPlasticViscoelasticFrictionSpecies;
        specieBig -> setDensity(particleDensityBig);
        specieBig -> setStiffnessAndRestitutionCoefficient(particleStiffnessBig, bigBigRestitutionCoeff, particleMassBig);
        // plastic-adhesive part
        specieBig -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
        specieBig -> setCohesionStiffness(particleCohesiveStiffness);
        specieBig -> setPenetrationDepthMax(particlePlasticityDepth);
//        specieBig -> setAdhesionStiffness(particleAdhesionStiffness);
//        specieBig -> setAdhesionForceMax(particleMaxAdhesiveForce);
        
        specieBig -> setSlidingFrictionCoefficient(bigBigSlidingFrictionCoeff);
        specieBig -> setSlidingStiffness(particleStiffnessBig*2.0/7.0);
        specieBig -> setSlidingDissipation(specieBig -> getDissipation()*2.0/7.0);
        
        specieBig -> setRollingFrictionCoefficient(bigBigRollingFrictionCoeff);
        specieBig -> setRollingStiffness(particleStiffnessBig*2.0/7.0);
        specieBig -> setRollingDissipation(specieBig -> getDissipation()*2.0/7.0);
        
        specieBig -> setTorsionFrictionCoefficient(bigBigTorsionFrictionCoeff);
        specieBig -> setTorsionStiffness(particleStiffnessBig*2.0/7.0);
        specieBig -> setTorsionDissipation(specieBig -> getDissipation()*2.0/7.0);
        speciesHandler.addObject(specieBig);
        
        // SMALL-SMALL
//        specieSmall = new LinearViscoelasticFrictionSpecies;
        specieSmall = new LinearPlasticViscoelasticFrictionSpecies;
        specieSmall -> setDensity(particleDensitySmall);
        specieSmall -> setStiffnessAndRestitutionCoefficient(particleStiffnessSmall, smallSmallRestitutionCoeff, particleMassSmall);
        // plastic-adhesive part
        specieSmall -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
        specieSmall -> setCohesionStiffness(particleCohesiveStiffness);
        specieSmall -> setPenetrationDepthMax(particlePlasticityDepth);
//        specieSmall -> setAdhesionStiffness(particleAdhesionStiffness);
//        specieSmall -> setAdhesionForceMax(particleMaxAdhesiveForce);
        
        specieSmall -> setSlidingFrictionCoefficient(smallSmallSlidingFrictionCoeff);
        specieSmall -> setSlidingStiffness(particleStiffnessSmall*2.0/7.0);
        specieSmall -> setSlidingDissipation(specieSmall -> getDissipation()*2.0/7.0);
        
        specieSmall -> setRollingFrictionCoefficient(smallSmallRollingFrictionCoeff);
        specieSmall -> setRollingStiffness(particleStiffnessSmall*2.0/7.0);
        specieSmall -> setRollingDissipation(specieSmall -> getDissipation()*2.0/7.0);
        
        specieSmall -> setTorsionFrictionCoefficient(smallSmallTorsionFrictionCoeff);
        specieSmall -> setTorsionStiffness(particleStiffnessSmall*2.0/7.0);
        specieSmall -> setTorsionDissipation(specieSmall -> getDissipation()*2.0/7.0);
        speciesHandler.addObject(specieSmall);
        
        // WALL-WALL    ( *** DENSITY AND MASS USED ARE THE MEAN OF THE RESPECTIVE PARTICLE BIG AND SMALL ONES, FRICTION SET TO ZERO, RESTITUTION TO 1 *** )
//        specieWall = new LinearViscoelasticFrictionSpecies;
        specieWall = new LinearPlasticViscoelasticFrictionSpecies;
        specieWall -> setDensity(0.5*(particleDensityBig + particleDensitySmall));
        specieWall -> setStiffnessAndRestitutionCoefficient(wallStiffness, 1.0, (particleMassBig + particleMassSmall)/2.0);
        // plastic-adhesive part
        specieWall -> setUnloadingStiffnessMax(0.0);
        specieWall -> setCohesionStiffness(0.0);
        specieWall -> setPenetrationDepthMax(0.0);
//        specieWall -> setAdhesionStiffness(0.0);
//        specieWall -> setAdhesionForceMax(0.0);
        
        specieWall -> setSlidingFrictionCoefficient(0.0);
        specieWall -> setSlidingStiffness(wallStiffness*2.0/7.0);
        specieWall -> setSlidingDissipation(specieWall -> getDissipation()*2.0/7.0);
        
        specieWall -> setRollingFrictionCoefficient(0.0);
        specieWall -> setRollingStiffness(wallStiffness*2.0/7.0);
        specieWall -> setRollingDissipation(specieWall -> getDissipation()*2.0/7.0);
        
        specieWall -> setTorsionFrictionCoefficient(0.0);
        specieWall -> setTorsionStiffness(wallStiffness*2.0/7.0);
        specieWall -> setTorsionDissipation(specieWall -> getDissipation()*2.0/7.0);
        speciesHandler.addObject(specieWall);
        
        // BIG-WALL
        speciesHandler.getMixedObject(specieBig, specieWall) -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessBig + wallStiffness), bigWallRestitutionCoeff, particleMassBig);
        // plastic-adhesive part
        speciesHandler.getMixedObject(specieBig, specieWall) -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
        speciesHandler.getMixedObject(specieBig, specieWall) -> setCohesionStiffness(particleCohesiveStiffness);
        speciesHandler.getMixedObject(specieBig, specieWall) -> setPenetrationDepthMax(particlePlasticityDepth);
//        speciesHandler.getMixedObject(specieBig, specieWall) -> setAdhesionStiffness(particleAdhesionStiffness);
//        speciesHandler.getMixedObject(specieBig, specieWall) -> setAdhesionForceMax(particleMaxAdhesiveForce);
        
        speciesHandler.getMixedObject(specieBig, specieWall) -> setSlidingFrictionCoefficient(bigWallSlidingFrictionCoeff);
        speciesHandler.getMixedObject(specieBig, specieWall) -> setSlidingStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
        speciesHandler.getMixedObject(specieBig, specieWall) -> setSlidingDissipation(speciesHandler.getMixedObject(specieBig, specieWall) -> getDissipation()*2.0/7.0);
        
        speciesHandler.getMixedObject(specieBig, specieWall) -> setRollingFrictionCoefficient(bigWallRollingFrictionCoeff);
        speciesHandler.getMixedObject(specieBig, specieWall) -> setRollingStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
        speciesHandler.getMixedObject(specieBig, specieWall) -> setRollingDissipation(speciesHandler.getMixedObject(specieBig, specieWall) -> getDissipation()*2.0/7.0);
        
        speciesHandler.getMixedObject(specieBig, specieWall) -> setTorsionFrictionCoefficient(bigWallTorsionFrictionCoeff);
        speciesHandler.getMixedObject(specieBig, specieWall) -> setTorsionStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
        speciesHandler.getMixedObject(specieBig, specieWall) -> setTorsionDissipation(speciesHandler.getMixedObject(specieBig, specieWall) -> getDissipation()*2.0/7.0);
        
        // SMALL-WALL
        speciesHandler.getMixedObject(specieSmall, specieWall) -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessSmall + wallStiffness), smallWallRestitutionCoeff, particleMassSmall);
        // plastic-adhesive part
        speciesHandler.getMixedObject(specieSmall, specieWall) -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
        speciesHandler.getMixedObject(specieSmall, specieWall) -> setCohesionStiffness(particleCohesiveStiffness);
        speciesHandler.getMixedObject(specieSmall, specieWall) -> setPenetrationDepthMax(particlePlasticityDepth);
//        speciesHandler.getMixedObject(specieSmall, specieWall) -> setAdhesionStiffness(particleAdhesionStiffness);
//        speciesHandler.getMixedObject(specieSmall, specieWall) -> setAdhesionForceMax(particleMaxAdhesiveForce);
        
        speciesHandler.getMixedObject(specieSmall, specieWall) -> setSlidingFrictionCoefficient(smallWallSlidingFrictionCoeff);
        speciesHandler.getMixedObject(specieSmall, specieWall) -> setSlidingStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
        speciesHandler.getMixedObject(specieSmall, specieWall) -> setSlidingDissipation(speciesHandler.getMixedObject(specieSmall, specieWall) -> getDissipation()*2.0/7.0);
        
        speciesHandler.getMixedObject(specieSmall, specieWall) -> setRollingFrictionCoefficient(smallWallRollingFrictionCoeff);
        speciesHandler.getMixedObject(specieSmall, specieWall) -> setRollingStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
        speciesHandler.getMixedObject(specieSmall, specieWall) -> setRollingDissipation(speciesHandler.getMixedObject(specieSmall, specieWall) -> getDissipation()*2.0/7.0);
        
        speciesHandler.getMixedObject(specieSmall, specieWall) -> setTorsionFrictionCoefficient(smallWallTorsionFrictionCoeff);
        speciesHandler.getMixedObject(specieSmall, specieWall) -> setTorsionStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
        speciesHandler.getMixedObject(specieSmall, specieWall) -> setTorsionDissipation(speciesHandler.getMixedObject(specieSmall, specieWall) -> getDissipation()*2.0/7.0);
        
        // BIG-SMALL
        speciesHandler.getMixedObject(specieBig, specieSmall) -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessBig + particleStiffnessSmall), bigSmallRestitutionCoeff, 0.5*(particleMassBig + particleMassSmall));
        // plastic-adhesive part
        speciesHandler.getMixedObject(specieBig, specieSmall) -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
        speciesHandler.getMixedObject(specieBig, specieSmall) -> setCohesionStiffness(particleCohesiveStiffness);
        speciesHandler.getMixedObject(specieBig, specieSmall) -> setPenetrationDepthMax(particlePlasticityDepth);
//        speciesHandler.getMixedObject(specieBig, specieSmall) -> setAdhesionStiffness(particleAdhesionStiffness);
//        speciesHandler.getMixedObject(specieBig, specieSmall) -> setAdhesionForceMax(particleMaxAdhesiveForce);
        
        speciesHandler.getMixedObject(specieBig, specieSmall) -> setSlidingFrictionCoefficient(bigSmallSlidingFrictionCoeff);
        speciesHandler.getMixedObject(specieBig, specieSmall) -> setSlidingStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
        speciesHandler.getMixedObject(specieBig, specieSmall) -> setSlidingDissipation(speciesHandler.getMixedObject(specieBig, specieSmall) -> getDissipation()*2.0/7.0);
        
        speciesHandler.getMixedObject(specieBig, specieSmall) -> setRollingFrictionCoefficient(bigSmallRollingFrictionCoeff);
        speciesHandler.getMixedObject(specieBig, specieSmall) -> setRollingStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
        speciesHandler.getMixedObject(specieBig, specieSmall) -> setRollingDissipation(speciesHandler.getMixedObject(specieBig, specieSmall) -> getDissipation()*2.0/7.0);
        
        speciesHandler.getMixedObject(specieBig, specieSmall) -> setTorsionFrictionCoefficient(bigSmallTorsionFrictionCoeff);
        speciesHandler.getMixedObject(specieBig, specieSmall) -> setTorsionStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
        speciesHandler.getMixedObject(specieBig, specieSmall) -> setTorsionDissipation(speciesHandler.getMixedObject(specieBig, specieSmall) -> getDissipation()*2.0/7.0);
        
        if (verbose)
        {
            if (biModal)
            {
                std::cout << "\tBIG-BIG stiffness and dissipation coefficients: " << specieBig -> getLoadingStiffness() << " " << specieBig -> getDissipation() << "\n";
                std::cout << "\tBIG-BIG friction coefficients: " << bigBigSlidingFrictionCoeff << " " << bigBigRollingFrictionCoeff << " " << bigBigTorsionFrictionCoeff << "\n";
                std::cout << "\tBIG-BIG tangential stiffnesses: " << specieBig -> getSlidingStiffness() << " " << specieBig -> getRollingStiffness() << " " << specieBig -> getTorsionStiffness() << "\n";
                std::cout << "\tBIG-BIG tangential dissipation coefficients: " << specieBig -> getSlidingDissipation() << " " << specieBig -> getRollingDissipation() << " " << specieBig -> getTorsionDissipation() << "\n";
                std::cout << "\tBIG-BIG collision time: " << std::setprecision(4) << specieBig -> getCollisionTime(particleMassBig) << "\n\n";
                
                std::cout << "\tSMALL-SMALL stiffness and dissipation coefficients: " << specieSmall -> getLoadingStiffness() << " " << specieSmall -> getDissipation() << "\n";
                std::cout << "\tSMALL-SMALL friction coefficients: " << smallSmallSlidingFrictionCoeff << " " << smallSmallRollingFrictionCoeff << " " << smallSmallTorsionFrictionCoeff << "\n";
                std::cout << "\tSMALL-SMALL tangential stiffnesses: " << specieSmall -> getSlidingStiffness() << " " << specieSmall -> getRollingStiffness() << " " << specieSmall -> getTorsionStiffness() << "\n";
                std::cout << "\tSMALL-SMALL tangential dissipation coefficients: " << specieSmall -> getSlidingDissipation() << " " << specieSmall -> getRollingDissipation() << " " << specieSmall -> getTorsionDissipation() << "\n";
                std::cout << "\tSMALL-SMALL collision time: " << std::setprecision(4) << specieSmall -> getCollisionTime(particleMassSmall) << "\n\n";
                
                std::cout << "\tBIG-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(specieBig, specieWall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getDissipation() << "\n";
                std::cout << "\tBIG-WALL friction coefficients: " << bigWallSlidingFrictionCoeff << " " << bigWallRollingFrictionCoeff << " " << bigWallTorsionFrictionCoeff << "\n";
                std::cout << "\tBIG-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(specieBig, specieWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getTorsionStiffness() << "\n";
                std::cout << "\tBIG-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(specieBig, specieWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getTorsionDissipation() << "\n";
                std::cout << "\tBIG-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(specieBig, specieWall) -> getCollisionTime(particleMassBig) << "\n\n";
                
                std::cout << "\tSMALL-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getDissipation() << "\n";
                std::cout << "\tSMALL-WALL friction coefficients: " << smallWallSlidingFrictionCoeff << " " << smallWallRollingFrictionCoeff << " " << smallWallTorsionFrictionCoeff << "\n";
                std::cout << "\tSMALL-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getTorsionStiffness() << "\n";
                std::cout << "\tSMALL-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getTorsionDissipation() << "\n";
                std::cout << "\tSMALL-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(specieSmall, specieWall) -> getCollisionTime(particleMassSmall) << "\n\n";
                
                std::cout << "\tBIG-SMALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(specieBig, specieSmall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieSmall) -> getDissipation() << "\n";
                std::cout << "\tBIG-SMALL friction coefficients: " << bigSmallSlidingFrictionCoeff << " " << bigSmallRollingFrictionCoeff << " " << bigSmallTorsionFrictionCoeff << "\n";
                std::cout << "\tBIG-SMALL tangential stiffnesses: " << speciesHandler.getMixedObject(specieBig, specieSmall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieSmall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieSmall) -> getTorsionStiffness() << "\n";
                std::cout << "\tBIG-SMALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(specieBig, specieSmall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(specieBig, specieSmall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(specieBig, specieSmall) -> getTorsionDissipation() << "\n";
                std::cout << "\tBIG-SMALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(specieBig, specieSmall) -> getCollisionTime(0.5*(particleMassBig + particleMassSmall)) << "\n\n";
                
                std::cout << "\tBIG-BIG collision time / TIME STEP: " << std::setprecision(4) << (specieBig -> getCollisionTime(particleMassBig))/getTimeStep() << "\n";
                std::cout << "\tSMALL-SMALL collision time / TIME STEP: " << std::setprecision(4) << (specieSmall -> getCollisionTime(particleMassSmall))/getTimeStep() << "\n";
                std::cout << "\tBIG-WALL collision time / TIME STEP: " << std::setprecision(4) << (speciesHandler.getMixedObject(specieBig, specieWall) -> getCollisionTime(particleMassBig))/getTimeStep() << "\n";
                std::cout << "\tSMALL-WALL collision time / TIME STEP: " << std::setprecision(4) << (speciesHandler.getMixedObject(specieSmall, specieWall) -> getCollisionTime(particleMassSmall))/getTimeStep() << "\n";
                std::cout << "\tBIG-SMALL collision time / TIME STEP: " << std::setprecision(4) << (speciesHandler.getMixedObject(specieBig, specieSmall) -> getCollisionTime(0.5*(particleMassBig + particleMassSmall)))/getTimeStep() << "\n\n";
            }
            else
            {
                std::cout << "\tPARTICLE stiffness and dissipation coefficients: " << specieBig -> getLoadingStiffness() << " " << specieBig -> getDissipation() << "\n";
                std::cout << "\tPARTICLE friction coefficients: " << bigBigSlidingFrictionCoeff << " " << bigBigRollingFrictionCoeff << " " << bigBigTorsionFrictionCoeff << "\n";
                std::cout << "\tPARTICLE tangential stiffnesses: " << specieBig -> getSlidingStiffness() << " " << specieBig -> getRollingStiffness() << " " << specieBig -> getTorsionStiffness() << "\n";
                std::cout << "\tPARTICLE tangential dissipation coefficients: " << specieBig -> getSlidingDissipation() << " " << specieBig -> getRollingDissipation() << " " << specieBig -> getTorsionDissipation() << "\n";
                std::cout << "\tPARTICLE collision time: " << std::setprecision(4) << specieBig -> getCollisionTime(particleMassBig) << "\n\n";
                
                std::cout << "\tPARTICLE-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(specieBig, specieWall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getDissipation() << "\n";
                std::cout << "\tPARTICLE-WALL friction coefficients: " << bigWallSlidingFrictionCoeff << " " << bigWallRollingFrictionCoeff << " " << bigWallTorsionFrictionCoeff << "\n";
                std::cout << "\tPARTICLE-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(specieBig, specieWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getTorsionStiffness() << "\n";
                std::cout << "\tPARTICLE-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(specieBig, specieWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getTorsionDissipation() << "\n";
                std::cout << "\tPARTICLE-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(specieBig, specieWall) -> getCollisionTime(particleMassBig) << "\n\n";
                
                std::cout << "\tPARTICLE collision time / TIME STEP: " << std::setprecision(4) << (specieBig -> getCollisionTime(particleMassBig))/getTimeStep() << "\n";
                std::cout << "\tPARTICLE-WALL collision time / TIME STEP: " << std::setprecision(4) << (speciesHandler.getMixedObject(specieBig, specieWall) -> getCollisionTime(particleMassBig))/getTimeStep() << "\n\n";
            }
        }
    }
    
    // compute the number of particles and the number of sets
    void computeNumberOfSetsAndParticles()
    {
        nSmallParticles = (int)((bulkPackingFractionPreCompression*casingVolume/particleVolumeSmall)*(smallToBigMassRatio*particleDensityBig/(particleDensitySmall + smallToBigMassRatio*particleDensityBig)));
        nBigParticles = (int)((bulkPackingFractionPreCompression*casingVolume - nSmallParticles*particleVolumeSmall)/particleVolumeBig);
 
        // the max number of sets is set to 20
        nSets = 20;
        
        // takes the smaller between nBigParticles and nSmallParticles
        // starts from 20 and checks if nXPerSet > 1
        // if not, nSets is decreased to 19 and the check is repeated, and so on
        if (nSmallParticles > nBigParticles) {while (nBigParticles/nSets < 1.0) {nSets--;};}
        else {while (nSmallParticles/nSets < 1.0) {nSets--;};}
        
        nSmallPerSet = (int)(nSmallParticles/nSets);
        nBigPerSet = (int)(nBigParticles/nSets);
        nParticlesPerSet = nSmallPerSet + nBigPerSet;
        
        if (verbose)
        {
            std::cout << "\tNumber of BIG particles needed: " << nBigParticles << "\n";
            std::cout << "\tNumber of SMALL particles needed: " << nSmallParticles << "\n";
            std::cout << "\tNumber of particle sets needed: " << nSets << "\n";
            std::cout << "\tNumber of BIG particles loaded per set: " << nBigPerSet << "\n";
            std::cout << "\tNumber of SMALL particles loaded per set: " << nSmallPerSet << "\n";
            std::cout << "\tNumber of TOTAL particles loaded per set: " << nParticlesPerSet << "\n";
            std::cout << "\tTotal number of BIG particles loaded: " << nBigPerSet*nSets << "\n";
            std::cout << "\tTotal number of SMALL particles loaded: " << nSmallPerSet*nSets << "\n";
            std::cout << "\tSanity check: NB/NS = " << nBigParticles/nSmallParticles << ", nB/nS = " << nBigPerSet/nSmallPerSet << ", (NB/NS)/(nB/nS) = " << (nBigParticles/nSmallParticles)/(nBigPerSet/nSmallPerSet) << "\n";
        }
    }
    
    // computes the filling region dimensions based on the mean particle sizes and the number of particles per set
    void computeInitializationRegionDimensions()
    {
        // compute the mean volume of each set of particles
        volumeOfParticlesPerSet = nBigPerSet*particleVolumeBig + nSmallPerSet*particleVolumeSmall;
        
        // computes the height of the initialization region assuming a packing fraction of 0.1
        // vol_set = 0.1*vol_cyl = 0.1*pi*r^2*h -> h = vol_set/(0.1*pi*r^2)
        initRegionHeight = 10.0*volumeOfParticlesPerSet/(constants::pi*pow(casingRadius,2.0));
        
        if (verbose)
        {
            std::cout << "\tVolume of particle set: " << volumeOfParticlesPerSet << "\n";
            std::cout << "\tInitialization region height: " << initRegionHeight << "\n";
        }
    }
    
    // computes the total volume of all the loaded particles
    void computeParticleTotalVolume()
    {
        particleTotalVolume = nSets*volumeOfParticlesPerSet;
        
        if (verbose)
        {
            std::cout << "\tTotal volume of particles: " << particleTotalVolume << "\n";
            std::cout << "\tSanity check (nSets*volumePerSet)/(bulkPackingFraction*casingVolume): " << particleTotalVolume/(bulkPackingFractionPreCompression*casingVolume) << "\n";
        }
    }

    // computes the particle bed height and the filling region positioning
    void computeHeights()
    {
        double settledBedHeight;
        settledBedHeight = particleTotalVolume/(bulkPackingFractionPreCompression*constants::pi*pow(casingRadius,2.0));
        
        // the total height of the casing is set to be bed_height + 1.5*init_region_height
        totalHeight = settledBedHeight + 1.5*initRegionHeight;
        
//        // UNCOMMENT FOR SINGLE SETTLING STEP
//        totalHeight = 0.5*casingHeight + 2.1*meanRadiusBig*(1.0 + dispersityBig)*numberOfLevels;
        
        if (verbose)
        {
            std::cout << "\tThe settled powder bed height is assumed to be: " << settledBedHeight << "\n";
            std::cout << "\tThe total height of the system is set to: " << totalHeight << "\n";
            std::cout << "\tSanity check: settledBedHeight/casingHeight = " << settledBedHeight/casingHeight << "\n";
        }
    }

    // sets the simulation domain and the walls
    void setBoundaries()
    {
        // for the X-Y values 1.1*casingRadius is used, for the Z boundaries the extrema + 0.1*casingRadius is used
        setXMin(-1.1*casingRadius);
        setYMin(-1.1*casingRadius);
        setZMin(-0.1*casingRadius);
        
        setXMax(1.1*casingRadius);
        setYMax(1.1*casingRadius);
        setZMax(totalHeight + 0.1*casingRadius);
        
        if (verbose)
        {
            std::cout << "\tSimulation domain minimum point: (" << getXMin() << " , " << getYMin() << " , " << getZMin() << ")\n";
            std::cout << "\tSimulation domain maximum point: (" << getXMax() << " , " << getYMax() << " , " << getZMax() << ")\n";
        }
    }

    // makes the geometric components
    void makeGeometry()
    {
        wallHandler.clear();
        
        // the basis of the compaction chamber
        basis.setSpecies(specieWall);
        basis.set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,0.0));
        wallHandler.copyAndAddObject(basis);
        
        // the roof of the compaction chamber
        roof.setSpecies(specieWall);
        roof.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,totalHeight));
        wallHandler.copyAndAddObject(roof);
        
        // the external casing
        AxisymmetricIntersectionOfWalls(casing);
        casing.setSpecies(specieWall);
        casing.setPosition(Vec3D(0.0,0.0,0.0));
        casing.setOrientation(Vec3D(0.0,0.0,1.0));
        casing.addObject(Vec3D(1.0,0.0,0.0),Vec3D(casingRadius,0.0,0.0));
        wallHandler.copyAndAddObject(casing);
    }
    
    // makes one particle set
    void makeParticleSet()
    {
        double x, y, z;
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        
        for (int i=0; i<nParticlesPerSet; i++)
        {
            if (i<nBigPerSet)
            {
                if (!setsInserted && !i) {p0.setRadius(meanRadiusBig*(1.0 + dispersityBig));}
                else {p0.setRadius(meanRadiusBig*(1.0 + dispersityBig*random.getRandomNumber(-1.0,1.0)));}
                p0.setSpecies(specieBig);
            }
            else
            {
                if (!setsInserted && !i) {p0.setRadius(meanRadiusSmall*(1.0 - dispersitySmall));}
                else {p0.setRadius(meanRadiusSmall*(1.0 + dispersitySmall*random.getRandomNumber(-1.0,1.0)));}
                p0.setSpecies(specieSmall);
            }
            
            do {
                x = (casingRadius - 1.1*p0.getRadius())*random.getRandomNumber(-1.0,1.0);
                y = (casingRadius - 1.1*p0.getRadius())*random.getRandomNumber(-1.0,1.0);
            } while (pow(x,2.0) + pow(y,2.0) > pow(casingRadius - 1.1*p0.getRadius(),2.0));
            z = totalHeight - p0.getRadius() - initRegionHeight*random.getRandomNumber(0.0,1.0);
            
            p0.setPosition(Vec3D(x, y, z));
            particleHandler.copyAndAddObject(p0);
        }
        
        setsInserted++;
        std::cout << "Inserted particle set n " << setsInserted << "\n";
        
        if (verbose)
        {
            std::cout << "\tTotal number of particles: " << particleHandler.getNumberOfObjects() << "\n";
            std::cout << "\tSanity check: setsInserted*nParticlesPerSet/numberOfParticlesInTheSystem = " << setsInserted*nParticlesPerSet/particleHandler.getNumberOfObjects() << "\n";
        }
    }
    
    // removes the particles lying above the casing height
    void particleTrimmer()
    {
        for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            if (particleHandler.getObject(i) -> getPosition().Z + particleHandler.getObject(i) -> getRadius() > casingHeight) particleHandler.removeObject(i);
        }
        
        if (verbose)
        {
            std::cout << "Trimming of the particle done. New max particle height: " << particleHandler.getHighestPositionComponentParticle(2) -> getPosition().Z << "\n";
        }
    }
    
    // creates the piston
    void makePiston()
    {
        // the piston initial height is created at 2*maximumParticleRadiusPossible over the highest particle
        pistonHeight = particleHandler.getHighestPositionComponentParticle(2) -> getPosition().Z + 2.0*meanRadiusBig*(1.0 + dispersityBig);
        
        // the compression piston
        piston.setSpecies(specieWall);
        piston.set(pistonHeight, casingRadius, 0.0, 0.0);
        pistonPointer = wallHandler.copyAndAddObject(piston);
    }
    
    // sets the piston velocity based on the smallest particle size
    void resetPistonVelocity()
    {
        // takes the smallest difference in height steps and uses that as a reference for the piston velocity, such that vPistonMax = dHmin/(100*timeStep)
        double dHmin;
        
        // the dHmin is calculated after every compression step to tune the velocity accordingly
        if (compressionStep == 0)
        {
            dHmin = pistonHeight - heightSteps[0];
        }
        else
        {
            dHmin = fabs(heightSteps[compressionStep] - heightSteps[compressionStep-1]);
        }
        
        // takes the smallest velocity between dHmin/(100*timeStep) and 0.0005*smallestParticleRatio/timeStep
        (dHmin/(100.0*getTimeStep()) < 0.0005*meanRadiusSmall*(1.0 - dispersitySmall)/getTimeStep()) ? (pistonVelocity = -dHmin/(100.0*getTimeStep())) : (pistonVelocity = -0.0005*meanRadiusSmall*(1.0 - dispersitySmall)/getTimeStep());
        
        pistonPointer -> setVelocity(pistonVelocity);
        
        if (verbose)
        {
            std::cout << "\tRESETTING PISTON VELOCITY\n";
            std::cout << "\tPiston velocity: " << pistonVelocity << "\n";
            std::cout << "\tDisplacement per time step: " << pistonVelocity*getTimeStep() << "\n";
        }
    }
    
    // performs the data analysis needed for teh problem
    void makeDataAnalysis()
    {
        // evaluation of mean coordination number
        meanCoordinationNumber = 0.0;
        
        for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
        }
        
        meanCoordinationNumber /= particleHandler.getNumberOfObjects();
        
        // evaluation of particle overlaps and pressure
        pistonPressure = 0.0;
        basePressure = 0.0;
        meanTotalRelativeOverlap = 0.0;
        maxTotalRelativeOverlap = 0.0;
        meanPistonRelativeOverlap = 0.0;
        maxPistonRelativeOverlap = 0.0;
        meanBaseRelativeOverlap = 0.0;
        maxBaseRelativeOverlap = 0.0;
        
        int totalInteractionCounter = 0;
        int pistonInteractionCounter = 0;
        int baseInteractionCounter = 0;
        
        for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
        {
            // mean and maximum total overlap computation
            meanTotalRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxTotalRelativeOverlap) maxTotalRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            
            // piston pressure and mean and maximum particle-piston overlap computation
            if ((*i) -> getI() -> getIndex() == pistonPointer -> getIndex())
            {
                pistonPressure -= ((*i) -> getForce()).Z;
                
                meanPistonRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxPistonRelativeOverlap) maxPistonRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                
                pistonInteractionCounter++;
            }
            
            // base pressure and mean and maximum particle-base overlap computation
            if ((*i) -> getI() -> getIndex() == basis.getIndex())
            {
                basePressure += ((*i) -> getForce()).Z;
                
                meanBaseRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxBaseRelativeOverlap) maxBaseRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                
                baseInteractionCounter++;
            }
            
            totalInteractionCounter++;
        }
        
        pistonPressure /= constants::pi*pow(casingRadius,2.0);
        basePressure /= constants::pi*pow(casingRadius,2.0);
        meanTotalRelativeOverlap /= totalInteractionCounter;
        meanPistonRelativeOverlap /= pistonInteractionCounter;
        meanBaseRelativeOverlap /= baseInteractionCounter;
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
            
            if (compressionStep < loadingCompressionInstances)
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

//     // unloads the piston
//     void unloadPiston()
//     {
//         // sets the velocity such that it takes 1.0 second to get to the maximum height
//         pistonVelocity = (totalHeight - pistonPointer -> getHeight())/1.0;
//         pistonPointer -> setVelocity(pistonVelocity);
//         
//         if (pistonPointer -> getHeight() < totalHeight) pistonPointer -> movePiston(getTimeStep());
//     }
    
    // creates the data output file and writes the first row
    void makeOutputFile()
    {
        std::ostringstream cdatName;
        std::cout.unsetf(std::ios::floatfield);
        cdatName << getName() << ".cdat";
        
        outputFile.open(cdatName.str(), std::ios::out);
        outputFile << "time \t Eel \t Ekin/Eel \t h \t v \t P \t Pbase \t Pmax \t P/Pmax \t |(P - Pmax)/Pmax| \t dTotMean \t dTotMax \t dPistonMean \t dPistonMax \t dBaseMean \t dBaseMax \t k1 \t k2 \t kC \t phi" << std::endl;
    }
    
    // writes the compression data to the output file
    void writeDataToOutptFile()
    {
        outputFile <<
        getTime() << "   " <<
        getElasticEnergy() << "   " <<
        getKineticEnergy()/getElasticEnergy() << "   " <<
        pistonPointer -> getHeight() << "   " <<
        pistonPointer -> getVelocity() << "   " <<
        pistonPressure << "   " <<
        basePressure << "   " <<
        targetPressure << "   " <<
        (pistonPressure + 0.01)/(targetPressure + 0.01) << "   " <<
        fabs((pistonPressure - targetPressure)/(targetPressure + 0.01)) << "   " <<
        meanTotalRelativeOverlap << "   " <<
        maxTotalRelativeOverlap << "   " <<
        meanPistonRelativeOverlap << "   " <<
        maxPistonRelativeOverlap << "   " <<
        meanBaseRelativeOverlap << "   " <<
        maxBaseRelativeOverlap << "   " <<
        particleStiffnessBig << "   " <<
        particleMaxUnloadingStiffness << "   " <<
        particleCohesiveStiffness << "   " <<
        particlePlasticityDepth << "   " <<
        std::endl;
    }
    
    // adjusts the particle stiffness according to the pressure to be met and refreshes the species parameters involved
    void adjustK1()
    {
        // the stiffness is updated every 50 time steps
        if (fmod(getTime() - timePlaceholder, 50.0*getTimeStep()) < getTimeStep())
        {
            // if the pressure is too low increases the stiffness and vice-versa
            if (pistonPressure < targetPressure)
            {
                particleStiffnessBig *= 1.001;
                particleStiffnessSmall *= 1.001;
            }
            else
            {
                particleStiffnessBig *= 0.999;
                particleStiffnessSmall *= 0.999;
            }
            if (particleStiffnessBig > particleMaxUnloadingStiffness) 
            {
                std::cout << "\n\nFATAL ERROR: K1 > K2MAX.";
                exit(0);
            }
            
            specieBig -> setStiffnessAndRestitutionCoefficient(particleStiffnessBig, bigBigRestitutionCoeff, particleMassBig);
            specieBig -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
            specieBig -> setSlidingStiffness(particleStiffnessBig*2.0/7.0);
            specieBig -> setSlidingDissipation(specieBig -> getDissipation()*2.0/7.0);
            specieBig -> setRollingStiffness(particleStiffnessBig*2.0/7.0);
            specieBig -> setRollingDissipation(specieBig -> getDissipation()*2.0/7.0);
            specieBig -> setTorsionStiffness(particleStiffnessBig*2.0/7.0);
            specieBig -> setTorsionDissipation(specieBig -> getDissipation()*2.0/7.0);
            
            specieSmall -> setStiffnessAndRestitutionCoefficient(particleStiffnessSmall, smallSmallRestitutionCoeff, particleMassSmall);
            specieSmall -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
            specieSmall -> setSlidingStiffness(particleStiffnessSmall*2.0/7.0);
            specieSmall -> setSlidingDissipation(specieSmall -> getDissipation()*2.0/7.0);
            specieSmall -> setRollingStiffness(particleStiffnessSmall*2.0/7.0);
            specieSmall -> setRollingDissipation(specieSmall -> getDissipation()*2.0/7.0);
            specieSmall -> setTorsionStiffness(particleStiffnessSmall*2.0/7.0);
            specieSmall -> setTorsionDissipation(specieSmall -> getDissipation()*2.0/7.0);
            
            speciesHandler.getMixedObject(specieBig, specieWall) -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessBig + wallStiffness), bigWallRestitutionCoeff, particleMassBig);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setSlidingStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setSlidingDissipation(speciesHandler.getMixedObject(specieBig, specieWall) -> getDissipation()*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setRollingStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setRollingDissipation(speciesHandler.getMixedObject(specieBig, specieWall) -> getDissipation()*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setTorsionStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setTorsionDissipation(speciesHandler.getMixedObject(specieBig, specieWall) -> getDissipation()*2.0/7.0);
            
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessSmall + wallStiffness), smallWallRestitutionCoeff, particleMassSmall);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setSlidingStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setSlidingDissipation(speciesHandler.getMixedObject(specieSmall, specieWall) -> getDissipation()*2.0/7.0);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setRollingStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setRollingDissipation(speciesHandler.getMixedObject(specieSmall, specieWall) -> getDissipation()*2.0/7.0);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setTorsionStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setTorsionDissipation(speciesHandler.getMixedObject(specieSmall, specieWall) -> getDissipation()*2.0/7.0);
            
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessBig + particleStiffnessSmall), bigSmallRestitutionCoeff, 0.5*(particleMassBig + particleMassSmall));
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setSlidingStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setSlidingDissipation(speciesHandler.getMixedObject(specieBig, specieSmall) -> getDissipation()*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setRollingStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setRollingDissipation(speciesHandler.getMixedObject(specieBig, specieSmall) -> getDissipation()*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setTorsionStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setTorsionDissipation(speciesHandler.getMixedObject(specieBig, specieSmall) -> getDissipation()*2.0/7.0);
        }
    }
    
    // adjusts the particle stiffness according to the pressure to be met and refreshes the species parameters involved
    void adjustKC()
    {
        // the stiffness is updated every 50 time steps
        if (fmod(getTime() - timePlaceholder, 50.0*getTimeStep()) < getTimeStep())
        {
            // if the pressure is too low decrease cohesion and vice-versa
            if (pistonPressure < targetPressure)
            {
                particleCohesiveStiffness *= 0.999;
            }
            else
            {
                particleCohesiveStiffness *= 1.001;
            }
            
            if (particleCohesiveStiffness < 0) 
            {
                std::cout << "\n\nFATAL ERROR: COHESIVE STIFFNESS < 0.";
                exit(0);
            }
            
            specieBig -> setCohesionStiffness(particleCohesiveStiffness);
            specieSmall -> setCohesionStiffness(particleCohesiveStiffness);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setCohesionStiffness(particleCohesiveStiffness);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setCohesionStiffness(particleCohesiveStiffness);
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setCohesionStiffness(particleCohesiveStiffness);
        }

    }
    
    // writes the calibration data to file
    void writeCalibrationData()
    {
        std::ostringstream kdataName;
        std::cout.unsetf(std::ios::floatfield);
        kdataName << getName() << ".kdata";
        
        calibrationDataFile.open(kdataName.str(), std::ios::out);
        calibrationDataFile << "nLoadingSteps \t nReloadingSteps \t nUnloadingSteps" << std::endl;
        calibrationDataFile << loadingCompressionInstances << "\t" << reloadingCompressionInstances << "\t" << unloadingCompressionInstances << std::endl;
        calibrationDataFile << "height \t pressure \t k1 \t k2 \t kC \t phi" << std::endl;
        for (int i = 0; i < loadingCompressionInstances; i++)
        {
            calibrationDataFile << heightSteps[i] << "\t" << pressureSteps[i] << "\t" << calibratedK1[i] << "\t" << calibratedK2[i] << "\t" << calibratedKC[i] << "\t" << calibratedPhi[i] << std::endl;
        }
        for (int i = 0; i < reloadingCompressionInstances; i++)
        {
            calibrationDataFile << heightSteps[i] << "\t" << pressureSteps[i] << "\t" << calibratedK1[i] << "\t" << calibratedK2[i] << "\t" << calibratedKC[i] << "\t" << calibratedPhi[i] << std::endl;
        }
        for (int i = 0; i < unloadingCompressionInstances; i++)
        {
            calibrationDataFile << heightSteps[i] << "\t" << pressureSteps[i] << "\t" << calibratedK1[i] << "\t" << calibratedK2[i] << "\t" << calibratedKC[i] << "\t" << calibratedPhi[i] << std::endl;
        }

        calibrationDataFile.close();
    }
    
    // allocates the calibration data to an external pointer for further usage
    void exportCalibrationData()
    {
        for (int i = 0; i < calibrationStiffnessArrayLength; i++)
        {
            pointerToCalibrationRoutineStiffnessArray[i] = calibratedK1[i];
        }
    }
    
 
//  ----- GLOBAL FUNCTIONS -----
    void printTime() const override
    {
        if (verbose)
        {
            if (stage < 4)
            {
                std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() << ", Eratio = " << std::setprecision(6) << std::left << std::setw(10) << getKineticEnergy()/getElasticEnergy() << std::endl;
                std::cout.flush();
            }
            else
            {
                std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() << ", Elastic Energy = " << std::setprecision(6) << std::left << std::setw(10) << getElasticEnergy() << ", Eratio = " << getKineticEnergy()/getElasticEnergy() << ", h = " << pistonPointer -> getHeight() << ", v = " << pistonPointer -> getVelocity() << "\nP = " << std::setprecision(6) << std::left << std::setw(10) << pistonPressure << ", Pbase = " << basePressure << " , Pmax = " << targetPressure << ", P/Pmax = " << (pistonPressure + 0.01)/(targetPressure + 0.01) << ", |(P - Pmax)/Pmax| = " << fabs((pistonPressure - targetPressure)/(targetPressure + 0.01)) << "\ncN = " << std::setprecision(6) << std::left << std::setw(10) << meanCoordinationNumber << ", meanTotOverlap = " << meanTotalRelativeOverlap << ", maxTotOverlap = " << maxTotalRelativeOverlap << ", meanPistonOverlap = " << meanPistonRelativeOverlap << ", maxPistonOverlap = " << maxPistonRelativeOverlap << ", meanBaseOverlap = " << meanBaseRelativeOverlap << ", maxBaseOverlap = " << maxBaseRelativeOverlap << "\nk1 = " << std::setprecision(6) << std::left << std::setw(10) << specieBig -> getLoadingStiffness() << ", k2 = " << specieBig -> getUnloadingStiffnessMax() << ", kC = " << specieBig -> getCohesionStiffness() << ", phi = " << specieBig -> getPenetrationDepthMax() << std::endl << std::endl;
                std::cout.flush();
            }
        }
    }

    
// ----- VARIABLES -----
    // particle intrinsic properties
    double particleDensityBig, particleDensitySmall;
    double meanRadiusBig, meanRadiusSmall;
    double dispersityBig, dispersitySmall;
    double particleVolumeBig, particleVolumeSmall;
    double particleMassBig, particleMassSmall;
    double particleTotalVolume;
    
    // particle interaction properties
    double wallStiffness;
    double particleStiffnessBig, particleStiffnessSmall;
    double bigWallRestitutionCoeff, smallWallRestitutionCoeff;
    double bigBigRestitutionCoeff, smallSmallRestitutionCoeff, bigSmallRestitutionCoeff;
    // particle-wall friction coefficients
    double bigWallSlidingFrictionCoeff, smallWallSlidingFrictionCoeff;
    double bigWallRollingFrictionCoeff, smallWallRollingFrictionCoeff;
    double bigWallTorsionFrictionCoeff, smallWallTorsionFrictionCoeff;
    // particle-particle friction coefficients
    double bigBigSlidingFrictionCoeff, smallSmallSlidingFrictionCoeff, bigSmallSlidingFrictionCoeff;
    double bigBigRollingFrictionCoeff, smallSmallRollingFrictionCoeff, bigSmallRollingFrictionCoeff;
    double bigBigTorsionFrictionCoeff, smallSmallTorsionFrictionCoeff, bigSmallTorsionFrictionCoeff;
    // particle plastic and adhesive coefficients   // *** CARE: NO SEPARATE CASES FOR BIG AND SMALL YET! ***
    double particleMaxUnloadingStiffness;
    double particleCohesiveStiffness;
    double particlePlasticityDepth;
//    double particleAdhesionStiffness;
//    double particleMaxAdhesiveForce;
    
    // static geometry related variables
    double casingRadius;
    double casingHeight;
    double casingVolume;
    double totalHeight;
    
    // piston related variables
    double pistonHeight;
    double pistonVelocity;
    double pistonPressure;
    double targetPressure;
    double targetHeight;
    
    // data points related variables
    double *pressureStepsLoading;
    double *heightStepsLoading;
    double *calibratedLoadingK1;
    double *calibratedLoadingK2;
    double *calibratedLoadingKC;
    double *calibratedLoadingPhi;
    double *pressureStepsReloading;
    double *heightStepsReloading;
    double *calibratedReloadingK1;
    double *calibratedReloadingK2;
    double *calibratedReloadingKC;
    double *calibratedReloadingPhi;
    double *pressureStepsUnloading;
    double *heightStepsUnloading;
    double *calibratedUnloadingK1;
    double *calibratedUnloadingK2;
    double *calibratedUnloadingKC;
    double *calibratedUnloadingPhi;
    
    // data analysis variables
    double basePressure;
    double meanCoordinationNumber;
    double meanTotalRelativeOverlap;
    double maxTotalRelativeOverlap;
    double meanPistonRelativeOverlap;
    double maxPistonRelativeOverlap;
    double meanBaseRelativeOverlap;
    double maxBaseRelativeOverlap;
    
    // calibration routine variables
    double *pointerToCalibrationRoutineStiffnessArray;
    int calibrationStiffnessArrayLength;
    
    // simulation variables
    int stage;
    int nSets, nParticlesPerSet;
    double nBigPerSet, nSmallPerSet;
    double nBigParticles, nSmallParticles;
    double smallToBigMassRatio;
    double bulkPackingFractionPreCompression;
    int compressionStep;
    double timePlaceholder;
    int loadingCompressionInstances;
    double loadingCompressionInstanceDuration;
    int reloadingCompressionInstances;
    double reloadingCompressionInstanceDuration;
    int unloadingCompressionInstances;
    double unloadingCompressionInstanceDuration;
    
    // initialization variables
    double initRegionHeight;
    double volumeOfParticlesPerSet;
    int setsInserted;
    
    // Mercury-specific variables
    LinearPlasticViscoelasticFrictionSpecies *specieBig, *specieSmall, *specieWall;
    InfiniteWall basis, roof;
    CompressionPiston piston, *pistonPointer;
    AxisymmetricIntersectionOfWalls casing;
    SphericalParticle p0;
    
    // global variables
    bool verbose;
    bool biModal;
    bool isCompressing;
    bool pressureTargetMet;
    bool heightTargetMet;
    bool restartedFile;
    double cdatOutputTimeInterval;
    double energyRatioTolerance;
    bool loadingBranchDone = false;
    bool reloadingBranchDone = false;
    bool unloadingBranchdone = false;
    std::ofstream outputFile;
    std::ofstream historyDataFile;
    std::ofstream calibrationDataFile;
    std::ostringstream fileName;
};

#endif
