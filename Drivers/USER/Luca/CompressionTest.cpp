/*
    *** COMPRESSION TEST ***
 Particles are placed inside a cylindrical casing and let settle under gravity.
 When the system is at rest a piston is loaded and starts compressing the particles.
 The values of pressure towards the piston, the compression length and ther quantities are printed at each time step.
*/

#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <math.h>
#include "CompressionPiston.h"
#include <fstream>

/*
 TODOs:
 - Provide for settling problem only
 - Provide a loading restart file from a settling problem
 - Clean everything
 - The maxPressure in teh piston wall is unused. delete it
 - Output the relative overlap
 - Fix the setMaxTime() function properly
 - Unify the notations: change from pistonPointer -> getHeight() to pistonHeight and put the height evaluation in computePressure()
 - Check and correct the issue in the if (getKineticEnergy()/getElasticEnergy() < 1.0e-4) condition during compression. (Or maybe leave is as a "soft" check AFTER the simulation!?)
*/
 
class CompressionTest : public Mercury3D
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
    
    void actionsAfterTimeStep() override
    {
        // particle insertion loop
        if (stage == 2)
        {
            // other rounds of particle insertion
            if (setsInserted < nSets && getTime() > 0.05*setsInserted)
            {
                makeParticleSet();
            }
            
            // end of insertion, begin of the settling
            if (setsInserted == nSets)
            {
                std::cout << "Particle insertion terminated.\nSettling the particles...\n";
                stage++;
            }
        }
        
        // end of the settling, resetting times, begin of the compression
        if (stage == 3 && getKineticEnergy()/getElasticEnergy() < 1.e-4)
        {
            std::cout << "Particle settling terminated.\nTrimming the particles...\n";
            particleTrimmer();
            std::cout << "DONE\nCreating the piston...\n";
            makePiston();
            std::cout << "DONE\nCreating output data file...\n";
            makeOutputFile();
            
            // resetting time and tMax
            std::cout << "DONE\nResetting time...\n";
            setTime(0.0);
            std::cout << "DONE\nNew tMax = " << getTimeMax() << "\n";
            
            if (pressureDriven)
            {
                std::cout << "Creating the pressure steps...\n";
                makePressureStepsArray();
                std::cout << "DONE\nStarting the compression stage...\n";
            
                // begin the first compression
                targetPressure = pressureSteps[compressionStep];
                std::cout << "Target pressure: " << targetPressure << "\n";
                resetPistonVelocity();
                pressureTargetMet = false;
            }
            
            if (dataDriven)
            {
                std::cout << "Starting the compression stage...\n";
                targetHeight = heightSteps[compressionStep];
                targetPressure = pressureSteps[compressionStep];
                std::cout << "Target height: " << targetHeight << "\n";
                std::cout << "Target pressure: " << targetPressure << "\n";
                resetPistonVelocity();
                heightTargetMet = false;
                pressureTargetMet = false;
            }
            
            timePlaceholder = getTime();
            stage++;
        }
        
        // computes teh pressure and writes the data in the output file
        if (stage >= 4)
        {
            computePressure();
            
            if (fmod(getTime(),0.005) < getTimeStep())
            {
                writeDataToOutptFile();
            }
        }

        // goes on until all the compression and relaxation instances are ultimated
        if (stage == 4)
        {
            moveCompressionPiston();
        }
        
        // writes teh stiffness capibration data to file if data driven
        if (stage == 5)
        {
            if (dataDriven)
            {
                writeCalibrationData();
            }
        }
    }
    
    void actionAfterSolve()
    {
        // closes the stream to the output file
        outputFile.close();
        
        // deletes the dynamically allocated arrays
        delete[] heightSteps;
        delete[] pressureSteps;
        delete[] calibratedStiffnesses;
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
    
    // set the target piston compressing pressure (OVERLOADED)
    void setPistonTargetPressure(double pressure)
    {
        initialTargetPressure = pressure;
        compressionInstances = 1;
        reCompressionInstances = 0;
        pressureIncrementPerInstance = 0.0;
        pressureDecreasePerRelaxationInstance = 0.0;
        
        pressureDriven = true;
        dataDriven = false;
        
        if (verbose)
        {
            std::cout << "\tInitial target pressure: " << initialTargetPressure << "\n\n";
        }
    }
    void setPistonTargetPressure(double pressure, int nInstances, double pressureIncrement, double instanceDuration)
    {
        initialTargetPressure = pressure;
        compressionInstances = nInstances;
        reCompressionInstances = 0;
        pressureIncrementPerInstance = pressureIncrement;
        pressureDecreasePerRelaxationInstance = 0.0;
        compressionInstanceDuration = instanceDuration;
        
        pressureDriven = true;
        dataDriven = false;
        
        if (verbose)
        {
            std::cout << "\tInitial target pressure: " << initialTargetPressure << "\n";
            std::cout << "\tNumber of compression instances: " << compressionInstances << "\n";
            std::cout << "\tPressure increment per instance: " << pressureIncrementPerInstance << "\n";
            std::cout << "\tMaximum pressure: " << initialTargetPressure + (compressionInstances - 1)*pressureIncrementPerInstance << "\n";
            std::cout << "\tDuration of compression instance: " << compressionInstanceDuration << "\n\n";
        }
    }
    void setPistonTargetPressure(double pressure, int nInstances, double pressureIncrement, int nRelaxInstances, double pressureDecrease, double instanceDuration)
    {
        initialTargetPressure = pressure;
        compressionInstances = nInstances;
        reCompressionInstances = nRelaxInstances;
        pressureIncrementPerInstance = pressureIncrement;
        pressureDecreasePerRelaxationInstance = pressureDecrease;
        compressionInstanceDuration = instanceDuration;
        
        pressureDriven = true;
        dataDriven = false;
        
        if (verbose)
        {
            std::cout << "\tInitial target pressure: " << initialTargetPressure << "\n";
            std::cout << "\tNumber of compression instances: " << compressionInstances << "\n";
            std::cout << "\tPressure increment per instance: " << pressureIncrementPerInstance << "\n";
            std::cout << "\tNumber of re-compression instances: " << reCompressionInstances << "\n";
            std::cout << "\tPressure decrease per instance: " << pressureDecreasePerRelaxationInstance << "\n";
            std::cout << "\tMaximum pressure: " << initialTargetPressure + (compressionInstances - 1)*pressureIncrementPerInstance << "\n";
            std::cout << "\tMinimum pressure during re-compression: " << initialTargetPressure + (compressionInstances - 1)*pressureIncrementPerInstance - nRelaxInstances*pressureDecreasePerRelaxationInstance << "\n";
            std::cout << "\tDuration of compression instance: " << compressionInstanceDuration << "\n\n";
        }
    }
    
    // sets the number of compression cycles to be performed
    void setNumberOfCompressionCycles(int nCycles)
    {
        compressionCycles = nCycles;
    }
    
    // reads the data points and initializes the height-pressure array
    void setDataPoints(int n, double *heightArrayPointer, double *pressureArrayPointer, double duration)
    {
        compressionInstances = n;
        compressionInstanceDuration = duration;
        heightSteps = new double[n];
        pressureSteps = new double[n];
        calibratedStiffnesses = new double[n];
        
        pressureDriven = false;
        dataDriven = true;
        
        for (int i = 0; i < n; i++)
        {
            heightSteps[i] = heightArrayPointer[i];
            pressureSteps[i] = pressureArrayPointer[i];
        }
        
        if (verbose)
        {
            std::cout << "\tHeight-pressure steps:\n";
            for (int i = 0 ; i < n; i++) std::cout << "\t( " << heightSteps[i] << " , " << pressureSteps[i] << " )\n";
            std::cout << std::endl;
        }
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
        specieBig = new LinearViscoelasticFrictionSpecies;
        specieBig -> setDensity(particleDensityBig);
        specieBig -> setStiffnessAndRestitutionCoefficient(particleStiffnessBig, bigBigRestitutionCoeff, particleMassBig);
        
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
        specieSmall = new LinearViscoelasticFrictionSpecies;
        specieSmall -> setDensity(particleDensitySmall);
        specieSmall -> setStiffnessAndRestitutionCoefficient(particleStiffnessSmall, smallSmallRestitutionCoeff, particleMassSmall);
        
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
        specieWall = new LinearViscoelasticFrictionSpecies;
        specieWall -> setDensity(0.5*(particleDensityBig + particleDensitySmall));
        specieWall -> setStiffnessAndRestitutionCoefficient(wallStiffness, 1.0, (particleMassBig + particleMassSmall)/2.0);
        
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
                std::cout << "\tBIG-BIG stiffness and dissipation coefficients: " << specieBig -> getStiffness() << " " << specieBig -> getDissipation() << "\n";
                std::cout << "\tBIG-BIG friction coefficients: " << bigBigSlidingFrictionCoeff << " " << bigBigRollingFrictionCoeff << " " << bigBigTorsionFrictionCoeff << "\n";
                std::cout << "\tBIG-BIG tangential stiffnesses: " << specieBig -> getSlidingStiffness() << " " << specieBig -> getRollingStiffness() << " " << specieBig -> getTorsionStiffness() << "\n";
                std::cout << "\tBIG-BIG tangential dissipation coefficients: " << specieBig -> getSlidingDissipation() << " " << specieBig -> getRollingDissipation() << " " << specieBig -> getTorsionDissipation() << "\n";
                std::cout << "\tBIG-BIG collision time: " << std::setprecision(4) << specieBig -> getCollisionTime(particleMassBig) << "\n\n";
                
                std::cout << "\tSMALL-SMALL stiffness and dissipation coefficients: " << specieSmall -> getStiffness() << " " << specieSmall -> getDissipation() << "\n";
                std::cout << "\tSMALL-SMALL friction coefficients: " << smallSmallSlidingFrictionCoeff << " " << smallSmallRollingFrictionCoeff << " " << smallSmallTorsionFrictionCoeff << "\n";
                std::cout << "\tSMALL-SMALL tangential stiffnesses: " << specieSmall -> getSlidingStiffness() << " " << specieSmall -> getRollingStiffness() << " " << specieSmall -> getTorsionStiffness() << "\n";
                std::cout << "\tSMALL-SMALL tangential dissipation coefficients: " << specieSmall -> getSlidingDissipation() << " " << specieSmall -> getRollingDissipation() << " " << specieSmall -> getTorsionDissipation() << "\n";
                std::cout << "\tSMALL-SMALL collision time: " << std::setprecision(4) << specieSmall -> getCollisionTime(particleMassSmall) << "\n\n";
                
                std::cout << "\tBIG-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(specieBig, specieWall) -> getStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getDissipation() << "\n";
                std::cout << "\tBIG-WALL friction coefficients: " << bigWallSlidingFrictionCoeff << " " << bigWallRollingFrictionCoeff << " " << bigWallTorsionFrictionCoeff << "\n";
                std::cout << "\tBIG-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(specieBig, specieWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getTorsionStiffness() << "\n";
                std::cout << "\tBIG-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(specieBig, specieWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getTorsionDissipation() << "\n";
                std::cout << "\tBIG-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(specieBig, specieWall) -> getCollisionTime(particleMassBig) << "\n\n";
                
                std::cout << "\tSMALL-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getStiffness() << " " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getDissipation() << "\n";
                std::cout << "\tSMALL-WALL friction coefficients: " << smallWallSlidingFrictionCoeff << " " << smallWallRollingFrictionCoeff << " " << smallWallTorsionFrictionCoeff << "\n";
                std::cout << "\tSMALL-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getTorsionStiffness() << "\n";
                std::cout << "\tSMALL-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getTorsionDissipation() << "\n";
                std::cout << "\tSMALL-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(specieSmall, specieWall) -> getCollisionTime(particleMassSmall) << "\n\n";
                
                std::cout << "\tBIG-SMALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(specieBig, specieSmall) -> getStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieSmall) -> getDissipation() << "\n";
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
                std::cout << "\tPARTICLE stiffness and dissipation coefficients: " << specieBig -> getStiffness() << " " << specieBig -> getDissipation() << "\n";
                std::cout << "\tPARTICLE friction coefficients: " << bigBigSlidingFrictionCoeff << " " << bigBigRollingFrictionCoeff << " " << bigBigTorsionFrictionCoeff << "\n";
                std::cout << "\tPARTICLE tangential stiffnesses: " << specieBig -> getSlidingStiffness() << " " << specieBig -> getRollingStiffness() << " " << specieBig -> getTorsionStiffness() << "\n";
                std::cout << "\tPARTICLE tangential dissipation coefficients: " << specieBig -> getSlidingDissipation() << " " << specieBig -> getRollingDissipation() << " " << specieBig -> getTorsionDissipation() << "\n";
                std::cout << "\tPARTICLE collision time: " << std::setprecision(4) << specieBig -> getCollisionTime(particleMassBig) << "\n\n";
                
                std::cout << "\tPARTICLE-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(specieBig, specieWall) -> getStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getDissipation() << "\n";
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
        nSmallParticles = (int)(bulkPackingFractionPreCompression*casingVolume/particleVolumeSmall)*(smallToBigMassRatio*particleDensityBig/(particleDensitySmall + smallToBigMassRatio*particleDensityBig));
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
        // the piston velocity is set such that displacement = 0.0005*smallestParticleRatio/timeStep
        pistonVelocity = -0.0005*meanRadiusSmall*(1.0 - dispersitySmall)/getTimeStep();
        pistonPointer -> setVelocity(pistonVelocity);
        
        if (verbose)
        {
            std::cout << "\tRESETTING PISTON VELOCITY\n";
            std::cout << "\tPiston velocity: " << pistonVelocity << "\n";
            std::cout << "\tDisplacement per time step: " << 0.0005*meanRadiusSmall*(1.0 - dispersitySmall) << "\n";
        }
    }
    
    // computes the pressure acting against the piston
    void computePressure()
    {
        pistonPressure = 0.0;
        basePressure = 0.0;
        pistonForce = 0.0;
        
        for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
        {
            if ((*i) -> getI() -> getIndex() == pistonPointer -> getIndex())
            {
//                particleRadius = particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius();
//                pistonPressure -= ((*i) -> getForce()).Z/(pow(particleRadius,2.0) - pow(particleRadius - (*i) -> getOverlap(),2.0));
//                pistonPressureUnderestimated -= ((*i) -> getForce()).Z;
                pistonForce -= ((*i) -> getForce()).Z;
            }
            if ((*i) -> getI() -> getIndex() == basis.getIndex())
            {
                basePressure += ((*i) -> getForce()).Z;
            }
        }
        
        pistonPressure = pistonForce/(constants::pi*pow(casingRadius,2.0));
        basePressure /= constants::pi*pow(casingRadius,2.0);
    }
    
    // moves the piston
    void moveCompressionPiston()
    {
        // checks if the system is under compression by checking the direction of the piston velocity
        if (pistonVelocity <= 0.0) isCompressing = true;
        else isCompressing = false;
        
        if (pressureDriven)
        {
            // reversing and halving the velocity if needed
            if ((isCompressing && pistonPressure > targetPressure) || (!isCompressing && pistonPressure < targetPressure))
            {
                pistonVelocity *= -0.5;
                pistonPointer -> setVelocity(pistonVelocity);
            }
            
            //                // halves the velocity if Ekin/Epot > 1.0e-4
            //                if (getKineticEnergy()/getElasticEnergy() < 1.0e-4)
            //                {
            //                    pistonVelocity *= 0.5;
            //                    pistonPointer -> setVelocity(pistonVelocity);
            //                }
            
            // actually moving the piston
            pistonPointer -> movePiston(getTimeStep());
            
            if (!pressureTargetMet && fabs((pistonPressure - targetPressure)/(targetPressure + 0.01)) < 1.0e-4)
            {
                std::cout << "Pressure target reached, now chillin' for " << compressionInstanceDuration << " seconds\n";
                timePlaceholder = getTime();
                pressureTargetMet = true;
            }
            
            if (pressureTargetMet && getTime() > timePlaceholder + compressionInstanceDuration)
            {
                compressionStep++;
                
                if (compressionStep < 2*(compressionInstances + reCompressionInstances)*compressionCycles)
                {
                    targetPressure = pressureSteps[compressionStep];
                    std::cout << "New target pressure " << targetPressure << "\n";
                    resetPistonVelocity();
                    pressureTargetMet = false;
                }
                else
                {
                    std::cout << "Compression cycles finished...\n";
                    setTimeMax(getTime() + compressionInstanceDuration);
                    stage++;
                }
            }
        }
        
        if (dataDriven)
        {
            // reversing and halving the velocity if needed
            if ((isCompressing && pistonPointer -> getHeight() < targetHeight) || (!isCompressing && pistonPointer -> getHeight() > targetHeight))
            {
                pistonVelocity *= -0.5;
                pistonPointer -> setVelocity(pistonVelocity);
            }
            
            //                // halves the velocity if Ekin/Epot > 1.0e-4
            //                // only works if compression step is 1 or bigger, since the particles usually move during step 0 and the velocity would end up being zero
            //                if (compressionStep > 0 && getKineticEnergy()/getElasticEnergy() < 1.0e-4)
            //                {
            //                    pistonVelocity *= 0.5;
            //                    pistonPointer -> setVelocity(pistonVelocity);
            //                }
            
            // actually moving the piston
            pistonPointer -> movePiston(getTimeStep());
            
            // when the height is met stop the piston
            if (!heightTargetMet && fabs((pistonPointer -> getHeight() - targetHeight)/(targetHeight + 0.01)) < 1.0e-4)
            {
                std::cout << "Height target reached. Start tuning the spring stiffness.\n";
                //                    pistonPointer -> setVelocity(0.0);
                heightTargetMet = true;
                timePlaceholder = getTime();
            }
            
            // when the height is met and the pressure is NOT met start tuning hte particle stiffness
            if (heightTargetMet && !pressureTargetMet && fabs((pistonPressure - targetPressure)/(targetPressure + 0.01)) > 1.0e-4)
            {
                adjustStiffnesses();
            }
            
            // when the pressure is met as well stops the stiffness tuning
            if (heightTargetMet && !pressureTargetMet && fabs((pistonPressure - targetPressure)/(targetPressure + 0.01)) < 1.0e-4)
            {
                std::cout << "Pressure target reached, now chillin' for " << compressionInstanceDuration << " seconds\n";
//                calibratedStiffnesses[compressionStep] = particleStiffnessBig;
                timePlaceholder = getTime();
                pressureTargetMet = true;
            }
            
//            // goes on adjusting the stiffness
//            if (heightTargetMet && pressureTargetMet && getTime() < timePlaceholder + compressionInstanceDuration)
//            {
//                adjustStiffnesses();
//            }
            
            // when both height and pressure are met then waits and updates the target values
            if (heightTargetMet && pressureTargetMet && getTime() > timePlaceholder + compressionInstanceDuration)
            {
                calibratedStiffnesses[compressionStep] = particleStiffnessBig;
                compressionStep++;
                
                if (compressionStep < compressionInstances)
                {
                    targetHeight = heightSteps[compressionStep];
                    targetPressure = pressureSteps[compressionStep];
                    std::cout << "New target height: " << targetHeight << "\n";
                    std::cout << "New target pressure: " << targetPressure << "\n";
                    resetPistonVelocity();
                    heightTargetMet = false;
                    pressureTargetMet = false;
                }
                else
                {
                    std::cout << "Compression cycles finished...\n";
                    setTimeMax(getTime() + compressionInstanceDuration);
                    stage++;
                }
            }
            
        }
    }
    
    // unloads the piston
    void unloadPiston()
    {
        // sets the velocity such that it takes 1.0 second to get to the maximum height
        pistonVelocity = (totalHeight - pistonPointer -> getHeight())/1.0;
        pistonPointer -> setVelocity(pistonVelocity);
        
        if (pistonPointer -> getHeight() < totalHeight) pistonPointer -> movePiston(getTimeStep());
    }
    
    // creates the data output file and writes the first row
    void makeOutputFile()
    {
        std::ostringstream cdatName;
        std::cout.unsetf(std::ios::floatfield);
        cdatName << getName() << ".cdat";
        
        outputFile.open(cdatName.str(), std::ios::out);
        outputFile << "time \t Eel \t Ekin/Eel \t h \t v \t F \t P \t Pbase \t Pmax \t P/Pmax \t |(P - Pmax)/Pmax| " << std::endl;
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
        pistonForce << "   " <<
        pistonPressure << "   " <<
        basePressure << "   " <<
        targetPressure << "   " <<
        (pistonPressure + 0.01)/(targetPressure + 0.01) << "   " <<
        fabs((pistonPressure - targetPressure)/(targetPressure + 0.01)) << "   " <<
        std::endl;
    }
    
    // generates the target pressures according to instances and cycles
    void makePressureStepsArray()
    {
        pressureSteps = new double[2*(compressionInstances + reCompressionInstances)*compressionCycles];
        
        for (int i = 0; i < compressionCycles; i++)
        {
            pressureSteps[2*i*(compressionInstances + reCompressionInstances)] = initialTargetPressure;
            pressureSteps[2*(i + 1)*(compressionInstances + reCompressionInstances) - 1] = 0.0;
            
            for (int j = 1; j < 2*(compressionInstances + reCompressionInstances)-1; j++)
            {
                if (j < compressionInstances) pressureSteps[j + 2*i*(compressionInstances + reCompressionInstances)] = pressureSteps[j - 1] + pressureIncrementPerInstance;
                else if (j < compressionInstances + reCompressionInstances) pressureSteps[j + 2*i*(compressionInstances + reCompressionInstances)] = pressureSteps[j - 1] - pressureDecreasePerRelaxationInstance;
                else if (j < compressionInstances + 2*reCompressionInstances) pressureSteps[j + 2*i*(compressionInstances + reCompressionInstances)] = pressureSteps[j - 1] + pressureDecreasePerRelaxationInstance;
                else pressureSteps[j + 2*i*(compressionInstances + reCompressionInstances)] = pressureSteps[j - 1] - pressureIncrementPerInstance;
            }
        }
        
        if (verbose)
        {
            std::cout << "\tPressure steps:\n\t";
            for (int i = 0 ; i < 2*(compressionInstances + reCompressionInstances)*compressionCycles; i++) std::cout << pressureSteps[i] << " ";
            std::cout << std::endl;
        }
    }
    
    // adjusts the particle stiffness according to the pressure to be met and refreshes the species parameters involved
    void adjustStiffnesses()
    {
        if (fmod(getTime() - timePlaceholder, 50.0*getTimeStep()) < getTimeStep())
        {
            // if the pressure is too low increases the stiffness and vice-versa
            if (pistonPressure < targetPressure) {particleStiffnessBig *= 1.001; particleStiffnessSmall += 1.001;}
            else {particleStiffnessBig *= 0.999; particleStiffnessSmall *= 0.999;}
            
            specieBig -> setStiffnessAndRestitutionCoefficient(particleStiffnessBig, bigBigRestitutionCoeff, particleMassBig);
            specieBig -> setSlidingStiffness(particleStiffnessBig*2.0/7.0);
            specieBig -> setSlidingDissipation(specieBig -> getDissipation()*2.0/7.0);
            specieBig -> setRollingStiffness(particleStiffnessBig*2.0/7.0);
            specieBig -> setRollingDissipation(specieBig -> getDissipation()*2.0/7.0);
            specieBig -> setTorsionStiffness(particleStiffnessBig*2.0/7.0);
            specieBig -> setTorsionDissipation(specieBig -> getDissipation()*2.0/7.0);
            
            specieSmall -> setStiffnessAndRestitutionCoefficient(particleStiffnessSmall, smallSmallRestitutionCoeff, particleMassSmall);
            specieSmall -> setSlidingStiffness(particleStiffnessSmall*2.0/7.0);
            specieSmall -> setSlidingDissipation(specieSmall -> getDissipation()*2.0/7.0);
            specieSmall -> setRollingStiffness(particleStiffnessSmall*2.0/7.0);
            specieSmall -> setRollingDissipation(specieSmall -> getDissipation()*2.0/7.0);
            specieSmall -> setTorsionStiffness(particleStiffnessSmall*2.0/7.0);
            specieSmall -> setTorsionDissipation(specieSmall -> getDissipation()*2.0/7.0);
            
            speciesHandler.getMixedObject(specieBig, specieWall) -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessBig + wallStiffness), bigWallRestitutionCoeff, particleMassBig);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setSlidingStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setSlidingDissipation(speciesHandler.getMixedObject(specieBig, specieWall) -> getDissipation()*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setRollingStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setRollingDissipation(speciesHandler.getMixedObject(specieBig, specieWall) -> getDissipation()*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setTorsionStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieWall) -> setTorsionDissipation(speciesHandler.getMixedObject(specieBig, specieWall) -> getDissipation()*2.0/7.0);
            
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessSmall + wallStiffness), smallWallRestitutionCoeff, particleMassSmall);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setSlidingStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setSlidingDissipation(speciesHandler.getMixedObject(specieSmall, specieWall) -> getDissipation()*2.0/7.0);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setRollingStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setRollingDissipation(speciesHandler.getMixedObject(specieSmall, specieWall) -> getDissipation()*2.0/7.0);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setTorsionStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
            speciesHandler.getMixedObject(specieSmall, specieWall) -> setTorsionDissipation(speciesHandler.getMixedObject(specieSmall, specieWall) -> getDissipation()*2.0/7.0);
            
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessBig + particleStiffnessSmall), bigSmallRestitutionCoeff, 0.5*(particleMassBig + particleMassSmall));
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setSlidingStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setSlidingDissipation(speciesHandler.getMixedObject(specieBig, specieSmall) -> getDissipation()*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setRollingStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setRollingDissipation(speciesHandler.getMixedObject(specieBig, specieSmall) -> getDissipation()*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setTorsionStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
            speciesHandler.getMixedObject(specieBig, specieSmall) -> setTorsionDissipation(speciesHandler.getMixedObject(specieBig, specieSmall) -> getDissipation()*2.0/7.0);
        }
//        std::cout << "Stiffnesses: " << specieBig -> getStiffness() << " " << specieSmall -> getStiffness() << " " << specieWall -> getStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieWall) -> getStiffness() << " " << speciesHandler.getMixedObject(specieSmall, specieWall) -> getStiffness() << " " << speciesHandler.getMixedObject(specieBig, specieSmall) -> getStiffness() << "\n";
    }
    
    // writes the calibration data to file
    void writeCalibrationData()
    {
        std::ostringstream kdataName;
        std::cout.unsetf(std::ios::floatfield);
        kdataName << getName() << ".kdata";
        
        calibrationDataFile.open(kdataName.str(), std::ios::out);
        calibrationDataFile << "height \t pressure \t stiffness" << std::endl;
        for (int i = 0; i < compressionInstances; i++)
        {
            calibrationDataFile << heightSteps[i] << "\t" << pressureSteps[i] << "\t" << calibratedStiffnesses[i] << std::endl;
        }

        calibrationDataFile.close();
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
                std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() << ", Elastic Energy = " << std::setprecision(6) << std::left << std::setw(10) << getElasticEnergy() << ", Eratio = " << std::setprecision(6) << std::left << std::setw(10) << getKineticEnergy()/getElasticEnergy() << ", h = " << std::setprecision(6) << std::left << std::setw(10) << pistonPointer -> getHeight() << ", v = " << std::setprecision(6) << std::left << std::setw(10) << pistonPointer -> getVelocity() << ", F = " << std::setprecision(6) << std::left << std::setw(10) << pistonForce << ", P = " << std::setprecision(6) << std::left << std::setw(10) << pistonPressure << ", Pbase = " << std::setprecision(6) << std::left << std::setw(10) << basePressure << " , Pmax = " << std::setprecision(6) << std::left << std::setw(10) << targetPressure << ", P/Pmax = " << std::setprecision(6) << std::left << std::setw(10) << (pistonPressure + 0.01)/(targetPressure + 0.01) << ", |(P - Pmax)/Pmax| = " << std::setprecision(6) << std::left << std::setw(10) << fabs((pistonPressure - targetPressure)/(targetPressure + 0.01)) << ", k = " << std::setprecision(6) << std::left << std::setw(10) << specieBig -> getStiffness() << std::endl;
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
    
    // static geometry related variables
    double casingRadius;
    double casingHeight;
    double casingVolume;
    double totalHeight;
    
    // piston related variables
    double pistonHeight;
    double pistonVelocity;
    double pistonForce;
    double pistonPressure;
    double basePressure;
    double targetPressure;
    double targetHeight;
    int compressionInstances;
    int reCompressionInstances;
    double pressureIncrementPerInstance;
    double pressureDecreasePerRelaxationInstance;
    double particleContactSurface;
    double *pressureSteps;
    double *heightSteps;
    double *calibratedStiffnesses;
    double initialTargetPressure;
    
    // simulation variables
    int stage;
    int nSets, nParticlesPerSet;
    double nBigPerSet, nSmallPerSet;
    double nBigParticles, nSmallParticles;
    double smallToBigMassRatio;
    double bulkPackingFractionPreCompression;
    int compressionStep;
    double timePlaceholder;
    double compressionInstanceDuration;
    int compressionCycles;
    
    // initialization variables
    double initRegionHeight;
    double volumeOfParticlesPerSet;
    int setsInserted;
    
    // Mercury-specific variables
    LinearViscoelasticFrictionSpecies *specieBig, *specieSmall, *specieWall;
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
    bool pressureDriven = false;
    bool dataDriven = false;
    std::ofstream outputFile;
    std::ofstream calibrationDataFile;
    std::ostringstream fileName;
};

int main(int argc, char *argv[])
{
    double particleRadius = 0.001; // the radius of the small particle before rescaling
    double scalingFactor = 2.5; // the scaling factor for the particle sizes and cutoffs
    
    // the height and pressure arrays
    const int nDataPoints = 20;
    double heightArray [nDataPoints] = {0.0970949, 0.0968136, 0.0955561, 0.0955572, 0.0955581, 0.0955582, 0.0955585, 0.0930522, 0.0930498, 0.0930478, 0.0930471, 0.0918122, 0.0909015, 0.0909015, 0.0909015, 0.0909015, 0.0909015, 0.0909015, 0.0909015, 0.0909015};
    double pressureArray [nDataPoints] = {1730.74, 1978.18, 3031.96, 3008.31, 3006.01, 3005.56, 3004.63, 5464.74, 5478.79, 5485.37, 5487.87, 7029.78, 7990.82, 7994.32, 7993.63, 7994.36, 7993.89, 7994.07, 7994.06, 7994.05};
    
    CompressionTest cTest;
    cTest.setName("CompressionTest");
    
    // sets simulation parameters
    cTest.setTimeStep(2.0e-5);
    cTest.setTimeMax(100.0);
    cTest.setGravity(Vec3D(0., 0., -9.81));
    cTest.setSystemDimensions(3);
    cTest.setVerbose(true);
    
    // sets the number of saved time steps such that the output is printed every 0.01s
    cTest.setSaveCount(0.01/cTest.getTimeStep());
    
    // sets the particle intrinsic properties
    cTest.setParticleDensity(2000.0, 2000.0);
    cTest.setParticleRadiusAndDispersity(scalingFactor*2.0*particleRadius, 0.10, scalingFactor*particleRadius, 0.10);
    
    // sets the small-to-big total mass ratio
    cTest.setTotalMassRatio(0.1);
    
    // sets the stiffnesses
    cTest.setWallStiffness(5000.0);
    cTest.setParticleStiffness(500.0); // 1395.47   1500.00
    
    // sets the particle-wall friction coefficients
    cTest.setParticleWallSlidingFrictionCoeff(0.3, 0.3);
    cTest.setParticleWallRollingFrictionCoeff(0.01, 0.01);
    cTest.setParticleWallTorsionFrictionCoeff(0.0, 0.0);
    
    // sets the particle-particle friction coefficients
    cTest.setParticleParticleSlidingFrictionCoeff(0.5, 0.8, 0.7);
    cTest.setParticleParticleRollingFrictionCoeff(0.05, 0.1, 0.75);
    cTest.setParticleParticleTorsionFrictionCoeff(0.0, 0.0, 0.0);
    
    // sets the particle restitution coefficients
    cTest.setParticleWallRestitutionCoeff(0.8, 0.3);
    cTest.setParticleParticleRestitutionCoeff(0.5, 0.2, 0.3);
    
    // sets the casing dimensions
    cTest.setCasingDimensions(0.04, 0.1);
    
    // sets the bulk packing fraction of the powder prior to compression
    cTest.setbulkPackingFractionPreCompression(0.70);
    
//    // sets the compression pressure
//    cTest.setPistonTargetPressure(3000.0, 3, 2500.0, 2, 2000.0, 0.5);
//    // sets the number of compression-relaxation cycles
//    cTest.setNumberOfCompressionCycles(1);
    
    cTest.setDataPoints(nDataPoints, heightArray, pressureArray, 0.5);
    
    // sets additional built-in arguments for the xballs visualization
    cTest.setXBallsAdditionalArguments("-h 800 -p 10 -o 200 -3dturn 3");
    
    cTest.solve();
    
    return 0;
}


//      JUNKYARD
/*
 // UNCOMMENT FOR SINGLE SETTLING STEP
 
 void computeNumberOfLevelsAndParticles()
 {
 // computes the total number of big ans small particles needed
 nSmallParticles = (int)(bulkPackingFractionPreCompression*casingVolume/particleVolumeSmall)*(smallToBigMassRatio*particleDensityBig/(particleDensitySmall + smallToBigMassRatio*particleDensityBig));
 nBigParticles = (int)((bulkPackingFractionPreCompression*casingVolume - nSmallParticles*particleVolumeSmall)/particleVolumeBig);
 
 // computes the number of particle rings per particle level
 numberOfRings = (int)(2.0*casingRadius/(2.1*meanRadiusBig*(1.0 + dispersityBig)));
 if (numberOfRings%2) {numberOfRings = (numberOfRings - 1)/2;}
 else {numberOfRings /= 2;}
 
 // computation of ring distance from centre and number of particles on each ring
 ringRadialDistance = new double[numberOfRings];
 particlesPerRing = new int[numberOfRings];
 int numberOfParticlesPerLevel = 0;
 
 for (int i = 0; i < numberOfRings; i++)
 {
 if (i == 0)
 {
 if (numberOfRings%2)
 {
 ringRadialDistance[0] = 0.0;
 particlesPerRing[0] = 1;
 }
 else
 {
 ringRadialDistance[0] = 2.05*meanRadiusBig*(1.0 + dispersityBig);
 particlesPerRing[0] = (int)(constants::pi*ringRadialDistance[0]/(1.05*meanRadiusBig*(1.0 + dispersityBig))) - 1;
 }
 }
 else
 {
 ringRadialDistance[i] = ringRadialDistance[i-1] + 2.1*meanRadiusBig*(1.0 + dispersityBig);
 particlesPerRing[i] = (int)(constants::pi*ringRadialDistance[i]/(1.05*meanRadiusBig*(1.0 + dispersityBig))) - 1;
 }
 
 numberOfParticlesPerLevel += particlesPerRing[i];
 }
 
 // computation of number of levels needed
 numberOfLevels = (int)((nSmallParticles + nBigParticles)/numberOfParticlesPerLevel) + 1;
 
 // computation of number of big/small particles per level
 if (nSmallParticles > nBigParticles)
 {
 nBigPerLevel = (int)(numberOfParticlesPerLevel*nBigParticles/(nSmallParticles + nBigParticles));
 nSmallPerLevel = numberOfParticlesPerLevel - nBigPerLevel;
 }
 else
 {
 nSmallPerLevel = (int)(numberOfParticlesPerLevel*nSmallParticles/(nSmallParticles + nBigParticles));
 nBigPerLevel = numberOfParticlesPerLevel - nSmallPerLevel;
 }
 }
 

// UNCOMMENT FOR SINGLE SETTLING STEP

 void makeParticles()
 {
 int smallCounter, bigCounter;
 double angularDistance;
 double randomSeed;
 p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
 
 for (int nLevel = 0; nLevel < numberOfLevels; nLevel++)
 {
 smallCounter = 0;
 bigCounter = 0;
 
 for (int nRing = 0; nRing < numberOfRings; nRing++)
 {
 angularDistance = 2.0*constants::pi/particlesPerRing[nRing];
 randomSeed = random.getRandomNumber(0.0,angularDistance);
 
 for (int nParticles = 0; nParticles < particlesPerRing[nRing]; nParticles++)
 {
 p0.setPosition(Vec3D(ringRadialDistance[nRing]*cos(nParticles*angularDistance + randomSeed), ringRadialDistance[nRing]*sin(nParticles*angularDistance + randomSeed), 0.5*casingHeight + 2.1*meanRadiusBig*(1.0 + dispersityBig)*nLevel));
 
 if (smallCounter < nSmallPerLevel)
 {
 if (bigCounter < nBigPerLevel)
 {
 if (random.getRandomNumber(0.0,1.0) < 0.5)
 {
 p0.setRadius(meanRadiusBig*(1.0 + dispersityBig*random.getRandomNumber(-1.0,1.0)));
 p0.setSpecies(specieBig);
 bigCounter++;
 }
 else
 {
 p0.setRadius(meanRadiusSmall*(1.0 + dispersitySmall*random.getRandomNumber(-1.0,1.0)));
 p0.setSpecies(specieSmall);
 smallCounter++;
 }
 }
 else
 {
 p0.setRadius(meanRadiusSmall*(1.0 + dispersitySmall*random.getRandomNumber(-1.0,1.0)));
 p0.setSpecies(specieSmall);
 smallCounter++;
 }
 }
 else
 {
 p0.setRadius(meanRadiusBig*(1.0 + dispersityBig*random.getRandomNumber(-1.0,1.0)));
 p0.setSpecies(specieBig);
 bigCounter++;
 }
 
 particleHandler.copyAndAddObject(p0);
 }
 }
 }
 
 delete[] ringRadialDistance;
 delete[] particlesPerRing;
 }
 


 //    // UNCOMMENT FOR SINGLE SETTLING STEP
 //    double *ringRadialDistance;
 //    int *particlesPerRing;
 //    int numberOfParticlesPerLevel;
 //    int numberOfRings;
 //    int numberOfLevels;
 //    int nSmallPerLevel, nBigPerLevel;


 */






