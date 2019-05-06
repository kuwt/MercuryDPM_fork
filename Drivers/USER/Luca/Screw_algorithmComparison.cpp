#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Walls/InfiniteWall.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include "Helicoid05.h"
#include <fstream>

/*
 Last update : 30.5.17
 
 ToDo:
 MAJOR
 
 MINOR

 WHATEVER
 
 */

class ScrewAlgorithmComparison : public Mercury3D
{
private:
    
    void setupInitialConditions()
    {
        stage = 1;
        
        // computes the mean mass and volume of the particles
        std::cout << "Computing particle mean mass and volume...\n";
        setParticleMassAndVolume();
        std::cout << "DONE!\n";
        
        // sets the species of the pbjects
        std::cout << "Setting up species...\n";
        setSpecies();
        std::cout << "DONE!\n";
        
        // sets the box around the screw
        std::cout << "Setting up the filling box dimensions...";
        setFillingBoxDimensions();
        std::cout << "\nFilling box dimensions: (" << fillingBoxMax.X - fillingBoxMin.X << " , " << fillingBoxMax.Y - fillingBoxMin.Y << " , " << fillingBoxMax.Z - fillingBoxMin.Z << ").\n";
        std::cout << "Filling box volume: " << fillingBoxVolume << ".\n";
        
        // sets the number of particles needed
        std::cout << "Setting up the number of particles needed (accounts for the shaft volume)...";
        setNumberOfParticles();
        std::cout << "\nMono-modal particle size distribution. Polydispersity ratio: " << particleDispersity << ".";
        std::cout << "\nMean particle radius: " << particleRadius << ".";
        std::cout << "\nNumber of particles needed: " << numberOfParticles << ".\n";
        
        // sets the lattice properties where the particles will be initially loaded
        std::cout << "Setting up the loading lattice properties...\n";
        setloadingBoxGridProperties();
        std::cout << "DONE!\n";
        
        // sets the size of the initial particle loading box on top of the filling box
        std::cout << "Setting up the loading box dimensions...";
        setLoadingBoxDimensions();
        std::cout << "\nLoading box dimensions: (" << loadingBoxMax.X - loadingBoxMin.X << " , " << loadingBoxMax.Y - loadingBoxMin.Y << " , " << loadingBoxMax.Z - loadingBoxMin.Z << ").\n";
        std::cout << "Loading lattice grid size: (" << loadingLatticeGridSizeX << " , " << loadingLatticeGridSizeY << " , " << loadingLatticeGridSizeZ << ").\n";
        std::cout << "Loading lattice sites distance: (" << loadingLatticeSitesDistanceHorizontal << " , " << loadingLatticeSitesDistanceVertical << " , " << loadingLatticeSitesDistanceHorizontal << ").\n";
        std::cout << "Loading lattice offsets: (" << loadingLatticeOffset.X << " , " << loadingLatticeOffset.Y << " , " << loadingLatticeOffset.Z << ").\n";
        
        // sets the simulation domain
        std::cout << "Setting up the simulation boundaries...";
        setBoundaries();
        std::cout << "\nSimulation boundaries: ((" << getXMin() << " , " << getXMax() << ") , (" << getYMin() << " , " << getYMax() << ") , (" << getZMin() << " , " << getZMax() << ")).\n";

        // creates the particles in the lattice
        std::cout << "Creating the particles... ";
        makeParticles();
        std::cout << "DONE!\n";
        
        // creates the external walls and the boundary condition
        std::cout << "Creating the boundary walls... ";
        makeBoundaries();
        std::cout << "DONE!\n";
        
        // creates the screw
        std::cout << "Creating the screw... ";
        makeScrew();
        std::cout << "DONE!\n";
    }
    
    void actionsBeforeTimeStep()
    {
        // stage 1: the particles are loaded. Once they are still the external casing is loaded and the particles outside removed
        if (stage == 1 && getKineticEnergy()/getElasticEnergy() < 1.e-4)
        {
            std::cout << "Initial number of particles: " << particleHandler.getNumberOfObjects() << ".\n";
            std::cout << "Creating the screw casing and cutting the particles outside... ";
            makeCasing();
            std::cout << "DONE!\n New number of particles: " << particleHandler.getNumberOfObjects() << ".\n";
            std::cout << "Relaxation of the system...\n";
            
            stage++;
        }
        
        // stage 2: the particles are cut according to the desired filling ratio
        if (stage == 2 && getKineticEnergy()/getElasticEnergy() < 1.e-4)
        {
            std::cout << "Cutting the particles according to the filling ratio...\n";
            makeFillingRatio();
            std::cout << "DONE!\nFinal number of particles: " << particleHandler.getNumberOfObjects() << ".\n";
            std::cout << "Relaxation of the system...\n";
            
            stage++;
        }
        
        // stage 3: the max time is reset according to the desired screw running time
        if (stage == 3 && getKineticEnergy()/getElasticEnergy() < 1.e-4)
        {
            std::cout << "Resetting time...\n";
            setTime(0.0);
            std::cout << "Resetting time max...\n";
            setTimeMax(numberOfvelocityIncrements*durationOfVelocityIncrement);
            std::cout << "Creating the particle list array...\n";
            makeParticleList();
            std::cout << "Creating the .cdat output file...\n";
            makeOutputFile();
            std::cout << "Starting the screw rotation...\n";
            
            stage++;
        }
        
        // stage 4: the screw is set in motion
        if (stage == 4)
        {
            makeRotation();
            makeVelocityIncrement();
            
            
            if (fmod(getTime(),cdatOutputTimeInterval) < getTimeStep())
            {
                makeDataAnalysis();
                writeDataToOutptFile();
                refreshParticleList();
            }
        }
    }
    
    void actionAfterSolve()
    {
        // closes the stream to the output file
        outputFile.close();
        
        // deletes the dynamically allocated memory
        delete[] velocitySteps;
        delete[] particleList;
    }
    
    
public:
//  ----- FUNCTIONS CALLED IN THE MAIN -----
// particle properties
    // sets particle density
    void setParticleDensity(double rho)
    {
        particleDensity = rho;
    }
    
    // sets particle radius and dispersity
    void setParticleRadiusAndDispersity(double r, double disp)
    {
        particleRadius = r;
        particleDispersity = disp;
    }
    
    
// screw properties
    // sets the screw geometry parameters
    void setScrewGeometry(double length, double casing, double radius, double shaft, double n, double thickness)
    {
        screwLength = length;
        screwCasingRadius = casing;
        screwBladeRadius = radius;
        screwShaftRadius = shaft;
        screwNumberOfTurns = n;
        screwThickness = thickness;
    }
    
    // sets the screw origin
    void setScrewOrigin(Vec3D origin)
    {
        screwOrigin = origin;
    }
    
    // sets the scew filling ratio
    void setScrewFillingRatio(double fR)
    {
        fillingRatio = fR;
    }
    
    
// interaction parameters
    // sets the particle and the wall stiffness
    void setStiffnesses(double kP, double kW)
    {
        particleStiffness = kP;
        wallStiffness = kW;
    }
    
    // sets the sliding friction coefficients of the species
    void setSlidingFrictionCoefficients(double ppsf, double pwsf)
    {
        particleParticleSlidingFriction = ppsf;
        particleWallSlidingFriction = pwsf;
    }
    
    // sets the rolling friction coefficients of the species
    void setRollingFrictionCoefficients(double pprf, double pwrf)
    {
        particleParticleRollingFriction = pprf;
        particleWallRollingFriction = pwrf;
    }
    
    // sets the torsion friction coefficients of the species
    void setTorsionFrictionCoefficients(double pptf, double pwtf)
    {
        particleParticleTorsionFriction = pptf;
        particleWallTorsionFriction = pwtf;
    }
    
    // sets the restitution coefficients of the species
    void setRestitutionCoefficients(double pprc, double pwrc)
    {
        particleParticleRestitutionCoefficient = pprc;
        particleWallRestitutionCoefficient = pwrc;
    }
    
    
// global parameters
    // sets the packing fraction of the particles after the casing insertion cut
    void setPackingFraction(bool pfSetting, double pf)
    {
        fixedPackingFraction = pfSetting;
        packingFraction = pf;
    }
    
    // sets the ratio between the filling box size and the screw casing diameter
    void setBoxToScrewWidthRatio(bool swrSetting, double ratio)
    {
        fixedBoxToScrewWidthRatio = swrSetting;
        boxToScrewWidthRatio = ratio;
    }
    
    // sets the interparticle distance inside of the loading lattice normalized by the mean radius
    void setInterParticleLoadingLatticeRelativeDistance(double relDist)
    {
        interParticleLoadingLatticeRelativeDistance = relDist;
        interParticleLoadingDistance = interParticleLoadingLatticeRelativeDistance*particleRadius*(1.0 + particleDispersity);
    }
    
    // generates the set of velocity steps and their duration
    void setVelocitySteps(int n, double *velocityArrayPointer, double duration)
    {
        numberOfvelocityIncrements = n;
        durationOfVelocityIncrement = duration;
        velocitySteps = new double[n];
        
        for (int i = 0; i < n; i++)
        {
            velocitySteps[i] = velocityArrayPointer[i];
        }
        
        if (verbose)
        {
            std::cout << "\tVelocity steps:\n\t";
            for (int i = 0 ; i < n; i++) std::cout << velocitySteps[i] << " ";
            std::cout << std::endl;
        }
    }
    
    
// development hacks
    // function to hack the number of loaded particles for debugging purpouses
    void developmentHackParticles(bool hack, int layers)
    {
        particleHack = hack;
        particleHackLayers = layers;
    }
    
    
//  ----- FUNCTIONS CALLED IN THE CLASS -----
// set functions
    // sets the filling box size
    void setFillingBoxDimensions()
    {
        if(!fixedBoxToScrewWidthRatio)  // automated setting
        {
            fillingBoxMin.X = screwOrigin.X - (screwCasingRadius + 4.0*particleRadius);
            fillingBoxMin.Y = screwOrigin.Y - (screwCasingRadius + 4.0*particleRadius);
            fillingBoxMin.Z = screwOrigin.Z;
            
            fillingBoxMax.X = screwOrigin.X + (screwCasingRadius + 4.0*particleRadius);
            fillingBoxMax.Y = screwOrigin.Y + (screwCasingRadius + 4.0*particleRadius);
            fillingBoxMax.Z = screwOrigin.Z + screwLength;
        }
        else    // manual setting of teh ratio
        {
            fillingBoxMin.X = screwOrigin.X - boxToScrewWidthRatio*screwCasingRadius;
            fillingBoxMin.Y = screwOrigin.Y - boxToScrewWidthRatio*screwCasingRadius;
            fillingBoxMin.Z = screwOrigin.Z;
            
            fillingBoxMax.X = screwOrigin.X + boxToScrewWidthRatio*screwCasingRadius;
            fillingBoxMax.Y = screwOrigin.Y + boxToScrewWidthRatio*screwCasingRadius;
            fillingBoxMax.Z = screwOrigin.Z + screwLength;
        }
        
        fillingBoxVolume = (fillingBoxMax.X - fillingBoxMin.X)*(fillingBoxMax.Y - fillingBoxMin.Y)*(fillingBoxMax.Z - fillingBoxMin.Z);
    }
    
    // sets the loading box parameters
    void setloadingBoxGridProperties()      // explain how this thing works
    {
        loadingLatticeSitesDistanceHorizontal = 2.0*particleRadius*(1.0 + particleDispersity) + interParticleLoadingDistance;
        loadingLatticeSitesDistanceVertical = 0.5*sqrt(3.0)*loadingLatticeSitesDistanceHorizontal;
        
        loadingLatticeGridSizeX = (int)((fillingBoxMax.X - fillingBoxMin.X - interParticleLoadingDistance)/loadingLatticeSitesDistanceHorizontal);
        loadingLatticeGridSizeZ = (int)((fillingBoxMax.Z - fillingBoxMin.Z - interParticleLoadingDistance)/loadingLatticeSitesDistanceHorizontal);
        
        loadingLatticeOffset.X = 0.5*(fillingBoxMax.X - fillingBoxMin.X - (loadingLatticeSitesDistanceHorizontal*loadingLatticeGridSizeX + interParticleLoadingDistance));
        loadingLatticeOffset.Y = 0.0;
        loadingLatticeOffset.Z = 0.5*(fillingBoxMax.Z - fillingBoxMin.Z - (loadingLatticeSitesDistanceHorizontal*loadingLatticeGridSizeZ + interParticleLoadingDistance));
        
        int nParticlesPerOddLatticePlane = loadingLatticeGridSizeX*loadingLatticeGridSizeZ;
        int nParticlesPerEvenLatticePlane = (loadingLatticeGridSizeX - 1)*(loadingLatticeGridSizeZ - 1);
        
        loadingLatticeGridSizeY = (int)(numberOfParticles/(nParticlesPerEvenLatticePlane + nParticlesPerOddLatticePlane));
        if (loadingLatticeGridSizeY*(nParticlesPerEvenLatticePlane+nParticlesPerOddLatticePlane) - numberOfParticles > nParticlesPerOddLatticePlane)
        {
            loadingLatticeGridSizeY = 2*(loadingLatticeGridSizeY + 1);
        }
        else
        {
            loadingLatticeGridSizeY = 2*loadingLatticeGridSizeY + 1;
        }
    }
    
    // sets the loading box size
    void setLoadingBoxDimensions()
    {
        loadingBoxMin.X = fillingBoxMin.X;
        loadingBoxMin.Y = fillingBoxMax.Y;
        loadingBoxMin.Z = fillingBoxMin.Z;
        
        loadingBoxMax.X = fillingBoxMax.X;
        loadingBoxMax.Y = fillingBoxMax.Y + loadingLatticeSitesDistanceVertical*(loadingLatticeGridSizeY + 1);
        loadingBoxMax.Z = fillingBoxMax.Z;
    }
    
    // sets the simulation domain
    void setBoundaries()
    {
        setXMin(fillingBoxMin.X);
        setYMin(fillingBoxMin.Y);
        setZMin(fillingBoxMin.Z);
        
        setXMax(loadingBoxMax.X);
        setYMax(loadingBoxMax.Y);
        setZMax(loadingBoxMax.Z);
    }
    
    // sets the particle mass and volume from density and radius
    void setParticleMassAndVolume()
    {
        particleVolume = 4.*constants::pi*pow(particleRadius,3.)/3.;
        particleMass = particleDensity*particleVolume;
    }
    
    // sets the objects species
    void setSpecies()
    {
        speciesHandler.clear();
        
        // particle
        specieParticle = new LinearViscoelasticFrictionSpecies;
        specieParticle -> setDensity(particleDensity);
        specieParticle -> setStiffnessAndRestitutionCoefficient(particleStiffness, particleParticleRestitutionCoefficient, particleMass);
        
        specieParticle -> setSlidingDissipation(specieParticle -> getDissipation()*2.0/7.0);
        specieParticle -> setSlidingStiffness(particleStiffness*2.0/7.0);
        specieParticle -> setSlidingFrictionCoefficient(particleParticleSlidingFriction);
        
        specieParticle -> setRollingDissipation(specieParticle -> getDissipation()*2.0/7.0);
        specieParticle -> setRollingStiffness(particleStiffness*2.0/7.0);
        specieParticle -> setRollingFrictionCoefficient(particleParticleRollingFriction);
        
        specieParticle -> setTorsionDissipation(specieParticle -> getDissipation()*2.0/7.0);
        specieParticle -> setTorsionStiffness(particleStiffness*2.0/7.0);
        specieParticle -> setTorsionFrictionCoefficient(particleParticleRollingFriction);
        speciesHandler.addObject(specieParticle);
        
        // wall
        specieWall = new LinearViscoelasticFrictionSpecies;
        specieWall -> setDensity(particleDensity);
        specieWall -> setStiffnessAndRestitutionCoefficient(wallStiffness, 1.0, particleMass);
        
        specieWall -> setSlidingDissipation(specieWall -> getDissipation()*2.0/7.0);
        specieWall -> setSlidingStiffness(wallStiffness*2.0/7.0);
        specieWall -> setSlidingFrictionCoefficient(0.0);
        
        specieWall -> setRollingDissipation(specieWall -> getDissipation()*2.0/7.0);
        specieWall -> setRollingStiffness(wallStiffness*2.0/7.0);
        specieWall -> setRollingFrictionCoefficient(0.0);
        
        specieWall -> setTorsionDissipation(specieWall -> getDissipation()*2.0/7.0);
        specieWall -> setTorsionStiffness(wallStiffness*2.0/7.0);
        specieWall -> setTorsionFrictionCoefficient(0.0);
        speciesHandler.addObject(specieWall);
        
        // particle-wall mixed
        auto specieMixedParticleWall = speciesHandler.getMixedObject(specieParticle, specieWall);
        specieMixedParticleWall -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffness + wallStiffness), particleWallRestitutionCoefficient, particleMass);
        
        specieMixedParticleWall -> setSlidingDissipation(specieMixedParticleWall -> getDissipation()*2.0/7.0);
        specieMixedParticleWall -> setSlidingStiffness(0.5*(particleStiffness + wallStiffness)*2.0/7.0);
        specieMixedParticleWall -> setSlidingFrictionCoefficient(particleWallSlidingFriction);
        
        specieMixedParticleWall -> setRollingDissipation(specieMixedParticleWall -> getDissipation()*2.0/7.0);
        specieMixedParticleWall -> setRollingStiffness(0.5*(particleStiffness + wallStiffness)*2.0/7.0);
        specieMixedParticleWall -> setRollingFrictionCoefficient(particleWallRollingFriction);
        
        specieMixedParticleWall -> setTorsionDissipation(specieMixedParticleWall -> getDissipation()*2.0/7.0);
        specieMixedParticleWall -> setTorsionStiffness(0.5*(particleStiffness + wallStiffness)*2.0/7.0);
        specieMixedParticleWall -> setTorsionFrictionCoefficient(particleWallTorsionFriction);
        
        // If the option is turned on, it prints all the above interaction parameters
        if (verbose)
        {
            std::cout << "\tPARTICLE stiffness and dissipation coefficients: " << specieParticle -> getStiffness() << " " << specieParticle -> getDissipation() << "\n";
            std::cout << "\tPARTICLE friction coefficients: " << particleParticleSlidingFriction << " " << particleParticleRollingFriction << " " << particleParticleTorsionFriction << "\n";
            std::cout << "\tPARTICLE tangential stiffnesses: " << specieParticle -> getSlidingStiffness() << " " << specieParticle -> getRollingStiffness() << " " << specieParticle -> getTorsionStiffness() << "\n";
            std::cout << "\tPARTICLE tangential dissipation coefficients: " << specieParticle -> getSlidingDissipation() << " " << specieParticle -> getRollingDissipation() << " " << specieParticle -> getTorsionDissipation() << "\n";
            std::cout << "\tPARTICLE collision time: " << std::setprecision(4) << specieParticle -> getCollisionTime(particleMass) << "\n\n";
            
            std::cout << "\tPARTICLE-WALL stiffness and dissipation coefficients: " << specieMixedParticleWall -> getStiffness() << " " << specieMixedParticleWall -> getDissipation() << "\n";
            std::cout << "\tPARTICLE-WALL friction coefficients: " << particleWallSlidingFriction << " " << particleWallRollingFriction << " " << particleWallTorsionFriction << "\n";
            std::cout << "\tPARTICLE-WALL tangential stiffnesses: " << specieMixedParticleWall -> getSlidingStiffness() << " " << specieMixedParticleWall -> getRollingStiffness() << " " << specieMixedParticleWall -> getTorsionStiffness() << "\n";
            std::cout << "\tPARTICLE-WALL tangential dissipation coefficients: " << specieMixedParticleWall -> getSlidingDissipation() << " " << specieMixedParticleWall -> getRollingDissipation() << " " << specieMixedParticleWall -> getTorsionDissipation() << "\n";
            std::cout << "\tPARTICLE-WALL collision time: " << std::setprecision(4) << specieMixedParticleWall -> getCollisionTime(particleMass) << "\n\n";
            
            std::cout << "\tPARTICLE-PARTICLE collision time / TIME STEP: " << std::setprecision(4) << (specieParticle -> getCollisionTime(particleMass))/getTimeStep() << "\n";
            std::cout << "\tPARTICLE-WALL collision time / TIME STEP: " << std::setprecision(4) << (specieMixedParticleWall -> getCollisionTime(particleMass))/getTimeStep() << "\n\n";
        }
    }
    
    // sets the number of particles
    void setNumberOfParticles() // safe constant * assumed packing fraction * (box volume - shaft volume) / mean particle volume
    {
        numberOfParticles = (int)(1.2*0.6*(fillingBoxVolume - constants::pi*screwLength*pow(screwShaftRadius,2.))/particleVolume);
    }
    
    
// make functions
    // makes the external walls and boundary condition
    void makeBoundaries()
    {
        wallHandler.clear();
        wall.setSpecies(specieWall);
        
        wall.set(Vec3D(-1.,0.,0.),Vec3D(getXMin(),0.,0.));
        wallHandler.copyAndAddObject(wall);
        wall.set(Vec3D(1.,0.,0.),Vec3D(getXMax(),0.,0.));
        wallHandler.copyAndAddObject(wall);
        
        wall.set(Vec3D(0,-1.,0.),Vec3D(0.,getYMin(),0.));
        wallHandler.copyAndAddObject(wall);
        wall.set(Vec3D(0.,1.,0.),Vec3D(0.,getYMax(),0.));
        wallHandler.copyAndAddObject(wall);
        
        boundary.set(Vec3D(0.,0.,1.), getZMin(), getZMax());
        boundaryHandler.copyAndAddObject(boundary);
    }
    
    // loads the particles in the loading lattice
    void makeParticles()
    {
        double xPos, yPos, zPos;
        
        p0.setSpecies(specieParticle);
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        
        if (particleHack) loadingLatticeGridSizeY = particleHackLayers;
        
        for (int i=0; i<loadingLatticeGridSizeY; i++)
        {
            yPos = loadingBoxMin.Y + loadingLatticeOffset.Y + i*loadingLatticeSitesDistanceVertical;
            
            if (!(i%2))
            {
                for (int j=0; j<loadingLatticeGridSizeX; j++)
                {
                    xPos = loadingBoxMin.X + loadingLatticeOffset.X + particleRadius*(1.0 + particleDispersity) + interParticleLoadingDistance + j*loadingLatticeSitesDistanceHorizontal;
                    for (int k=0; k<loadingLatticeGridSizeZ; k++)
                    {
                        zPos = loadingBoxMin.Z + loadingLatticeOffset.Z + particleRadius*(1.0 + particleDispersity) + interParticleLoadingDistance + k*loadingLatticeSitesDistanceHorizontal;
                        
                        p0.setPosition(Vec3D(xPos,yPos,zPos));
                        p0.setRadius(particleRadius*(1.0 + particleDispersity*random.getRandomNumber(-1.,1.)));
                        particleHandler.copyAndAddObject(p0);
                    }
                }
            }
            else
            {
                for (int j=0; j<loadingLatticeGridSizeX-1; j++)
                {
                    xPos = loadingBoxMin.X + loadingLatticeOffset.X + 0.5*interParticleLoadingDistance + (j+1)*loadingLatticeSitesDistanceHorizontal;
                    for (int k=0; k<loadingLatticeGridSizeZ-1; k++)
                    {
                        zPos = loadingBoxMin.Z + loadingLatticeOffset.Z + 0.5*interParticleLoadingDistance + (k+1)*loadingLatticeSitesDistanceHorizontal;
                        
                        p0.setPosition(Vec3D(xPos,yPos,zPos));
                        p0.setRadius(particleRadius*(1.0 + particleDispersity*random.getRandomNumber(-1.,1.)));
                        particleHandler.copyAndAddObject(p0);
                    }
                }
            }
        }
    }

    // makes the screw
    void makeScrew()
    {
        // screw parameters
        helicoid.setSpecies(specieWall);
        helicoid.set(screwOrigin, screwLength, screwBladeRadius, screwShaftRadius, screwNumberOfTurns, screwAngularVelocity, screwThickness, true);
        helicoidPointer = wallHandler.copyAndAddObject(helicoid);
    }
    
    // makes the screw casing and cuts the particles accordingly
    void makeCasing()
    {
        // if the final packing fraction is not user defined the particles with r >= screwCasingRadius - .95*particleRadius are eliminated
        if(!fixedPackingFraction)
        {
            for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
            {
                if (pow(particleHandler.getObject(i) -> getPosition().X - screwOrigin.X,2.) +
                    pow(particleHandler.getObject(i) -> getPosition().Y - screwOrigin.Y,2.) >
                    pow(screwCasingRadius - .95*(particleHandler.getObject(i) -> getRadius()),2.)) particleHandler.removeObject(i);
            }
        }
        else    // the particles are eliminate according to the desired final packing fraction
        {
            screwFreeVolume = constants::pi*(screwLength*pow(screwCasingRadius,2.0) - screwThickness*sqrt(1.0 + pow(screwLength/(constants::pi*(screwBladeRadius + screwShaftRadius)), 2))*pow(screwBladeRadius,2.0) - (screwLength - screwThickness)*pow(screwShaftRadius,2.0));
            radialDistanceCuttingLength = getXMax();
            
            double totalVolumeOfParticles;
            
            // loop to find the threshold radius for cutting the particles
            do
            {
                // initializes the volume of the particles to 0
                totalVolumeOfParticles = 0.0;
                
                // loops over all the particles and computes the volume of the ones inside concentric shells
                for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
                {
                    if (pow(particleHandler.getObject(i) -> getPosition().X, 2.) + pow(particleHandler.getObject(i) -> getPosition().Y, 2.) < pow(radialDistanceCuttingLength,2.0)) totalVolumeOfParticles += pow(particleHandler.getObject(i) -> getRadius(), 3);
                }
                totalVolumeOfParticles *= 4.0*constants::pi/3.0;
                
                // updates the radial cutting length every time the packing fraction is higher than the desired one
                radialDistanceCuttingLength -= 0.1*(1.0 - particleDispersity)*particleRadius;
                
            } while (totalVolumeOfParticles > packingFraction*screwFreeVolume);

            // actually cuts the particles outside the cutting length
            for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
            {
                if (pow(particleHandler.getObject(i) -> getPosition().X - screwOrigin.X,2.) +
                    pow(particleHandler.getObject(i) -> getPosition().Y - screwOrigin.Y,2.) >
                    pow(radialDistanceCuttingLength,2.)) particleHandler.removeObject(i);
            }
        }

        casing.setSpecies(specieWall);
        casing.setPosition(screwOrigin);
        casing.setOrientation(Vec3D(0.0,0.0,1.0));
        casing.addObject(Vec3D(1.0,0.0,0.0),Vec3D(screwCasingRadius,0.0,0.0));
        casing.setAngularVelocity(Vec3D(0.,0.,0.));
        wallHandler.copyAndAddObject(casing);
    }
    
    // cuts the particles according to the filling ratio
    void makeFillingRatio()
    {
        for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            if (particleHandler.getObject(i) -> getPosition().Y - screwOrigin.Y > (2.0*fillingRatio-1.0)*screwCasingRadius) particleHandler.removeObject(i);
        }
    }
    
    // rotates the screw
    void makeRotation()
    {
        // the actual rotation of the blade
        helicoidPointer -> rotate(getTimeStep());
        
        // applies the proper angular velocity to the screw blade (for a correct collision computation)
        // IMPORTANT: the sign of the angular velocity depends on the rotation verse as well as if the screw is right/left handed
        helicoidPointer -> setAngularVelocity(Vec3D(0.,0.,-screwAngularVelocity));
        helicoidPointer -> setOrientation(Vec3D(0.0,0.0,1.0));
    }
    
    // increases the velocity
    void makeVelocityIncrement()
    {
        // steps to the next target screw velocity
        if (controlVariableForVelocityIncrement < numberOfvelocityIncrements && getTime() > durationOfVelocityIncrement*controlVariableForVelocityIncrement)
        {
            screwAngularVelocity = velocitySteps[controlVariableForVelocityIncrement];
            helicoidPointer -> setOmega(screwAngularVelocity);
            std::cout << "\nNew screw velocity: " << screwAngularVelocity << "\n";
            controlVariableForVelocityIncrement++;
        }
    }
    
    // creates the data array for volumetric throughput computation and allocates the particles info into it
    void makeParticleList()
    {
        particleList = new double[particleHandler.getNumberOfObjects()];
        
        for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--) {particleList[i] = particleHandler.getObject(i) -> getPosition().Z;}
    }
    
    // refreshes the particle info inside the particleList
    void refreshParticleList()
    {
        for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--) {particleList[i] = particleHandler.getObject(i) -> getPosition().Z;}
    }
    
    // performs the data analysis needed for teh problem
    void makeDataAnalysis()
    {
        // evaluation of mean particle position, mean particle velocity and volumetric throughput
        meanParticlePosition.setZero();
        meanParticleVelocity.setZero();
        volumetricThroughput = 0.0;
        meanCoordinationNumber = 0.0;
        
        for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            // particle position (r, /phi, z)
            meanParticlePosition.X += sqrt(pow(particleHandler.getObject(i) -> getPosition().X,2.0) + pow(particleHandler.getObject(i) -> getPosition().Y,2.0));
            meanParticlePosition.Y += atan2(particleHandler.getObject(i) -> getPosition().Y, particleHandler.getObject(i) -> getPosition().X);
            if (meanParticlePosition.Y < 0.0) meanParticlePosition.Y += 2.0*constants::pi;
            meanParticlePosition.Z += particleHandler.getObject(i) -> getPosition().Z;
            
            // particle velocity d(r, /phi, z)/dt
            meanParticleVelocity.X += ((particleHandler.getObject(i) -> getPosition().X)*(particleHandler.getObject(i) -> getVelocity().X) + (particleHandler.getObject(i) -> getPosition().Y)*(particleHandler.getObject(i) -> getVelocity().Y))/meanParticlePosition.X;
            meanParticleVelocity.Y += ((particleHandler.getObject(i) -> getPosition().Y)*(particleHandler.getObject(i) -> getVelocity().X) - (particleHandler.getObject(i) -> getPosition().X)*(particleHandler.getObject(i) -> getVelocity().Y))/pow(meanParticlePosition.X,2.0);
            meanParticleVelocity.Z += particleHandler.getObject(i) -> getVelocity().Z;
            
            // volumetric throughput
            if (fabs(particleHandler.getObject(i) -> getPosition().Z - particleList[i]) > 0.5*screwLength)
            {
                (particleHandler.getObject(i) -> getPosition().Z - particleList[i]) < 0.0 ? (volumetricThroughput += pow(particleHandler.getObject(i) -> getRadius(),3.0)) : (volumetricThroughput -= pow(particleHandler.getObject(i) -> getRadius(),3.0));
            }
            
            // coordination number
            meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
        }
        
        meanParticlePosition /= particleHandler.getNumberOfObjects();
        meanParticleVelocity /= particleHandler.getNumberOfObjects();
        volumetricThroughput *= 4.0*constants::pi/(3.0*cdatOutputTimeInterval);
        meanCoordinationNumber /= particleHandler.getNumberOfObjects();
        
        // evaluation of the total torque acting on teh screw and of the mean/max particle relative overlap
        meanRelativeOverlap = 0.0;
        maxRelativeOverlap = 0.0;
        totalTorque = 0.0;
        int interactionCounter = 0;
        Vec3D axialPoint;
        axialPoint.X = 0.0;
        axialPoint.Y = 0.0;
        Vec3D translatedForce;
        Vec3D translatedContactPoint;
        
        for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
        {
            // mean and maximum overlap computation
            meanRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (meanRelativeOverlap > maxRelativeOverlap) maxRelativeOverlap = meanRelativeOverlap;
            
            // torque computation
            if ((*i) -> getI() -> getIndex() == helicoidPointer -> getIndex())
            {
                axialPoint.Z = ((*i) -> getContactPoint()).Z;
                
                translatedContactPoint = ((*i) -> getContactPoint()) - axialPoint;
                translatedForce = ((*i) -> getForce()) - axialPoint;
                
                totalTorque += Vec3D::cross(translatedContactPoint, translatedForce).getLength();
            }
            
            interactionCounter++;
        }
        
        meanRelativeOverlap /= interactionCounter;
    }
    
// additional output data
    // creates the data output file and writes the first row
    void makeOutputFile()
    {
        std::ostringstream cdatName;
        std::cout.unsetf(std::ios::floatfield);
        cdatName << getName() << ".cdat";
        
        outputFile.open(cdatName.str(), std::ios::out);
        outputFile << "1. time \t 2. screw_angle \t 3. omega \t 4. Rmean \t 5. PhiMean \t 6. zMean \t 7. vRmean \t 8. vPhiMean \t 9. vZmean \t 10. dV/dt \t 11. torque \t 12. meanCoordinationNumber \t 13. meanOverlap \t 14. maxOverlap " << std::endl;
    }
    
    // writes the compression data to the output file
    void writeDataToOutptFile()
    {
        outputFile <<
        getTime() << "   " <<
        helicoidPointer -> getOffset() << "   " <<
        screwAngularVelocity << "   " <<
        meanParticlePosition.X << "   " <<
        meanParticlePosition.Y << "   " <<
        meanParticlePosition.Z << "   " <<
        meanParticleVelocity.X << "   " <<
        meanParticleVelocity.Y << "   " <<
        meanParticleVelocity.Z << "   " <<
        volumetricThroughput << "   " <<
        totalTorque << "   " <<
        meanCoordinationNumber << "   " <<
        meanRelativeOverlap << "   " <<
        maxRelativeOverlap << "   " <<
        std::endl;
    }
    
    
//  ----- GLOBAL FUNCTIONS -----
    void printTime() const override
    {
        std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() << ", Eratio = " << std::setprecision(6) << std::left << std::setw(10) << getKineticEnergy()/getElasticEnergy() << std::endl;
        std::cout.flush();
    }
    

//  ----- CLASS VARIABLES -----
    // regions and boundaries parameters
    double fillingBoxVolume;
    bool fixedBoxToScrewWidthRatio;
    double boxToScrewWidthRatio;
    double interParticleLoadingLatticeRelativeDistance;
    double interParticleLoadingDistance;
    Vec3D loadingBoxMin;
    Vec3D loadingBoxMax;
    Vec3D fillingBoxMin;
    Vec3D fillingBoxMax;
    int loadingLatticeGridSizeX, loadingLatticeGridSizeY, loadingLatticeGridSizeZ;
    double loadingLatticeSitesDistanceHorizontal;
    double loadingLatticeSitesDistanceVertical;
    double interLatticeSmallParticlesDistance;
    Vec3D loadingLatticeOffset;
    double radialDistanceCuttingLength;
    
    // particle properties
    int numberOfParticles;
    double particleDensity;
    double particleRadius;
    double particleDispersity;
    double particleMass;
    double particleVolume;
    
    // interaction parameters
    double particleStiffness, wallStiffness;
    double particleParticleSlidingFriction;
    double particleParticleRollingFriction;
    double particleParticleTorsionFriction;
    double particleWallSlidingFriction;
    double particleWallRollingFriction;
    double particleWallTorsionFriction;
    double particleParticleRestitutionCoefficient;
    double particleWallRestitutionCoefficient;
    
    // screw properties
    double screwLength;
    double screwCasingRadius;
    double screwBladeRadius;
    double screwShaftRadius;
    double screwNumberOfTurns;
    double screwThickness;
    double screwFreeVolume;
    Vec3D screwOrigin;
    
    // screw operating parameters
    double screwAngularVelocity;
    double *velocitySteps;
    double fillingRatio;
    double packingFraction;
    double screwRunningTime;
    int numberOfvelocityIncrements;
    double durationOfVelocityIncrement;
    
    // data analysis variables
    double cdatOutputTimeInterval;
    Vec3D meanParticlePosition;
    Vec3D meanParticleVelocity;
    double volumetricThroughput;
    double totalTorque;
    double meanCoordinationNumber;
    double meanRelativeOverlap;
    double maxRelativeOverlap;
    double *particleList;
    
    // global paramenters
    int stage;
    bool verbose;
    bool fixedPackingFraction;
    int controlVariableForVelocityIncrement = 0;
    std::ofstream outputFile;
    
    // Mercury specific pointers
    BaseParticle p0;
    InfiniteWall wall;
    PeriodicBoundary boundary;
    AxisymmetricIntersectionOfWalls casing;
    Helicoid05 helicoid;
    Helicoid05 *helicoidPointer;
    LinearViscoelasticFrictionSpecies *specieParticle, *specieWall, *specieMixedParticleWall;
    
    // development hack variables
    bool particleHack;
    int particleHackLayers;
};


int main(int argc, char *argv[])
{
    // the particle diameter, needed to set both particle size and screw geometry
    double particleDiameter = 0.001;
    
    // the velocity steps array
    const int nVelocitySteps = 5;
    double velocityArray [nVelocitySteps] = {0.5*constants::pi, constants::pi, 2.0*constants::pi, 5.0*constants::pi, 10.0*constants::pi};
    
    // Problem setup
    ScrewAlgorithmComparison problem;
    problem.setGravity(Vec3D(0.,-9.81,0.));
    problem.setSystemDimensions(3);
    
    // time scales
    problem.setTimeStep(5.0e-5);
    problem.setTimeMax(30.0);
    problem.setSaveCount(0.01/problem.getTimeStep());
    problem.cdatOutputTimeInterval = 0.005;
    
    // sets particle PSD and density
    problem.setParticleRadiusAndDispersity(5*0.5*particleDiameter, 0.1);
    problem.setParticleDensity(2000.);
    
    // sets the stiffnesses
    problem.setStiffnesses(1000.0, 1000.0);
    
    // sets interaction parameters
    problem.setRestitutionCoefficients(0.8, 0.8);    // p-p and p-w restitution coefficients
    problem.setSlidingFrictionCoefficients(0.3, 0.3);   // p-p and p-w sliding friction coefficients
    problem.setRollingFrictionCoefficients(0.05, 0.05);   // p-p and p-w rolling friction coefficients
    problem.setTorsionFrictionCoefficients(0.0, 0.0);   // p-p and p-w torsion friction coefficients
    
    // length, casing radius, blade radius, shaft radius, number of turns, thickness
    problem.setScrewGeometry(27.0*particleDiameter, 13.5*particleDiameter, 10.0*particleDiameter, 4.0*particleDiameter, 1.0, particleDiameter);
    problem.setScrewOrigin(Vec3D(0.,0.,0.));
    
    // screw operational parameters
    problem.setScrewFillingRatio(0.5);
    problem.setVelocitySteps(nVelocitySteps, velocityArray, 15.0);
    
    // initial loading parameters
    problem.setPackingFraction(false, 0.65);   // if true forces the second value to be the packing fraction after the casing insertion
    problem.setBoxToScrewWidthRatio(true, 1.05); // if true manually sets the box-to-screw size ratio
    problem.setInterParticleLoadingLatticeRelativeDistance(0.1);   // relative to particle radius
    problem.developmentHackParticles(false, 20);
    
    // If true prints the interaction coefficients for every specie
    problem.verbose = true;

    // name setting
    std::ostringstream name;
    std::cout.unsetf(std::ios::floatfield);
//    name << "AlgorithmComparison/0.50/radius_" << 0.5*particleDiameter << "_fillingLevel_0.50";
    name << "Whatever";
    problem.setName(name.str());
//    name << "screw_" << radiusRatio << "_pf_" << std::fixed << std::setprecision(2) << packingFraction << "_" << std::fixed << std::setprecision(2) << sfPP << "_" << std::fixed << std::setprecision(2) << sfPW << "_0.01_0.01_nI_forValidation";
//    problem.setName(name.str());
//    problem.setName("Folder/screw_algorithmComparison");
    problem.setXBallsAdditionalArguments("-h 800 -p 10 -o 200 -3dturn 1");
    
//    problem.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
//    problem.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    
    problem.solve(argc, argv);
    
    return 0;
}


























