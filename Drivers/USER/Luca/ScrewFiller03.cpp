#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Walls/InfiniteWall.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include "Helicoid03.h"

/*
 Last update : 19.11.16
 Development stage : ALPHA
 
 ToDo:
 MAJOR
  - fix the restart file
 - provide separate entities for big and small particles!!! The interaction with the big particles is wrong otherwise! (Wrong interaction radius).
 
 MINOR

 
 WHATEVER
 - explain that the particles needed and teh particles actually loaded are not the same: p_loaded >= p_needed
 - implement a parameter study cycle in the main
 
 */

class ScrewFiller03 : public Mercury3D
{
private:
    
    void setupInitialConditions()
    {
        stage = 1;
        
        // computes the mean mass and volume of the particles
        std::cout << "Computing particle mean mass and volume... ";
        setParticleMassAndVolume();
        std::cout << "DONE!\n";
        
        // sets the species of the pbjects
        std::cout << "Setting up species... ";
        setSpecies();
        std::cout << "DONE!\n";
        
        // sets the box around the screw
        std::cout << "Setting up the filling box dimensions... ";
        setFillingBoxDimensions();
        std::cout << "\nFilling box dimensions: (" << fillingBoxMax.X - fillingBoxMin.X << " , " << fillingBoxMax.Y - fillingBoxMin.Y << " , " << fillingBoxMax.Z - fillingBoxMin.Z << ").\n";
        std::cout << "Filling box volume: " << fillingBoxVolume << ".\n";
        
        // sets the number of particles needed
        std::cout << "Setting up the number of particles needed (accounts for the shaft volume)... ";
        setNumberOfParticles();
        if (!biModalDispersity)
        {
            std::cout << "\nMono-modal particle size distribution. Polydispersity ratio: " << particleDispersity << ".";
            std::cout << "\nMean particle radius: " << particleRadius << ".";
            std::cout << "\nNumber of particles needed: " << numberOfParticles << ".\n";
        }
        else
        {
            std::cout << "\nBi-modal particle size distribution. Polydispersity ratio: " << particleDispersity << ".";
            std::cout << "\nMean small particle radius: " << particleRadius << ", mean big particle radius: " << 2.0*particleRadius << ".";
            std::cout << "\nTotal number of particles needed: " << numberOfParticles << ".";
            std::cout << "\nNumber of small particles needed: " << numberOfParticles - numberOfCubicLattices << ", number of big particles needed: " << numberOfCubicLattices << ".\n";
        }
        
        // sets the lattice properties where the particles will be initially loaded
        std::cout << "Setting up the loading lattice properties... ";
        setloadingBoxGridProperties();
        std::cout << "DONE!\n";
        
        // sets the size of the initial particle loading box on top of the filling box
        std::cout << "Setting up the loading box dimensions... ";
        setLoadingBoxDimensions();
        std::cout << "\nLoading box dimensions: (" << loadingBoxMax.X - loadingBoxMin.X << " , " << loadingBoxMax.Y - loadingBoxMin.Y << " , " << loadingBoxMax.Z - loadingBoxMin.Z << ").\n";
        std::cout << "Loading lattice grid size: (" << loadingLatticeGridSizeX << " , " << loadingLatticeGridSizeY << " , " << loadingLatticeGridSizeZ << ").\n";
        std::cout << "Loading lattice sites distance: (" << loadingLatticeSitesDistanceHorizontal << " , " << loadingLatticeSitesDistanceVertical << " , " << loadingLatticeSitesDistanceHorizontal << ").\n";
        std::cout << "Loading lattice offsets: (" << loadingLatticeOffset.X << " , " << loadingLatticeOffset.Y << " , " << loadingLatticeOffset.Z << ").\n";
        
        // sets the simulation domain
        std::cout << "Setting up the simulation boundaries... ";
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
            std::cout << "Starting the screw rotation...\n";
            std::cout << "Resetting the maximum time...\n";
            setTimeMax(screwRunningTime + getTime());
            
            stage++;
        }
        
        // stage 4: the screw is set in motion
        if (stage == 4)
        {
            makeRotation();
        }
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
    
    // sets the bi-modal dispersity true or false
    void setParticleBiModalDispersity(bool biDisp)
    {
        biModalDispersity = biDisp;
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
    
    // sets the scew velocity
    void setScrewVelocity(double omega)
    {
        screwAngularVelocity = omega;
    }
    
    // sets the scew filling ratio
    void setScrewFillingRatio(double fR)
    {
        fillingRatio = fR;
    }
    
    // sets the scew running time
    void setScrewRunningTime(double time)
    {
        screwRunningTime = time;
    }
    
    
// interaction parameters
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
    // sets the collision time
    void setCollisionTime(double ct)
    {
        collisionTime = ct;
    }
    
    // sets the ratio between the filling box size and the screw casing diameter
    void setBoxToScrewWidthRatio(double ratio)
    {
        boxToScrewWidthRatio = ratio;
    }
    
    // sets the interparticle distance inside of the loading lattice normalized by the mean radius
    void setInterParticleLoadingLatticeRelativeDistance(double relDist)
    {
        interParticleLoadingLatticeRelativeDistance = relDist;
        interParticleLoadingDistance = interParticleLoadingLatticeRelativeDistance*particleRadius*(1.0 + particleDispersity);
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
        fillingBoxMin.X = screwOrigin.X - boxToScrewWidthRatio*screwCasingRadius;
        fillingBoxMin.Y = screwOrigin.Y - boxToScrewWidthRatio*screwCasingRadius;
        fillingBoxMin.Z = screwOrigin.Z;
        
        fillingBoxMax.X = screwOrigin.X + boxToScrewWidthRatio*screwCasingRadius;
        fillingBoxMax.Y = screwOrigin.Y + boxToScrewWidthRatio*screwCasingRadius;
        fillingBoxMax.Z = screwOrigin.Z + screwLength;
        
        fillingBoxVolume = (fillingBoxMax.X - fillingBoxMin.X)*(fillingBoxMax.Y - fillingBoxMin.Y)*(fillingBoxMax.Z - fillingBoxMin.Z);
    }
    
    // sets the loading box parameters
    void setloadingBoxGridProperties()      // explain how this thing works
    {
        if (!biModalDispersity) // mono-modal
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
        else    // bi-modal
        {
            interLatticeSmallParticlesDistance = 2.0*(3.0*particleRadius*(1.0 + particleDispersity) + interParticleLoadingDistance)/sqrt(3.0);
            double cubicLatticeSize = interLatticeSmallParticlesDistance + 2.0*particleRadius*(1.0 + particleDispersity) + interParticleLoadingDistance;
            
            loadingLatticeSitesDistanceHorizontal = cubicLatticeSize;
            loadingLatticeSitesDistanceVertical = cubicLatticeSize;
            
            loadingLatticeGridSizeX = (int)((fillingBoxMax.X - fillingBoxMin.X)/cubicLatticeSize);
            loadingLatticeGridSizeZ = (int)((fillingBoxMax.Z - fillingBoxMin.Z)/cubicLatticeSize);
            
            loadingLatticeGridSizeY = (int)(numberOfCubicLattices/(loadingLatticeGridSizeX*loadingLatticeGridSizeZ)) + 1;
            
            loadingLatticeOffset.X = 0.5*(fillingBoxMax.X - fillingBoxMin.X - cubicLatticeSize*loadingLatticeGridSizeX);
            loadingLatticeOffset.Y = 0.0;
            loadingLatticeOffset.Z = 0.5*(fillingBoxMax.Z - fillingBoxMin.Z - cubicLatticeSize*loadingLatticeGridSizeZ);
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
        specieParticle -> setCollisionTimeAndRestitutionCoefficient(collisionTime, particleParticleRestitutionCoefficient, particleMass);
        
        specieParticle -> setSlidingDissipation(specieParticle -> getDissipation()*2.0/7.0);
        specieParticle -> setSlidingStiffness(specieParticle -> getStiffness()*2.0/7.0);
        specieParticle -> setSlidingFrictionCoefficient(particleParticleSlidingFriction);
        
        specieParticle -> setRollingDissipation(specieParticle -> getDissipation()*2.0/7.0);
        specieParticle -> setRollingStiffness(specieParticle -> getStiffness()*2.0/7.0);
        specieParticle -> setRollingFrictionCoefficient(particleParticleRollingFriction);
        
        specieParticle -> setTorsionDissipation(specieParticle -> getDissipation()*2.0/7.0);
        specieParticle -> setTorsionStiffness(specieParticle -> getStiffness()*2.0/7.0);
        specieParticle -> setTorsionFrictionCoefficient(particleParticleTorsionFriction);
        speciesHandler.addObject(specieParticle);
        
        // wall
        specieWall = new LinearViscoelasticFrictionSpecies;
        specieWall -> setDensity(particleDensity);
        specieWall -> setCollisionTimeAndRestitutionCoefficient(collisionTime, particleWallRestitutionCoefficient, particleMass);
        
        specieWall -> setSlidingDissipation(specieWall -> getDissipation()*2.0/7.0);
        specieWall -> setSlidingStiffness(specieWall -> getStiffness()*2.0/7.0);
        specieWall -> setSlidingFrictionCoefficient(particleWallSlidingFriction);
        
        specieWall -> setRollingDissipation(specieWall -> getDissipation()*2.0/7.0);
        specieWall -> setRollingStiffness(specieWall -> getStiffness()*2.0/7.0);
        specieWall -> setRollingFrictionCoefficient(particleWallRollingFriction);
        
        specieWall -> setTorsionDissipation(specieWall -> getDissipation()*2.0/7.0);
        specieWall -> setTorsionStiffness(specieWall -> getStiffness()*2.0/7.0);
        specieWall -> setTorsionFrictionCoefficient(particleWallTorsionFriction);
        speciesHandler.addObject(specieWall);
        
        
        // particle-wall mixed
        auto specieMixedParticleWall = speciesHandler.getMixedObject(specieParticle, specieWall);
        specieMixedParticleWall -> setCollisionTimeAndRestitutionCoefficient(collisionTime, particleWallRestitutionCoefficient, particleMass, particleMass);
        
        specieMixedParticleWall -> setSlidingDissipation(specieMixedParticleWall -> getDissipation()*2.0/7.0);
        specieMixedParticleWall -> setSlidingStiffness(specieMixedParticleWall -> getStiffness()*2.0/7.0);
        specieMixedParticleWall -> setSlidingFrictionCoefficient(particleWallSlidingFriction);
        
        specieMixedParticleWall -> setRollingDissipation(specieMixedParticleWall -> getDissipation()*2.0/7.0);
        specieMixedParticleWall -> setRollingStiffness(specieMixedParticleWall -> getStiffness()*2.0/7.0);
        specieMixedParticleWall -> setRollingFrictionCoefficient(particleWallRollingFriction);
        
        specieMixedParticleWall -> setTorsionDissipation(specieMixedParticleWall -> getDissipation()*2.0/7.0);
        specieMixedParticleWall -> setTorsionStiffness(specieMixedParticleWall -> getStiffness()*2.0/7.0);
        specieMixedParticleWall -> setTorsionFrictionCoefficient(particleWallTorsionFriction);
    }
    
    // sets the number of particles
    void setNumberOfParticles() // safe constant * assumed packing fraction * (box volume - shaft volume) / mean particle volume
    {
        if (!biModalDispersity) // mono-modal
        {
            numberOfParticles = (int)(1.2*0.6*(fillingBoxVolume - constants::pi*screwLength*pow(screwShaftRadius,2.))/particleVolume);
            numberOfCubicLattices = 0;
        }
        else    // bi-modal
        {
            numberOfCubicLattices = (int)(1.2*0.6*(fillingBoxVolume - constants::pi*screwLength*pow(screwShaftRadius,2.))/(16.0*particleVolume));
            numberOfParticles = 9*numberOfCubicLattices;
        }
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
        
        if (!biModalDispersity) // mono-modal lattice
        {
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
        else    // bi-modal lattice
        {
            for (int i=0; i<loadingLatticeGridSizeY; i++)
            {
                yPos = loadingBoxMin.Y + loadingLatticeOffset.Y + (i + 0.5)*loadingLatticeSitesDistanceVertical;
                
                for (int j=0; j<loadingLatticeGridSizeX; j++)
                {
                    xPos = loadingBoxMin.X + loadingLatticeOffset.X + (j + 0.5)*loadingLatticeSitesDistanceHorizontal;
                    for (int k=0; k<loadingLatticeGridSizeZ; k++)
                    {
                        zPos = loadingBoxMin.Z + loadingLatticeOffset.Z + (k + 0.5)*loadingLatticeSitesDistanceHorizontal;
                        
                        // big particle placement
                        p0.setPosition(Vec3D(xPos,yPos,zPos));
                        p0.setRadius(2.0*particleRadius*(1.0 + particleDispersity*random.getRandomNumber(-1.,1.)));
                        particleHandler.copyAndAddObject(p0);
                        
                        // small particles placement along the lattice sites
                        for (int l=0; l<2; l++)
                        {
                            for (int m=0; m<2; m++)
                            {
                                for (int n=0; n<2; n++)
                                {
                                    p0.setPosition(Vec3D(xPos + (l - 0.5)*interLatticeSmallParticlesDistance,yPos + (m - 0.5)*interLatticeSmallParticlesDistance,zPos + (n - 0.5)*interLatticeSmallParticlesDistance));
                                    p0.setRadius(particleRadius*(1.0 + particleDispersity*random.getRandomNumber(-1.,1.)));
                                    particleHandler.copyAndAddObject(p0);
                                }
                            }
                        }
                    }
                }
            }
            
        }
    }
    
    // makes the screw
    void makeScrew()
    {
        shaft = new AxisymmetricIntersectionOfWalls();
        shaft -> setSpecies(specieWall);
        shaft -> setPosition(screwOrigin);
        shaft -> setOrientation(Vec3D(0.,0.,1.));
        shaft -> addObject(Vec3D(-1.,0.,0.), Vec3D(screwShaftRadius,0.,0.));
        shaft -> setAngularVelocity(Vec3D(0.,0.,0.));
        wallHandler.addObject(shaft);
        
//        // Correct screw implementation (although rotation not working)
//        helicoid.setSpecies(specieWall);
//        helicoid.set(screwOrigin, screwLength, screwBladeRadius, screwNumberOfTurns, screwAngularVelocity, screwThickness);
//        wallHandler.copyAndAddObject(helicoid);
        
        // Pointer implementation of the screw (wrong one but does the trick)
        helicoid = wallHandler.copyAndAddObject(Helicoid03());
        helicoid -> set(screwOrigin, screwLength, screwBladeRadius, screwNumberOfTurns, screwAngularVelocity, screwThickness);
        helicoid -> setSpecies(specieWall);
    }
    
    // makes the screw casing
    void makeCasing()
    {
        for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            if (pow(particleHandler.getObject(i) -> getPosition().X - screwOrigin.X,2.) +
                pow(particleHandler.getObject(i) -> getPosition().Y - screwOrigin.Y,2.) >
                pow(screwCasingRadius - .95*(particleHandler.getObject(i) -> getRadius()),2.)) particleHandler.removeObject(i);
        }
        
        casing = new AxisymmetricIntersectionOfWalls();
        casing -> setSpecies(specieWall);
        casing -> setPosition(screwOrigin);
        casing -> setOrientation(Vec3D(0.,0.,1.));
        casing -> addObject(Vec3D(1.,0.,0.), Vec3D(screwCasingRadius,0.,0.));
        casing -> setAngularVelocity(Vec3D(0.,0.,0.));
        wallHandler.addObject(casing);
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
        helicoid -> move_time(getTimeStep());
        shaft -> setAngularVelocity(Vec3D(0.,0.,screwAngularVelocity));
        shaft -> setOrientation(Vec3D(0.0,0.0,1.0));
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
    int numberOfCubicLattices;
    double interLatticeSmallParticlesDistance;
    Vec3D loadingLatticeOffset;
    
    // particle properties
    int numberOfParticles;
    double particleDensity;
    double particleRadius;
    double particleDispersity;
    double particleMass;
    double particleVolume;
    bool biModalDispersity;
    
    // interaction parameters
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
    Vec3D screwOrigin;
    
    // screw operating parameters
    double screwAngularVelocity;
    double fillingRatio;
    double screwRunningTime;
    
    // global paramenters
    int stage;
    double collisionTime;
    
    // Mercury specific pointers
    BaseParticle p0;
    InfiniteWall wall;
    PeriodicBoundary boundary;
    AxisymmetricIntersectionOfWalls *shaft, *casing;
    Helicoid03 *helicoid;
    LinearViscoelasticFrictionSpecies *specieParticle, *specieWall, *specieMixedParticleWall;
    
    // development hack variables
    bool particleHack;
    int particleHackLayers;
};


int main(int argc, char *argv[])
{
    // Problem setup
    ScrewFiller03 problem;
    problem.setName("ScrewFiller03");

    problem.setGravity(Vec3D(0.,-9.81,0.));
    problem.setSystemDimensions(3);
    
    problem.setCollisionTime(0.001);
    problem.setTimeStep(0.00001);
    problem.setTimeMax(20.0);
    
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(1000, problem.getTimeMax(), problem.getTimeStep()));
    
    problem.setParticleRadiusAndDispersity(0.009,0.1);
    problem.setParticleDensity(2000.);
    problem.setSlidingFrictionCoefficients(0.2, 0.5);   // p-p and p-w sliding friction coefficients
    problem.setRollingFrictionCoefficients(0.1, 0.2);   // p-p and p-w rolling friction coefficients
    problem.setTorsionFrictionCoefficients(0.01, 0.1);   // p-p and p-w torsion friction coefficients
    problem.setRestitutionCoefficients(0.8, 0.9);    // p-p and p-w restitution coefficients
    
    problem.setParticleBiModalDispersity(false);    // switches the bi-modal size dispersity on or off
    
    // length, casing radius, blade radius, shaft radius, number of turns, thickness
    problem.setScrewGeometry(.4, .185, .185, .05, 1.0, 0.03);
    problem.setScrewOrigin(Vec3D(0.,0.,0.));
    problem.setScrewVelocity(constants::pi);
    problem.setScrewFillingRatio(.35);
    problem.setScrewRunningTime(20.0);  // this will overwrite the value set by setTimeMax once the screw starts turning
    
    problem.setBoxToScrewWidthRatio(1.5);
    problem.setInterParticleLoadingLatticeRelativeDistance(0.1);   // relative to particle radius
    
    problem.developmentHackParticles(false, 8);

//    problem.setXBallsAdditionalArguments("-solidf -v0");
    problem.setXBallsAdditionalArguments("-h 800 -p 10 -o 200 -3dturn 1");
    
    problem.solve(argc, argv);
    
    return 0;
}


























