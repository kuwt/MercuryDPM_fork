/*
    *** COMPRESSION TEST - SETTLING ***
 Particles are placed inside a cylindrical casing and let settle under gravity.
 When the system is at rest the restart file is created and teh program quits.
*/

#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <math.h>
#include <fstream>

/*
 TODOs:
 
*/
 
class CompressionTest_Settling : public Mercury3D
{
private:
    
    void setupInitialConditions()
    {
        stage = 1;
        setsInserted = 0;
        
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
    
    void actionsAfterTimeStep()
    {
        // particle insertion loop
        if (stage == 2)
        {
            // other rounds of particle insertion
            if (setsInserted < nSets + 1 && getTime() > 0.05*setsInserted)
            {
                makeParticleSet();
            }
             // when all teh particles are loaded, cut the ones exceeding the max height
            if (setsInserted == nSets + 1 && getTime() > 0.05*setsInserted)
            {
                std::cout << "Particle insertion terminated.\nRemoving the exceeding particles...\n";
                removeParticles();
                std::cout << "DONE.\nSettling the particles...";
                stage++;
            }
        }
    }
    
    void actionAfterSolve()
    {
        
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
        auto specieBigWall = speciesHandler.getMixedObject(specieBig, specieWall);
        specieBigWall -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessBig + wallStiffness), bigWallRestitutionCoeff, particleMassBig);
        
        specieBigWall -> setSlidingFrictionCoefficient(bigWallSlidingFrictionCoeff);
        specieBigWall -> setSlidingStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
        specieBigWall -> setSlidingDissipation(specieBigWall -> getDissipation()*2.0/7.0);
        
        specieBigWall -> setRollingFrictionCoefficient(bigWallRollingFrictionCoeff);
        specieBigWall -> setRollingStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
        specieBigWall -> setRollingDissipation(specieBigWall -> getDissipation()*2.0/7.0);
        
        specieBigWall -> setTorsionFrictionCoefficient(bigWallTorsionFrictionCoeff);
        specieBigWall -> setTorsionStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
        specieBigWall -> setTorsionDissipation(specieBigWall -> getDissipation()*2.0/7.0);
        
        // SMALL-WALL
        auto specieSmallWall = speciesHandler.getMixedObject(specieSmall, specieWall);
        specieSmallWall -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessSmall + wallStiffness), smallWallRestitutionCoeff, particleMassSmall);
        
        specieSmallWall -> setSlidingFrictionCoefficient(smallWallSlidingFrictionCoeff);
        specieSmallWall -> setSlidingStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
        specieSmallWall -> setSlidingDissipation(specieSmallWall -> getDissipation()*2.0/7.0);
        
        specieSmallWall -> setRollingFrictionCoefficient(smallWallRollingFrictionCoeff);
        specieSmallWall -> setRollingStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
        specieSmallWall -> setRollingDissipation(specieSmallWall -> getDissipation()*2.0/7.0);
        
        specieSmallWall -> setTorsionFrictionCoefficient(smallWallTorsionFrictionCoeff);
        specieSmallWall -> setTorsionStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
        specieSmallWall -> setTorsionDissipation(specieSmallWall -> getDissipation()*2.0/7.0);
        
        // BIG-SMALL
        auto specieBigSmall = speciesHandler.getMixedObject(specieBig, specieSmall);
        specieBigSmall -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessBig + particleStiffnessSmall), bigSmallRestitutionCoeff, 0.5*(particleMassBig + particleMassSmall));
        
        specieBigSmall -> setSlidingFrictionCoefficient(bigSmallSlidingFrictionCoeff);
        specieBigSmall -> setSlidingStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
        specieBigSmall -> setSlidingDissipation(specieBigSmall -> getDissipation()*2.0/7.0);
        
        specieBigSmall -> setRollingFrictionCoefficient(bigSmallRollingFrictionCoeff);
        specieBigSmall -> setRollingStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
        specieBigSmall -> setRollingDissipation(specieBigSmall -> getDissipation()*2.0/7.0);
        
        specieBigSmall -> setTorsionFrictionCoefficient(bigSmallTorsionFrictionCoeff);
        specieBigSmall -> setTorsionStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
        specieBigSmall -> setTorsionDissipation(specieBigSmall -> getDissipation()*2.0/7.0);
        
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
                
                std::cout << "\tBIG-WALL stiffness and dissipation coefficients: " << specieBigWall -> getStiffness() << " " << specieBigWall -> getDissipation() << "\n";
                std::cout << "\tBIG-WALL friction coefficients: " << bigWallSlidingFrictionCoeff << " " << bigWallRollingFrictionCoeff << " " << bigWallTorsionFrictionCoeff << "\n";
                std::cout << "\tBIG-WALL tangential stiffnesses: " << specieBigWall -> getSlidingStiffness() << " " << specieBigWall -> getRollingStiffness() << " " << specieBigWall -> getTorsionStiffness() << "\n";
                std::cout << "\tBIG-WALL tangential dissipation coefficients: " << specieBigWall -> getSlidingDissipation() << " " << specieBigWall -> getRollingDissipation() << " " << specieBigWall -> getTorsionDissipation() << "\n";
                std::cout << "\tBIG-WALL collision time: " << std::setprecision(4) << specieBigWall -> getCollisionTime(particleMassBig) << "\n\n";
                
                std::cout << "\tSMALL-WALL stiffness and dissipation coefficients: " << specieSmallWall -> getStiffness() << " " << specieSmallWall -> getDissipation() << "\n";
                std::cout << "\tSMALL-WALL friction coefficients: " << smallWallSlidingFrictionCoeff << " " << smallWallRollingFrictionCoeff << " " << smallWallTorsionFrictionCoeff << "\n";
                std::cout << "\tSMALL-WALL tangential stiffnesses: " << specieSmallWall -> getSlidingStiffness() << " " << specieSmallWall -> getRollingStiffness() << " " << specieSmallWall -> getTorsionStiffness() << "\n";
                std::cout << "\tSMALL-WALL tangential dissipation coefficients: " << specieSmallWall -> getSlidingDissipation() << " " << specieSmallWall -> getRollingDissipation() << " " << specieSmallWall -> getTorsionDissipation() << "\n";
                std::cout << "\tSMALL-WALL collision time: " << std::setprecision(4) << specieSmallWall -> getCollisionTime(particleMassSmall) << "\n\n";
                
                std::cout << "\tBIG-SMALL stiffness and dissipation coefficients: " << specieBigSmall -> getStiffness() << " " << specieBigSmall -> getDissipation() << "\n";
                std::cout << "\tBIG-SMALL friction coefficients: " << bigSmallSlidingFrictionCoeff << " " << bigSmallRollingFrictionCoeff << " " << bigSmallTorsionFrictionCoeff << "\n";
                std::cout << "\tBIG-SMALL tangential stiffnesses: " << specieBigSmall -> getSlidingStiffness() << " " << specieBigSmall -> getRollingStiffness() << " " << specieBigSmall -> getTorsionStiffness() << "\n";
                std::cout << "\tBIG-SMALL tangential dissipation coefficients: " << specieBigSmall -> getSlidingDissipation() << " " << specieBigSmall -> getRollingDissipation() << " " << specieBigSmall -> getTorsionDissipation() << "\n";
                std::cout << "\tBIG-SMALL collision time: " << std::setprecision(4) << specieBigSmall -> getCollisionTime(0.5*(particleMassBig + particleMassSmall)) << "\n\n";
            }
            else
            {
                std::cout << "\tPARTICLE stiffness and dissipation coefficients: " << specieBig -> getStiffness() << " " << specieBig -> getDissipation() << "\n";
                std::cout << "\tPARTICLE friction coefficients: " << bigBigSlidingFrictionCoeff << " " << bigBigRollingFrictionCoeff << " " << bigBigTorsionFrictionCoeff << "\n";
                std::cout << "\tPARTICLE tangential stiffnesses: " << specieBig -> getSlidingStiffness() << " " << specieBig -> getRollingStiffness() << " " << specieBig -> getTorsionStiffness() << "\n";
                std::cout << "\tPARTICLE tangential dissipation coefficients: " << specieBig -> getSlidingDissipation() << " " << specieBig -> getRollingDissipation() << " " << specieBig -> getTorsionDissipation() << "\n";
                std::cout << "\tPARTICLE collision time: " << std::setprecision(4) << specieBig -> getCollisionTime(particleMassBig) << "\n\n";
                
                std::cout << "\tPARTICLE-WALL stiffness and dissipation coefficients: " << specieBigWall -> getStiffness() << " " << specieBigWall -> getDissipation() << "\n";
                std::cout << "\tPARTICLE-WALL friction coefficients: " << bigWallSlidingFrictionCoeff << " " << bigWallRollingFrictionCoeff << " " << bigWallTorsionFrictionCoeff << "\n";
                std::cout << "\tPARTICLE-WALL tangential stiffnesses: " << specieBigWall -> getSlidingStiffness() << " " << specieBigWall -> getRollingStiffness() << " " << specieBigWall -> getTorsionStiffness() << "\n";
                std::cout << "\tPARTICLE-WALL tangential dissipation coefficients: " << specieBigWall -> getSlidingDissipation() << " " << specieBigWall -> getRollingDissipation() << " " << specieBigWall -> getTorsionDissipation() << "\n";
                std::cout << "\tPARTICLE-WALL collision time: " << std::setprecision(4) << specieBigWall -> getCollisionTime(particleMassBig) << "\n\n";
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
    
    // removes the particles laying above the casing height
    void removeParticles()
    {
        for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            if (particleHandler.getObject(i) -> getPosition().Z + particleHandler.getObject(i) -> getRadius() > casingHeight) particleHandler.removeObject(i);
        }
    }

 
//  ----- GLOBAL FUNCTIONS -----
    void printTime() const override
    {
        std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() << ", Eratio = " << std::setprecision(6) << std::left << std::setw(10) << getKineticEnergy()/getElasticEnergy() << std::endl;
        std::cout.flush();
    }
    
    bool continueSolve() const override
    {
        if (stage == 3 && getKineticEnergy()<1.0e-5*getElasticEnergy()) return false;
        return true;
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
    
    // simulation variables
    int stage;
    int nSets, nParticlesPerSet;
    double nBigPerSet, nSmallPerSet;
    double nBigParticles, nSmallParticles;
    double smallToBigMassRatio;
    double bulkPackingFractionPreCompression;

    // initialization variables
    double initRegionHeight;
    double volumeOfParticlesPerSet;
    int setsInserted;
    
    // Mercury-specific variables
    LinearViscoelasticFrictionSpecies *specieBig, *specieSmall, *specieWall;
    LinearViscoelasticFrictionSpecies *specieBigWall, *specieSmallWall, *specieBigSmall;
    InfiniteWall basis, roof;
    AxisymmetricIntersectionOfWalls casing;
    BaseParticle p0;
    
    // global variables
    bool verbose;
    bool biModal;
};

int main(int argc, char *argv[])
{
    double scalingFactor = 2.0; // the scaling factor for the particle sizes and cutoffs
    
    CompressionTest_Settling ctSettle;
    ctSettle.setName("CompressionTest_Settling");
    
    // sets simulation parameters
    ctSettle.setTimeStep(1.0e-5);
    ctSettle.setTimeMax(30.0);
    ctSettle.setGravity(Vec3D(0., 0., -9.81));
    ctSettle.setSystemDimensions(3);
    ctSettle.setVerbose(true);
    
    // sets the number of saved time steps such that the output is printed every 0.01s
    ctSettle.setSaveCount(0.01/ctSettle.getTimeStep());
    
    // sets the particle intrinsic properties
    ctSettle.setParticleDensity(2000.0, 2000.0);
    ctSettle.setParticleRadiusAndDispersity(scalingFactor*0.002, 0.10, scalingFactor*0.001, 0.10);
    
    // sets the small-to-big total mass ratio
    ctSettle.setTotalMassRatio(0.1);
    
    // sets the stiffnesses
    ctSettle.setWallStiffness(2000.0);
    ctSettle.setParticleStiffness(1000.0);
    
    // sets the particle-wall friction coefficients
    ctSettle.setParticleWallSlidingFrictionCoeff(0.3, 0.3);
    ctSettle.setParticleWallRollingFrictionCoeff(0.01, 0.01);
    ctSettle.setParticleWallTorsionFrictionCoeff(0.0, 0.0);
    
    // sets the particle-particle friction coefficients
    ctSettle.setParticleParticleSlidingFrictionCoeff(0.5, 0.8, 0.7);
    ctSettle.setParticleParticleRollingFrictionCoeff(0.05, 0.1, 0.75);
    ctSettle.setParticleParticleTorsionFrictionCoeff(0.0, 0.0, 0.0);
    
    // sets the particle restitution coefficients
    ctSettle.setParticleWallRestitutionCoeff(0.8, 0.3);
    ctSettle.setParticleParticleRestitutionCoeff(0.5, 0.2, 0.3);
    
    // sets the casing dimensions
    ctSettle.setCasingDimensions(0.04, 0.1);
    
    // sets the bulk packing fraction of the powder prior to compression
    ctSettle.setbulkPackingFractionPreCompression(0.70);
    
    // sets additional built-in arguments for the xballs visualization
    ctSettle.setXBallsAdditionalArguments("-h 800 -p 10 -o 200 -3dturn 3");
    
    ctSettle.solve();
    
    return 0;
}







