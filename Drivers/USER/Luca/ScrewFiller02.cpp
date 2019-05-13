
#include "Mercury3D.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Helicoid03.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "sstream"

/*
 Last update : 6.09.16 at 18:00
 Development stage : ALPHA
 
 ToDo:
 MAJOR
 - correct the setTimeMax stuff, and check if it quits the sim beforetime
 - bug in the particle cutting: wrong number of objects (in the mono-modal case)
    probably the mistake is here 
        nLoadingCycles = (int)(nTotParticles/nParticlesPerLoadingCycle);
    since, according to the rounding, I can take different nLoadingCycles
    which might lead to nLoadingCycles*nParticlesPerLoadingCycle != nTotParticles
    *** I THINK NOW IT'S WORKING ***
 - change the particle trimmer function so it uses the handler.getnumberof objects
 - check if the bi-disperse case works
 
 MINOR
 - check for the estimated time and correct it (the estimated time is far longer than the actual one)
 - write teh checks that exit the execution if parameters make no sense
 - ask if special things have to be accounted for when writing the restart file
 
 WHATEVER
 - create a .h and .cc of this stuff
 - write the comments in doxigen style and standardise them
 - provide a talkative version
 - write a small thing that outputs the spatial distribution of the big/small particles to check if isotropic (cut in cake style)
 
 DEVELOPMENT

 */

class ScrewFiller : public Mercury3D
{
    
public:
//  --- STANDARD FUNCTIONS ---
    
    void setupInitialConditions()
    {
        
        std::cout << "\nWeighting the particles...";
        computeMassSphere();
        std::cout << "\nWeighting DONE!\nDefining species...";
        setupSpecies();
        std::cout << "\nSpecies definition DONE!\nSetting up boundaries...";
        
        wallHandler.clear();
        setupBoundaries();
        std::cout << "\nBoundaries DONE!\nSetting up walls...";
        setupWalls();
        std::cout << "\nWalls DONE!\nSetting up screw...";
        setupScrew();
        std::cout << "\nScrew DONE!\nUsing dirty algebra to compute numbers and times...";
        
        computeParticleNumber();
        computeParticlesPerLoadingCycle();
        std::cout << "\nWe need " << nLoadingCycles*nParticlesPerLoadingCycle << " particles to fill the box.";
        std::cout << "\nWe'll run " << nLoadingCycles << " insertion cycles of " << nParticlesPerLoadingCycle << " particles each.";
        settlingTimeEstimator();
        setTimeMax(getEstimatedTime() + 2.0);
        std::cout << "\nThe settling process should take around " << loadingTotalTime << " seconds of physical time.\n";
        
        // calls the first loading cycle
        if (stage == 0)
        {
            
            std::cout << "\nLoading cycle n. " << loadCyclesCounter + 1 << "\n";
            particleHandler.clear();
            particleLoader();
            loadCyclesCounter++;
            
            stage++;
            
        }
        
    }
    
    void actionsBeforeTimeStep()
    {
        
        // calls the remaining loading cycles
        if (stage == 1 && (int)(getTime()/loadingCycleTime) == loadCyclesCounter /*&& loadCyclesCounter < nLoadingCycles*/)
        {
            
            std::cout << "\nLoading cycle n. " << loadCyclesCounter + 1 << "\n";
            particleLoader();
            loadCyclesCounter++;
            
            if (loadCyclesCounter == nLoadingCycles)
            {
                
                std::cout << "\nParticle loading terminated.\nParticle settle begun...\n";
                stage++;
                
            }
            
        }
        
        // let's check the trimming in another file, to determine the correct cutting condition
        // remember to modify the continueSolve function
        
//        // after the loading time as soon as the particles are still loads the cylinder and cuts the particles outside
//        if (stage == 2 && getKineticEnergy()/getElasticEnergy() < 1.e-4)
//        {
//            
//            std::cout << "\nCutting the particles and loading the casing...\n";
//            particleTrimmer(nLoadingCycles*nParticlesPerLoadingCycle);
//            setupCasing();
//            
//            stage++;
//            
//        }
        
        // after the cut, whenever the particles are still, the restart file is written and the program quits

    }
    
    
//  --- FUNCTIONS CALLED IN THE MAIN ---
    
    // sets the small particle radius and dispersity
    void setParticleRadiusAndDPolyspersity(double radiusSmall, double dispSmall)
    {
        
        particleRadiusSmall = radiusSmall;
        particleRadiusBig = radiusSmall;
        particleDispersitySmall = dispSmall;
        particleDispersityBig = dispSmall;
        biDisperseSystem = false;
        
    }
    
    // sets the small particle radius and dispersity
    void setParticleRadiusAndDPolyspersity(double radiusSmall, double dispSmall, double dispBig)
    {
        
        particleRadiusSmall = radiusSmall;
        particleRadiusBig = 2.0*radiusSmall;
        particleDispersitySmall = dispSmall;
        particleDispersityBig = dispBig;
        biDisperseSystem = true;
        
    }
    
    // sets the density of the particles
    void setParticleDensity(double density)
    {
        
        particleDensity = density;
        
    }
    
    // sets the sliding friction coefficients
    void setSlidingFrictionCoefficients(double ppsf, double pwsf)
    {
        
        particleParticleSlidingFriction = ppsf;
        particleWallSlidingFriction = pwsf;
        
    }
    
    // sets the rolling friction coefficients
    void setRollingFrictionCoefficients(double pprf, double pwrf)
    {
        
        particleParticleRollingFriction = pprf;
        particleWallRollingFriction = pwrf;
        
    }
    
    // sets the torsion friction coefficients
    void setTorsionFrictionCoefficients(double pptf, double pwtf)
    {
        
        particleParticleTorsionFriction = pptf;
        particleWallTorsionFriction = pwtf;
        
    }
    
    // sets the coefficient of restitution
    void setRestitutionCoefficients(double pprc, double pwrc)
    {
        
        particleParticleRestitutionCoefficient = pprc;
        particleWallRestitutionCoefficient = pwrc;
        
    }
    
    // sets the collision time
    void setCollisionTime(double ct)
    {
        collisionTime = ct;
    }
    
    // sets the geometry variables of the screw
    void setScrewGeometry(double length, double casing, double radius, double shaft, double n, double thickness)
    {
        
        screwLength = length;
        screwCasingRadius = casing;
        screwBladeRadius = radius;
        screwShaftRadius = shaft;
        screwNumberOfTurns = n;
        screwThickness = thickness;
        
    }
    
    // sets the initial point of the screw's axis
    void setScrewOrigin(Vec3D origin)
    {
        screwOrigin = origin;
    }
    
    // sets the angular velocity of the screw
    void setScrewVelocity(double omega)
    {
        
        screwAngularVelocity = omega;
        
    }
    

//  --- FUNCTIONS USED BY THE CLASS ---
    
    // sets up the species
    void setupSpecies()
    {
        
        speciesHandler.clear();
        
        // small pure
        speciesSmall = new LinearViscoelasticFrictionSpecies;
        speciesSmall -> setDensity(particleDensity);
        speciesSmall -> setCollisionTimeAndRestitutionCoefficient(collisionTime, particleParticleRestitutionCoefficient, particleMassSmall);
        
        speciesSmall -> setSlidingDissipation(speciesSmall -> getDissipation()*2.0/7.0);
        speciesSmall -> setSlidingStiffness(speciesSmall -> getStiffness()*2.0/7.0);
        speciesSmall -> setSlidingFrictionCoefficient(particleParticleSlidingFriction);
        
        speciesSmall -> setRollingDissipation(speciesSmall -> getDissipation()*2.0/7.0);
        speciesSmall -> setRollingStiffness(speciesSmall -> getStiffness()*2.0/7.0);
        speciesSmall -> setRollingFrictionCoefficient(particleParticleRollingFriction);
        
        speciesSmall -> setTorsionDissipation(speciesSmall -> getDissipation()*2.0/7.0);
        speciesSmall -> setTorsionStiffness(speciesSmall -> getStiffness()*2.0/7.0);
        speciesSmall -> setTorsionFrictionCoefficient(particleParticleTorsionFriction);
        speciesHandler.addObject(speciesSmall);
        
        
        // big pure
        speciesBig = new LinearViscoelasticFrictionSpecies;
        speciesBig -> setDensity(particleDensity);
        speciesBig -> setCollisionTimeAndRestitutionCoefficient(collisionTime, particleParticleRestitutionCoefficient, particleMassBig);
        
        speciesBig -> setSlidingDissipation(speciesBig -> getDissipation()*2.0/7.0);
        speciesBig -> setSlidingStiffness(speciesBig -> getStiffness()*2.0/7.0);
        speciesBig -> setSlidingFrictionCoefficient(particleParticleSlidingFriction);
        
        speciesBig -> setRollingDissipation(speciesBig -> getDissipation()*2.0/7.0);
        speciesBig -> setRollingStiffness(speciesBig -> getStiffness()*2.0/7.0);
        speciesBig -> setRollingFrictionCoefficient(particleParticleRollingFriction);
        
        speciesBig -> setTorsionDissipation(speciesBig -> getDissipation()*2.0/7.0);
        speciesBig -> setTorsionStiffness(speciesBig -> getStiffness()*2.0/7.0);
        speciesBig -> setTorsionFrictionCoefficient(particleParticleTorsionFriction);
        speciesHandler.addObject(speciesBig);
        
        
        // screw pure
        speciesScrew = new LinearViscoelasticFrictionSpecies;
        speciesScrew -> setDensity(particleDensity);
        speciesScrew -> setCollisionTimeAndRestitutionCoefficient(collisionTime, particleWallRestitutionCoefficient, particleMassSmall);
        
        speciesScrew -> setSlidingDissipation(speciesScrew -> getDissipation()*2.0/7.0);
        speciesScrew -> setSlidingStiffness(speciesScrew -> getStiffness()*2.0/7.0);
        speciesScrew -> setSlidingFrictionCoefficient(particleWallSlidingFriction);
        
        speciesScrew -> setRollingDissipation(speciesScrew -> getDissipation()*2.0/7.0);
        speciesScrew -> setRollingStiffness(speciesScrew -> getStiffness()*2.0/7.0);
        speciesScrew -> setRollingFrictionCoefficient(particleWallRollingFriction);
        
        speciesScrew -> setTorsionDissipation(speciesScrew -> getDissipation()*2.0/7.0);
        speciesScrew -> setTorsionStiffness(speciesScrew -> getStiffness()*2.0/7.0);
        speciesScrew -> setTorsionFrictionCoefficient(particleWallTorsionFriction);
        speciesHandler.addObject(speciesScrew);

        
        // small - screw mixed
        auto speciesMixedSmallScrew = speciesHandler.getMixedObject(speciesSmall, speciesScrew);
        speciesMixedSmallScrew -> setCollisionTimeAndRestitutionCoefficient(collisionTime, particleWallRestitutionCoefficient, particleMassSmall, particleMassSmall);
        
        speciesMixedSmallScrew -> setSlidingDissipation(speciesMixedSmallScrew -> getDissipation()*2.0/7.0);
        speciesMixedSmallScrew -> setSlidingStiffness(speciesMixedSmallScrew -> getStiffness()*2.0/7.0);
        speciesMixedSmallScrew -> setSlidingFrictionCoefficient(particleWallSlidingFriction);
        
        speciesMixedSmallScrew -> setRollingDissipation(speciesMixedSmallScrew -> getDissipation()*2.0/7.0);
        speciesMixedSmallScrew -> setRollingStiffness(speciesMixedSmallScrew -> getStiffness()*2.0/7.0);
        speciesMixedSmallScrew -> setRollingFrictionCoefficient(particleWallRollingFriction);
        
        speciesMixedSmallScrew -> setTorsionDissipation(speciesMixedSmallScrew -> getDissipation()*2.0/7.0);
        speciesMixedSmallScrew -> setTorsionStiffness(speciesMixedSmallScrew -> getStiffness()*2.0/7.0);
        speciesMixedSmallScrew -> setTorsionFrictionCoefficient(particleWallTorsionFriction);

        
        // small - big mixed
        auto speciesMixedSmallBig = speciesHandler.getMixedObject(speciesSmall, speciesBig);
        speciesMixedSmallBig -> setCollisionTimeAndRestitutionCoefficient(collisionTime, particleParticleRestitutionCoefficient, particleMassSmall, particleMassBig);
        
        speciesMixedSmallBig -> setSlidingDissipation(speciesMixedSmallBig -> getDissipation()*2.0/7.0);
        speciesMixedSmallBig -> setSlidingStiffness(speciesMixedSmallBig -> getStiffness()*2.0/7.0);
        speciesMixedSmallBig -> setSlidingFrictionCoefficient(particleParticleSlidingFriction);
        
        speciesMixedSmallBig -> setRollingDissipation(speciesMixedSmallBig -> getDissipation()*2.0/7.0);
        speciesMixedSmallBig -> setRollingStiffness(speciesMixedSmallBig -> getStiffness()*2.0/7.0);
        speciesMixedSmallBig -> setRollingFrictionCoefficient(particleParticleRollingFriction);
        
        speciesMixedSmallBig -> setTorsionDissipation(speciesMixedSmallBig -> getDissipation()*2.0/7.0);
        speciesMixedSmallBig -> setTorsionStiffness(speciesMixedSmallBig -> getStiffness()*2.0/7.0);
        speciesMixedSmallBig -> setTorsionFrictionCoefficient(particleParticleTorsionFriction);
        
        
        // big - screw mixed
        auto speciesMixedBigScrew = speciesHandler.getMixedObject(speciesBig, speciesScrew);
        speciesMixedBigScrew -> setCollisionTimeAndRestitutionCoefficient(collisionTime, particleWallRestitutionCoefficient, particleMassBig, particleMassSmall);
        
        speciesMixedBigScrew -> setSlidingDissipation(speciesMixedBigScrew -> getDissipation()*2.0/7.0);
        speciesMixedBigScrew -> setSlidingStiffness(speciesMixedBigScrew -> getStiffness()*2.0/7.0);
        speciesMixedBigScrew -> setSlidingFrictionCoefficient(particleWallSlidingFriction);
        
        speciesMixedBigScrew -> setRollingDissipation(speciesMixedBigScrew -> getDissipation()*2.0/7.0);
        speciesMixedBigScrew -> setRollingStiffness(speciesMixedBigScrew -> getStiffness()*2.0/7.0);
        speciesMixedBigScrew -> setRollingFrictionCoefficient(particleWallRollingFriction);
        
        speciesMixedBigScrew -> setTorsionDissipation(speciesMixedBigScrew -> getDissipation()*2.0/7.0);
        speciesMixedBigScrew -> setTorsionStiffness(speciesMixedBigScrew -> getStiffness()*2.0/7.0);
        speciesMixedBigScrew -> setTorsionFrictionCoefficient(particleWallTorsionFriction);
        
    }
    
    // defines the boundaries for the simulation
    void setupBoundaries()
    {
        
        setXMax(1.5*screwCasingRadius);
        setXMin(-1.5*screwCasingRadius);
        setYMax(4.0*screwCasingRadius);
        setYMin(-1.5*screwCasingRadius);
        setZMax(screwLength);
        setZMin(0.0);
        
    }
    
    // creates the walls according to the boundaries
    void setupWalls()
    {
        
        wallHandler.clear();
        wall.setSpecies(speciesScrew);
        
        wall.set(Vec3D(-1.,0.,0.),Vec3D(getXMin(),0.,0.));
        wallHandler.copyAndAddObject(wall);
        wall.set(Vec3D(1.,0.,0.),Vec3D(getXMax(),0.,0.));
        wallHandler.copyAndAddObject(wall);
        
        wall.set(Vec3D(0,-1.,0.),Vec3D(0.,getYMin(),0.));
        wallHandler.copyAndAddObject(wall);
        wall.set(Vec3D(0.,1.,0.),Vec3D(0.,getYMax(),0.));
        wallHandler.copyAndAddObject(wall);
        
        wall.set(Vec3D(0.,0.,-1.),Vec3D(0.,0.,getZMin()));
        wallHandler.copyAndAddObject(wall);
        wall.set(Vec3D(0.,0.,1.),Vec3D(0.,0.,getZMax()));
        wallHandler.copyAndAddObject(wall);
        
        boundary.set(Vec3D(0.,0.,1.), getZMin(), getZMax());
        boundaryHandler.copyAndAddObject(boundary);
        
    }
    
    // creates the screw
    void setupScrew()
    {
        
        shaft = new AxisymmetricIntersectionOfWalls();
        shaft -> setSpecies(speciesScrew);
        shaft -> setPosition(screwOrigin);
        shaft -> setOrientation(Vec3D(0.,0.,1.));
        shaft -> addObject(Vec3D(-1.,0.,0.), Vec3D(screwShaftRadius,0.,0.));
        shaft -> setAngularVelocity(Vec3D(0.,0.,0.));
        wallHandler.addObject(shaft);
        
        //helicoid = new Helicoid03();
        auto helicoid = wallHandler.copyAndAddObject(Helicoid03());
        helicoid -> setSpecies(speciesScrew);
        helicoid -> set(screwOrigin, screwLength, screwBladeRadius, screwNumberOfTurns, screwAngularVelocity, screwThickness);
        
    }
    
    // computes the needed number of particles
    // todo : explain how this works
    void computeParticleNumber()
    {
        
        double boxVolume = screwLength*(9.*pow(screwCasingRadius,2) - constants::pi*pow(screwShaftRadius,2));
        double boxToSmallParticleVolumeRatio = 0.60*3.*boxVolume/(4.*constants::pi*pow(particleRadiusSmall,3.));
        
        if (biDisperseSystem) boxToSmallParticleVolumeRatio *= (9./16.);
        
        nTotParticles = (int)(boxToSmallParticleVolumeRatio) + 9. - fmod((int)(boxToSmallParticleVolumeRatio),9.);
        
        nBigParticles = (biDisperseSystem ? nTotParticles/9. : 0.);
        nSmallParticles = nTotParticles - nBigParticles;
        
    }
    
    // computes the mass of a sphere given its radius and its density
    void computeMassSphere()
    {
        
        particleMassSmall = 4.*constants::pi*pow(particleRadiusSmall,3.)/3.*particleDensity;
        particleMassBig = 4.*constants::pi*pow(particleRadiusBig,3.)/3.*particleDensity;
        
    }
    
    // computes the number of particles to be loaded in each loading cycle
    // todo : explain how this works (evaluates volume, assumes PF to be .25, rounds to 9)
    void computeParticlesPerLoadingCycle()
    {
        
        double loadingRegionVolume = 3.*screwLength*pow(screwCasingRadius,2.);
        double loadingRegionToSmallParticleVolumeRatio = 0.25*3.*loadingRegionVolume/(4.*constants::pi*pow(particleRadiusSmall,3.));
        
        if (biDisperseSystem) loadingRegionToSmallParticleVolumeRatio *= (9./16.);
        
        nParticlesPerLoadingCycle = ((int)(loadingRegionToSmallParticleVolumeRatio) + 9. - fmod((int)(loadingRegionToSmallParticleVolumeRatio),9.));
        nLoadingCycles = (int)(nTotParticles/nParticlesPerLoadingCycle) + 1;
        
    }
    
    // loads the particles in the system
    void particleLoader()
    {
        
        p.setVelocity(Vec3D(0.,0.,0.));
        
        for (int i=0; i <nParticlesPerLoadingCycle; i++)
        {
            
            if (!fmod(i,9) && biDisperseSystem)
            {
                
                p.setSpecies(speciesBig);
                p.setRadius(particleRadiusBig*(1.+random.getRandomNumber(-particleDispersityBig, particleDispersityBig)));
                
            }   else
            {
                
                p.setSpecies(speciesSmall);
                p.setRadius(particleRadiusSmall*(1.+random.getRandomNumber(-particleDispersitySmall, particleDispersitySmall)));
                
            }
            
            p.setPosition(Vec3D(
                                random.getRandomNumber(getXMin()+p.getRadius(), getXMax()-p.getRadius()),
                                random.getRandomNumber(2.0*screwCasingRadius + p.getRadius(), getYMax()-p.getRadius()),
                                random.getRandomNumber(getZMin()+p.getRadius(), getZMax()-p.getRadius())
                                ));
            particleHandler.copyAndAddObject(p);
            
        }

    }
    
    // estimates the physical time needed for the particle settling
    // todo : explain how this works (evaluates time from s = s0 + v0*t - .5*g*t^2 assuming s0=v0=0 and s=1.5boxSize)
    void settlingTimeEstimator()
    {
        
        loadingCycleTime = 2.0*sqrt(2.*1.5*screwCasingRadius/9.81);
        loadingTotalTime = loadingCycleTime * nLoadingCycles + 2.0;
        
    }
    
    // returns the estimated physical time duration of the settling process
    double getEstimatedTime()
    {
        
        return loadingTotalTime;
        
    }
    
    // checks if the particles from 0 to nP lay outside the screw casing and, if so, removes them
    void particleTrimmer(int nP)
    {
        
        for (int i=nP-1; i>=0; i--)
        {
            
            if (pow(particleHandler.getObject(i) -> getPosition().X - screwOrigin.X,2.) +
                pow(particleHandler.getObject(i) -> getPosition().Y - screwOrigin.Y,2.) >
                pow(screwCasingRadius - .75*(particleHandler.getObject(i) -> getRadius()),2.)) particleHandler.removeObject(i);
                
        }
        
    }
    
    // creates the screw casing
    void setupCasing()
    {
        
        casing = new AxisymmetricIntersectionOfWalls();
        casing -> setSpecies(speciesScrew);
        casing -> setPosition(screwOrigin);
        casing -> setOrientation(Vec3D(0.,0.,1.));
        casing -> addObject(Vec3D(1.,0.,0.), Vec3D(screwCasingRadius,0.,0.));
        casing -> setAngularVelocity(Vec3D(0.,0.,0.));
        wallHandler.addObject(casing);
        
    }
    
    // overrides the solve check to exit when at stage 3 and all particles are still
    bool continueSolve() const override
    {
        
        // after the cut, whenever the particles are still writes the restart file and quits
        if (stage == 2 && getKineticEnergy()/getElasticEnergy() < 1.e-4)
        {
            std::cout << "\nWriting the restart file and quitting...\n";
            return false;
            
        }
        
        return true;
        
    }
    
    // the counter needed to monitor the loading iterations
    int loadCyclesCounter = 0;
    
    // the various stages of the loading process
    int stage = 0;
    
    
private:
    // particle properties
    double particleRadiusSmall;
    double particleRadiusBig;
    double particleDispersitySmall;
    double particleDispersityBig;
    double particleMassSmall;
    double particleMassBig;
    
    double particleDensity;
    
    double particleParticleSlidingFriction;
    double particleParticleRollingFriction;
    double particleParticleTorsionFriction;
    double particleWallSlidingFriction;
    double particleWallRollingFriction;
    double particleWallTorsionFriction;
    
    double particleParticleRestitutionCoefficient;
    double particleWallRestitutionCoefficient;
    
    int nTotParticles;
    int nSmallParticles;
    int nBigParticles;
    int nParticlesPerLoadingCycle;
    int nLoadingCycles;
    
    bool biDisperseSystem;
    
    // screw properties
    double screwLength;
    double screwCasingRadius;
    double screwBladeRadius;
    double screwShaftRadius;
    double screwNumberOfTurns;
    double screwThickness;
    double screwAngularVelocity;
    
    Vec3D screwOrigin;

    // system variables
    double collisionTime;
    double loadingCycleTime;
    double loadingTotalTime;
    
    int numberOfParticleLoadingCycles;
    
    // physical objects
    SphericalParticle p;
    
    InfiniteWall wall;
    AxisymmetricIntersectionOfWalls *shaft;
    AxisymmetricIntersectionOfWalls *casing;
    //Helicoid03 *helicoid;
    PeriodicBoundary boundary;
    
    LinearViscoelasticFrictionSpecies *speciesSmall;
    LinearViscoelasticFrictionSpecies *speciesBig;
    LinearViscoelasticFrictionSpecies *speciesScrew;
    
};



int main(int argc UNUSED, char *argv[] UNUSED)
{
    
    ScrewFiller filler;
//    problem.autoNumber();
    filler.setName("ScrewFiller");
    
    filler.setGravity(-9.81*Vec3D(0.,1.,0.));
    filler.setSystemDimensions(3);
    
    filler.setCollisionTime(0.005);
    filler.setTimeStep(0.00001);
    filler.setSaveCount(0.01/filler.getTimeStep());
//    filler.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
//    filler.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    
    // set file type
    
    // radius of small particles, their dispersity and, if present, dispersity of the big ones
    // overloaded function : if the third argument is present then the bool biDispersity is turned to true
    filler.setParticleRadiusAndDPolyspersity(0.01, 0.1);
    
    filler.setParticleDensity(2000.);
    filler.setSlidingFrictionCoefficients(0.2, 0.5);   // p-p and p-w sliding friction coefficients
    filler.setRollingFrictionCoefficients(0.1, 0.2);   // p-p and p-w rolling friction coefficients
    filler.setTorsionFrictionCoefficients(0.01, 0.1);   // p-p and p-w torsion friction coefficients
    filler.setRestitutionCoefficients(0.8, 0.9);    // p-p and p-w restitution coefficients
    
    // length, casing radius, blade radius, shaft radius, number of turns, thickness
    filler.setScrewGeometry(.3, .195, .175, .075, 1.0, 0.08);
    filler.setScrewOrigin(Vec3D(0.,0.,0.));
    filler.setScrewVelocity(0.);
    
    filler.solve();
    
    return 0;
    
}




/*  --- JUNKYARD ---
 
 --- FUNCTIONS CALLED IN THE MAIN ---
 void setScrewOperationalParameters(...)    // todo
 
 --- MAIN ---

 
*/







