
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
 - bug in the particle cutting: wrong number of objects (in the mono-modal case)
    probably the mistake is here 
        nLoadingCycles = (int)(nTotParticles/nParticlesPerLoadingCycle);
    since, according to the rounding, I can take different nLoadingCycles
    which might lead to nLoadingCycles*nParticlesPerLoadingCycle != nTotParticles
    *** I THINK NOW IT'S WORKING ***
 
 MINOR
 - write teh checks that exit the execution if parameters make no sense
 - ask if special things have to be accounted for when writing the restart file
 
 WHATEVER
 - create a .h and .cc of this stuff
 - write the comments in doxigen style and standardise them
 - provide a talkative version
 - write a small thing that outputs the spatial distribution of the big/small particles to check if isotropic (cut in cake style)
 
 DEVELOPMENT

 */

class ScrewTrimmer : public Mercury3D
{
    
public:
//  --- STANDARD FUNCTIONS ---
    
    void setupInitialConditions() override
    {
        
        wallHandler.clear();
        setupBoundaries();
        std::cout << "\nBoundaries DONE!\nSetting up walls...";
        setupWalls();
        std::cout << "\nWalls DONE!\nSetting up screw...";
        setupScrew();
        std::cout << "\nScrew DONE!\nCutting the particles...";
        particleTrimmer();
        std::cout << "\nCutting DONE!\nPlacing the external casing...";
        //ToDo
        setupCasing();
        std::cout << "\nCasing DONE!\nWaiting 2 seconds for the particles to settle...";
        
    }
    
    void actionsBeforeTimeStep() override
    {

    }
    
    
//  --- FUNCTIONS CALLED IN THE MAIN ---
    
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
    
    // checks if the particles from 0 to nP lay outside the screw casing and, if so, removes them
    void particleTrimmer()
    {
        
        for (int i=particleHandler.getNumberOfObjects()-1; i >=0; i--)
        {
            
            if (pow(particleHandler.getObject(i) -> getPosition().X - screwOrigin.X,2.) +
                pow(particleHandler.getObject(i) -> getPosition().Y - screwOrigin.Y,2.) >
                pow(screwCasingRadius - .95*(particleHandler.getObject(i) -> getRadius()),2.)) particleHandler.removeObject(i);
            
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
    
    // overrides the solve check to exit when time > 2.0s and all particles are still
    bool continueSolve() const override
    {
        
        // after the cut, whenever the particles are still writes the restart file and quits
        if (getTime() >= 2. && getKineticEnergy()/getElasticEnergy() < 1.e-4)
        {
            std::cout << "\nWriting the restart file and quitting...\n";
            return false;
            
        }
        
        return true;
        
    }
    
private:
    // particle properties
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



int main(int argc, char *argv[])
{
    
    ScrewTrimmer trimmer;
    trimmer.setName("ScrewTrimmer");
    
    std::string fileName;
    
    if (argc>1)
    {
        
        fileName = argv[1];
        std::cout << "\nRestarting the system from " << fileName << "\n";
        trimmer.readRestartFile(fileName);
        
    }
    else
    {
        
        std::cout << "\nNo restart file specified. Exiting...\n";
        exit(0);
        
    }
    
    trimmer.setGravity(-9.81*Vec3D(0.,1.,0.));
    trimmer.setSystemDimensions(3);
    
    trimmer.setCollisionTime(0.005);
    trimmer.setTimeStep(0.00001);
    trimmer.setSaveCount(0.01/trimmer.getTimeStep());
//    trimmer.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
//    trimmer.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    
    trimmer.setParticleDensity(2000.);
    trimmer.setSlidingFrictionCoefficients(0.2, 0.5);   // p-p and p-w sliding friction coefficients
    trimmer.setRollingFrictionCoefficients(0.1, 0.2);   // p-p and p-w rolling friction coefficients
    trimmer.setTorsionFrictionCoefficients(0.01, 0.1);   // p-p and p-w torsion friction coefficients
    trimmer.setRestitutionCoefficients(0.8, 0.9);    // p-p and p-w restitution coefficients
    
    // loads the screw
    // *** IT MUST MATCH THE ONE IN THE FILLER FILE! ***
    // length, casing radius, blade radius, shaft radius, number of turns, thickness
    trimmer.setScrewGeometry(.3, .195, .175, .075, 1.0, 0.08);
    trimmer.setScrewOrigin(Vec3D(0.,0.,0.));
    
    trimmer.solve();
    
    return 0;
    
}




/*  --- JUNKYARD ---
 
 --- FUNCTIONS CALLED IN THE MAIN ---
 
 --- MAIN ---

 
*/







