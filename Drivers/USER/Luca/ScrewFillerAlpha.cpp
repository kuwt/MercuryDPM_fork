
#include "Mercury3D.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Walls/InfiniteWallWithHole.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Helicoid03.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "sstream"

/*
 ToDo:
 - should check for ovelap during the insertion (not really useful if blade AND particles are inflated)
 - inflate the particle radius
 - inflate the thickness also
 - correct the magic number in the radius change
 
 - correct the particle-screw_side interaction
 - adjust material properties
 - clean everything and make it consistent
 - fix the collision rule
 
 - reorder everything and eliminate the comments
 - put the right material properties and species
 - put the variables in the protected space and initialize them form the main
 */

class ScrewFiller : public Mercury3D {
    
private:
    const double rMax = .15;
    const double rMin = .05;
    const double length = 0.5;
    const double thickness = 0.06;
    const double screwCasingRadius = 1.0*rMax;
    
    const int nP = 1;
    // fix the particle diameter as half of the blade radius
    double particleRadius = 0.25*(rMax - rMin);

    int cuttingCycle = 0;
    int stage = 1;
    
    bool biDispersity = false;
    
    void setupInitialConditions() override
    {
        
        std::cout << "\nThe shaft rescaled radius is " << .5*(1. - rMin/screwCasingRadius) << ".\n";
        // gravity, particle radius
//        setGravity(-9.81*Vec3D(0.0, sqrt(3.0), 1.0)/2.0);
        setGravity(-9.81*Vec3D(0.0, 1.0, 0.0));
        
        // set problem geometry
        setXMax(screwCasingRadius);
        setXMin(-screwCasingRadius);
        
        setZMax(length);
        setZMin(0.0);
        
//        setYMax(1.2*5.0*rMax);
        setYMax(screwCasingRadius);
        setYMin(-screwCasingRadius);
        
        // set problem species
        speciesHandler.clear();
        
        particleMass = 4.*constants::pi*pow(particleRadius,3.)/3.*particleDensity;
        
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
        
        // periodic boundary in the screw's axis direction
        b0 = boundaryHandler.copyAndAddObject(PeriodicBoundary());
        b0->set(Vec3D(0,0,1), getZMin(), getZMax());
        
        xMinWall = wallHandler.copyAndAddObject(InfiniteWall());
        xMinWall->set(Vec3D(-1,0,0),Vec3D(getXMin(),0,0));
        
        xMaxWall = wallHandler.copyAndAddObject(InfiniteWall());
        xMaxWall->set(Vec3D(1,0,0),Vec3D(getXMax(),0,0));
        
        yMinWall = wallHandler.copyAndAddObject(InfiniteWall());
        yMinWall->set(Vec3D(0,-1,0),Vec3D(0,getYMin(),0));
        
        yMaxWall = wallHandler.copyAndAddObject(InfiniteWall());
        yMaxWall->set(Vec3D(0,1,0),Vec3D(0,getYMax(),0));

        // inner cylinder (the screw shaft) [I guess I have to use the axisymmetric wall]
        screwShaft = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        screwShaft->setPosition(Vec3D(0,0,0));
        screwShaft->setOrientation(Vec3D(0,0,1));
        screwShaft->addObject(Vec3D(-1,0,0),Vec3D(rMin,0,0));
        screwShaft->setSpecies(specieParticle);
        
        // creation of the screw and setting of its properties
        screw = wallHandler.copyAndAddObject(Helicoid03());
        screw -> set(Vec3D(0., 0., 0.), length, rMax, 1.0, constants::pi, thickness);
        screw -> setSpecies(specieWall);
        
        // external case
        screwCase = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        screwCase->setPosition(Vec3D(0,0,0));
        screwCase->setOrientation(Vec3D(0,0,1));
        screwCase->addObject(Vec3D(1,0,0),Vec3D(screwCasingRadius,0,0));
        screwCase->setSpecies(specieParticle);
        
        // particle creation
        // particle radius is uniformly randomly picked depending on the dispersity
        // dispersity of 1 means radius +- 1*radius
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0.setSpecies(specieParticle);
        p0.setRadius(particleRadius);
        p0.setPosition(Vec3D(0.0, -(rMax - particleRadius), particleRadius));
        particleHandler.copyAndAddObject(p0);
        
//        for (int i = 0; i < nP; i++)
//        {
//            
//            if (biDispersity)
//            {
//                
//                if (fmod(i,9.0))
//                {
//                    p0.setRadius(particleRadius*random.getRandomNumber(1.0-dispersity,1.0+dispersity));
//                }
//                else
//                {
//                    p0.setRadius(2.0*particleRadius*random.getRandomNumber(1.0-dispersity,1.0+dispersity));
//                }
//                
//            }
//            else
//            {
//                p0.setRadius(particleRadius*random.getRandomNumber(1.0-dispersity,1.0+dispersity));
//            }
//            
//            p0.setPosition(Vec3D(random.getRandomNumber(getXMin(),getXMax()),random.getRandomNumber(1.2*rMax,getYMax()),random.getRandomNumber(.0, length)));
//            particleHandler.copyAndAddObject(p0);
//            
//        }
        
    }
    
    //! [CST:beforetime]
    void actionsBeforeTimeStep() override
    {
        screwShaft->setOrientation(Vec3D(0.0,0.0,1.0));
        
        // settles for 2 seconds then starts with the rotation
        if (getTime() > 2.0)
        {
            screw->move_time(getTimeStep());
            screwShaft->setAngularVelocity(Vec3D(0.0,0.0,constants::pi));
            
            if (fmod(getTime() - 2.0, frictionObservationTime) < getTimeStep() && getTime() > frictionObservationTime)
            {
                particleWallSlidingFriction += 0.1;
                
                std::cout << "\nFriction increment. New value of particle-blade sliding friction: " << particleWallSlidingFriction << "\n";
            }
        }
        
        
//        // rotates the helicoid
//        if (getTime() > 5.0)
//        {
////            screw.incrementOffset(getTimeStep()*constants::pi);
//            screw->move_time(getTimeStep());
//            screwShaft->setAngularVelocity(Vec3D(0.0,0.0,constants::pi));
//        }
        
    }

public:
//    double particleRadius;
    double particleMass;
    double particleDensity;
    
    double frictionObservationTime;
    
    double collisionTime;
    
    Mdouble dispersity = 0.0;
    Mdouble fillingRatio = .50;
    
    LinearViscoelasticFrictionSpecies *specieParticle, *specieWall, *specieMixedParticleWall;
    SphericalParticle p0;
    Helicoid03 *screw;
    InfiniteWall *xMinWall, *xMaxWall;
    InfiniteWall *yMinWall, *yMaxWall;
    AxisymmetricIntersectionOfWalls* screwCase;
    AxisymmetricIntersectionOfWalls* screwShaft;
    PeriodicBoundary* b0;
    
    double particleParticleSlidingFriction;
    double particleParticleRollingFriction;
    double particleParticleTorsionFriction;
    double particleWallSlidingFriction;
    double particleWallRollingFriction;
    double particleWallTorsionFriction;
    double particleParticleRestitutionCoefficient;
    double particleWallRestitutionCoefficient;
};


int main(int argc UNUSED, char *argv[] UNUSED)
{
    
    ScrewFiller problem;
    
    problem.setSystemDimensions(3);
    
    problem.setName("FrictionTest_slid_1.0_var_roll_0.0_0.0_tors_0.0_0.0");
    
    problem.frictionObservationTime = 10;
    
    problem.collisionTime = 0.001;
    problem.setTimeStep(0.001/50.0);
    problem.setTimeMax(11*problem.frictionObservationTime + 2.0);   // 10 friction increments: from 0.0 t0 2.0 included
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(10000, problem.getTimeMax(), problem.getTimeStep()));
    
//    problem.particleRadius = 0.015;
    problem.particleDensity = 2000.0;
    
    problem.particleParticleRestitutionCoefficient = 0.8;
    problem.particleWallRestitutionCoefficient = 0.8;
    
    problem.particleParticleSlidingFriction = 1.0;
    problem.particleWallSlidingFriction = 0.0;
    
    problem.particleParticleRollingFriction = 0.0;
    problem.particleWallRollingFriction = 0.0;
    
    problem.particleParticleTorsionFriction = 0.0;
    problem.particleWallTorsionFriction = 0.0;
    
    // actually solving the problem
    problem.solve();
    
//    // run the following for parameter study
//    double mu;
//    
//    for (int i = 1; i <= 11; i++)
//    {
//        ScrewFiller problem;
//        
//        problem.setSystemDimensions(3);
//        problem.setTimeStep(0.008/50.0);
//        problem.setTimeMax(10.0);
//        
//        mu = 0.02*(i-1);
//        problem.mu = mu;
//        std::ostringstream oss;
//        oss << "ScrewFillerAlpha_" << mu;
//        
//        // set some basic problem properties
//        problem.setName(oss.str());
//        
//        problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(1000, problem.getTimeMax(), problem.getTimeStep()));
//        
//        
//        // actually solving the problem
//        problem.solve();
//    }
    
}


// the stuff below was fo cutting the particles outside of the case and keeping only one
//        // 1st cutting cycle: cuts outside the casing and creates the casing
//        if (getTime() > 3.0 && cuttingCycle == 0)
//        {
//            for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
//            {
//                if (pow(particleHandler.getObject(i) -> getPosition().X,2.) +
//                    pow(particleHandler.getObject(i) -> getPosition().Y,2.) >
//                    pow(screwCasingRadius - .95*(particleHandler.getObject(i) -> getRadius()),2.)) particleHandler.removeObject(i);
//            }
//            screwCase = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
//            screwCase->setPosition(Vec3D(0,0,0));
//            screwCase->setOrientation(Vec3D(0,0,1));
//            screwCase->addObject(Vec3D(1,0,0),Vec3D(screwCasingRadius,0,0));
//            screwCase->setSpecies(specieParticle);
//
//            cuttingCycle++;
//        }
//
//        // 2nd cutting cycle: cuts to achieve the desided filling ratio and/or number of particles
//        if (getTime() > 4.0 && cuttingCycle == 1)
//        {
//            for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
//            {
//                if (particleHandler.getObject(i) -> getPosition().Y - 0. > (2.0*fillingRatio-1.0)*screwCasingRadius) particleHandler.removeObject(i);
//            }
//            for (int i=particleHandler.getNumberOfObjects()-1; i>0; i--) {particleHandler.removeObject(i);}
//
//            cuttingCycle++;
//        }
//
//        if (cuttingCycle == 2)
//        {
//            setGravity(Vec3D(0.0,0.0,-9.81));
//            cuttingCycle++;
//        }

