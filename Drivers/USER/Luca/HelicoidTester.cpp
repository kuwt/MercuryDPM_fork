
#include "Mercury3D.h"
#include "Particles/BaseParticle.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Walls/InfiniteWallWithHole.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Helicoid03.h"
#include "Helicoid04.h"
#include "Helicoid03bis.h"
#include "Helicoid05.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "sstream"
#include "Shaft.h"

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

class HelicoidTester : public Mercury3D {
    
private:
    const double rMax = .15;
    const double rMin = .05;
    const double length = 0.5;
    const double thickness = 0.06;
    const double screwCasingRadius = 1.2*rMax;
    
    // should be N*particleRadius = rMax - rMin = 0.1  with  N > 6 (i.e. 3 particles are fitting along the blade width)
    // now N = 6 -> pR = 0.0167
    const double particleRadius = 0.010;
    
    // nP = 3500 fits for pR = 0.0135; rescale with a factor new_nP = old_nP*pow(old_pR/new_pR,3)
//    const int nP = (int)(3500*pow(0.0135/particleRadius,3.0));
    const int nP = 12000; // 2500

    int cuttingCycle = 0;
    
    bool biDispersity = false;
    
    void setupInitialConditions()
    {
        
        std::cout << "\nThe shaft rescaled radius is " << .5*(1. - rMin/screwCasingRadius) << ".\n";
        // gravity, particle radius
        setGravity(Vec3D(0.0,-9.81,0.0));
        
        // set problem geometry
        setXMax(1.2*screwCasingRadius);
        setXMin(-1.2*screwCasingRadius);
        
        setZMax(length);
        setZMin(0.0);
        
        setYMax(1.2*5.0*rMax);
        setYMin(-1.2*screwCasingRadius);
        
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
        screwShaft->setSpecies(specieParticle);
        screwShaft->setPosition(Vec3D(0,0,0));
        screwShaft->setOrientation(Vec3D(0,0,1));
        screwShaft->addObject(Vec3D(-1,0,0),Vec3D(rMin,0,0));
        
        // creation of the screw and setting of its properties
        screw = wallHandler.copyAndAddObject(Helicoid04());
//        screw -> set(Vec3D(0., 0., 0.), length, rMax, rMin, 1.0, constants::pi, thickness, true);
        screw -> set(Vec3D(0., 0., 0.*length), length, rMax, rMin, 1.0, constants::pi, thickness);
        screw -> setSpecies(specieWall);
        
        // particle creation
        // particle radius is uniformly randomly picked depending on the dispersity
        // dispersity of 1 means radius +- 1*radius
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0.setSpecies(specieParticle);
        
        for (int i = 0; i < nP; i++)
        {
            
            if (biDispersity)
            {
                
                if (fmod(i,9.0))
                {
                    p0.setRadius(particleRadius*random.getRandomNumber(1.0-dispersity,1.0+dispersity));
                }
                else
                {
                    p0.setRadius(2.0*particleRadius*random.getRandomNumber(1.0-dispersity,1.0+dispersity));
                }
                
            }
            else
            {
                p0.setRadius(particleRadius*random.getRandomNumber(1.0-dispersity,1.0+dispersity));
            }
            
            p0.setPosition(Vec3D(random.getRandomNumber(getXMin(),getXMax()),random.getRandomNumber(1.2*rMax,getYMax()),random.getRandomNumber(.0, length)));
            particleHandler.copyAndAddObject(p0);
            
        }
        
//        p0.setRadius(particleRadius*random.getRandomNumber(1.0-dispersity,1.0+dispersity));
//        p0.setPosition(Vec3D(rMin + 1.2*particleRadius,0.0,0.75*length + 1.2*particleRadius));
//        particleHandler.copyAndAddObject(p0);
        
        packingFractionScalingCoefficient = 0.25*3.0*(length*(pow(screwCasingRadius, 2) - pow(rMin, 2)) - thickness*(1.0 + pow(length/(constants::pi*(rMax + rMin)), 2))*(pow(rMax, 2) - pow(rMin, 2)));
        
        std::cout << "\nThe packing fraction scaling coefficient is " << packingFractionScalingCoefficient << "\n";
    }
    
    //! [CST:beforetime]
    void actionsBeforeTimeStep()
    {
        
        //        wallHandler.getObject(4)->setOrientation(Vec3D(0,0,1));
        screw->setOrientation(Vec3D(0,0,1));
        //        wallHandler.getObject(5)->setOrientation(Vec3D(0,0,1));
        
                screwShaft->setOrientation(Vec3D(0.0,0.0,1.0));
        
        // 1st cutting cycle: cuts outside the casing and creates the casing
        if (getTime() > 3.0 && cuttingCycle == 0)
        {
            //            double volParticles;
            //
            //            // computes the volume inside concentric cylinders until it meets the desired final packing fraction
            //            do
            //            {
            //                particleEliminationLoopCounter++;
            //
            //                // initializes the volume of the particles
            //                volParticles = 0.0;
            //
            //                // loops over all the particles and computes the volume of the ones inside concentric shells
            //                for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
            //                {
            //                    if (pow(particleHandler.getObject(i) -> getPosition().X, 2.) + pow(particleHandler.getObject(i) -> getPosition().Y, 2.) < pow(0.75*rMax + 0.1*particleEliminationLoopCounter*(1.0 - dispersity)*particleRadius, 2.)) volParticles += pow(particleHandler.getObject(i) -> getRadius(), 3);
            //                }
            //
            //                std::cout << "\n" << volParticles;
            //            } while (volParticles < packingFractionScalingCoefficient*packingFraction);
            //
            //            // cuts to achieve a desired packing fraction with the totally filled configuration
            //            for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
            //            {
            //                if (pow(particleHandler.getObject(i) -> getPosition().X, 2.) +
            //                    pow(particleHandler.getObject(i) -> getPosition().Y, 2.) >
            //                    pow(0.75*rMax + 0.1*particleEliminationLoopCounter*(1.0 - dispersity)*particleRadius ,2.)) particleHandler.removeObject(i);
            //            }
            //
            //            std::cout << "\nRadius_cut/casing_radius" << (0.75*rMax + 0.1*particleEliminationLoopCounter*(1.0 - dispersity)*particleRadius)/screwCasingRadius << "\n";
            
            // cuts to achieve the desired filling ratio
            for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
            {
                if (pow(particleHandler.getObject(i) -> getPosition().X,2.) +
                    pow(particleHandler.getObject(i) -> getPosition().Y,2.) >
                    pow(1.05*screwCasingRadius - .95*(particleHandler.getObject(i) -> getRadius()),2.)) particleHandler.removeObject(i);
            }
            
            screwCase = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
            screwCase->setPosition(Vec3D(0,0,0));
            screwCase->setOrientation(Vec3D(0,0,1));
            screwCase->addObject(Vec3D(1,0,0),Vec3D(screwCasingRadius,0,0));
            screwCase->setSpecies(specieParticle);
            
            cuttingCycle++;
        }
        
        // 2nd cutting cycle: cuts to achieve the desided filling ratio and/or number of particles
        if (getTime() > 4.0 && cuttingCycle == 1)
        {
            for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
            {
                if (particleHandler.getObject(i) -> getPosition().Y - 0. > (2.0*fillingRatio-1.0)*screwCasingRadius) particleHandler.removeObject(i);
            }
            
            cuttingCycle++;
        }
        
        if (cuttingCycle == 2)
        {
//            screw->rotate(getTimeStep());
            screw->move_time(getTimeStep());
            screw->setOrientation(Vec3D(0,0,1));
            screw->setAngularVelocity(Vec3D(0.0,0.0,-constants::pi));
            
            screwShaft->setOrientation(Vec3D(0,0,1));
            screwShaft->setAngularVelocity(Vec3D(0.0,0.0,-constants::pi));
            
//            if (screwRightHandeness)
//            {
//                screw->setAngularVelocity(Vec3D(0.0,0.0,-constants::pi));
//            }
//            else
//            {
//                screw->setAngularVelocity(Vec3D(0.0,0.0,constants::pi));
//            }
        }
    }

public:
    double particleMass;
    double particleDensity;
    
    double collisionTime;
    
    double packingFraction;
    double packingFractionScalingCoefficient;
    int particleEliminationLoopCounter = 0;
    
    Mdouble dispersity = 0.0;
    Mdouble fillingRatio;
    
    LinearViscoelasticFrictionSpecies *specieParticle, *specieWall, *specieMixedParticleWall;
    BaseParticle p0;
    Helicoid04 *screw;
//    Helicoid03bis *screw;
    InfiniteWall *xMinWall, *xMaxWall;
    InfiniteWall *yMinWall, *yMaxWall;
    AxisymmetricIntersectionOfWalls* screwCase;
    AxisymmetricIntersectionOfWalls* screwShaft;
    PeriodicBoundary* b0;
    
    bool screwRightHandeness;
    
//    Shaft* ssshaft;
    
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
    
    HelicoidTester problem;
    
    problem.setSystemDimensions(3);
    
//    problem.setName("HelicoidTester_Helicoid03bis");
    problem.setName("HelicoidTester_Helicoid04_forPARTICLES2017");
    
    problem.collisionTime = 0.002;
    problem.setTimeStep(0.002/50.0*2.5);
    problem.setTimeMax(10.0);
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimestep(1000, problem.getTimeMax(), problem.getTimeStep()));
    
    problem.particleDensity = 2000.0;
    
    problem.screwRightHandeness = true;
    
    problem.packingFraction = 0.7;
    problem.fillingRatio = 1.0;
    
    problem.particleParticleRestitutionCoefficient = 0.8;
    problem.particleWallRestitutionCoefficient = 0.8;
    
    problem.particleParticleSlidingFriction = 0.3;
    problem.particleWallSlidingFriction = 0.3;
    
    problem.particleParticleRollingFriction = 0.0;
    problem.particleWallRollingFriction = 0.0;
    
    problem.particleParticleTorsionFriction = 0.0;
    problem.particleWallTorsionFriction = 0.0;
    
    problem.setXBallsAdditionalArguments("-h 800 -p 10 -o 200 -3dturn 1");
    
    // actually solving the problem
    problem.solve();
    
}

