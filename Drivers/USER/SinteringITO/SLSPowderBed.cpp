//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Mercury3D.h"
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <Species/HertzianViscoelasticFrictionSpecies.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <map>
#include <random>
using constants::pi;
using mathsFunc::cubic;
enum class RollerType : unsigned char {
    NO_ROLLER = 0,
    FORWARD_ROTATING_ROLLER = 1,
    COUNTER_ROTATING_ROLLER = 2,
    FORWARD_ROTATING_ROLLER_WITH_BLADE = 3
};

/**
 * Simulates the compaction of an SLS powder bed by different methods, e.g. roller compaction
 *
 * The basic setup (defined in the constructor) consists of
 * - an empty rectangular domain periodic in x and y direction, and constrained by a base plate in z-direction
 * - a species containing the material properties (elastic modulus, shear modulus, friction, restitution, material density)
 * - a particle size distribution, either Gaussian or uniform in size, with given mean and std deviation.
 * The constructor then calls solve twice, first to create a basal layer of fixed particles,
 * then to create a layer of particles. Both layers are created by raining down particles.
 *
 * The user can then add certain features:
 * - a forward-rotating roller (parameters: roller radius)
 * - a counter-rotating roller (parameters: roller radius)
 * - a forward-rotating roller with a blade (parameters: roller radius, roller-blade distance, blade thickness)
 * For each compactor, a position function needs to be specified.
 */
#include "Mercury3D.h"
#include "Species/HertzianViscoelasticFrictionSpecies.h"

class PowderBed : public Mercury3D
{
public:

    PowderBed (Mdouble domainLength, Mdouble domainWidth, Mdouble domainDepth,
               ParticleSpecies& particleSpecies, BaseSpecies& wallSpecies, Mdouble relativeBasalLayerThickness=1.7)
        : relativeBasalLayerThickness_(relativeBasalLayerThickness)
    {
        //set name, gravity, output options
        setName("SLSPowderBed");
        setGravity({0,0,-9.8});
        setXBallsAdditionalArguments("-solidf -v0");
        logger(INFO,"Name of output file: %",getName());

        logger(INFO,"Defining domain as [0,%]x[-%,%]x[0,%]",domainLength,0.5*domainWidth,domainDepth,domainDepth);
        setMin({0,-0.5*domainWidth,0});
        setMax({domainLength,0.5*domainWidth,domainDepth});

        logger(INFO,"Adding contact laws");
        speciesHandler.copyAndAddObject(particleSpecies); //particle-paricle interactions
        speciesHandler.copyAndAddObject(particleSpecies); //wall-wall interactions (unused)
        speciesHandler.getMixedObject(0,1)->mixAll(&wallSpecies,&wallSpecies); //particle-wall interactions

        logger(INFO,"Setting up periodic boundary conditions in x and y");
        PeriodicBoundary p;
        p.set({1,0,0},getXMin(),getXMax());
        boundaryHandler.copyAndAddObject(p);
        p.set({0,1,0},getYMin(),getYMax());
        boundaryHandler.copyAndAddObject(p);

        logger(INFO,"Creating a base wall at z=0");
        InfiniteWall w;
        w.setSpecies(speciesHandler.getObject(0));
        w.set({0,0,-1},{0,0,getZMin()});
        wallHandler.copyAndAddObject(w);
    }

    void setGaussianDistribution (Mdouble meanRadius, Mdouble stdRadius) {
        //set up random number generator
        std::random_device rd;
        std::mt19937 gen(rd());

        //set up normal distribution
        std::normal_distribution<> d(meanRadius, stdRadius);

        //move down basal plate by relativeBasalLayerThickness_*particle radius
        setZMin(-relativeBasalLayerThickness_*meanRadius);
        wallHandler.getObject(0)->setPosition({0,0,getZMin()});

        //the volume to be added
        Mdouble addVolume = 0.65*(getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin());

        //add particles until the volume to be added is zero
        logger(INFO,"Adding particles ...");
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(meanRadius);
        Mdouble fillHeight = getZMin();
        while (addVolume>0) {
            Mdouble x = random.getRandomNumber(getXMin(), getXMax());
            Mdouble y = random.getRandomNumber(getYMin(), getYMax());
            Mdouble z = random.getRandomNumber(getZMin(), fillHeight);
            p.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p)) {
                particleHandler.copyAndAddObject(p);
                addVolume -= p.getVolume();
                do {
                    p.setRadius(d(gen));
                } while (p.getRadius()<0.2*meanRadius || p.getRadius()>1.8*meanRadius); //reject too small or large radii
                if (particleHandler.getNumberOfObjects()%100==0) std::cout << '.' << std::flush;
                if (particleHandler.getNumberOfObjects()%1000==0) std::cout << ' ';
                if (particleHandler.getNumberOfObjects()%10000==0) std::cout << addVolume << '\n';
            } else {
                fillHeight += 0.01*meanRadius; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());
    }

    void setUniformDistribution (Mdouble minRadius, Mdouble maxRadius) {
        //set up random number generator
        std::random_device rd;
        std::mt19937 gen(rd());

        //set up normal distribution
        std::uniform_real_distribution<> d(minRadius, maxRadius);

        //move down basal plate by relativeBasalLayerThickness_*particle radius
        Mdouble meanRadius = 0.5*(minRadius + maxRadius);
        setZMin(-relativeBasalLayerThickness_*meanRadius);
        wallHandler.getObject(0)->setPosition({0,0,getZMin()});

        //the volume to be added
        Mdouble addVolume = 0.65*(getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin());

        //add particles until the volume to be added is zero
        logger(INFO,"Adding particles ...");
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(meanRadius);
        Mdouble fillHeight = getZMin();
        while (addVolume>0) {
            Mdouble x = random.getRandomNumber(getXMin(), getXMax());
            Mdouble y = random.getRandomNumber(getYMin(), getYMax());
            Mdouble z = random.getRandomNumber(getZMin(), fillHeight);
            p.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p)) {
                particleHandler.copyAndAddObject(p);
                addVolume -= p.getVolume();
                do {
                    p.setRadius(d(gen));
                } while (p.getRadius()<0.1*meanRadius || p.getRadius()>0.9*meanRadius); //reject too small or large radii
                if (particleHandler.getNumberOfObjects()%100==0) std::cout << '.' << std::flush;
                if (particleHandler.getNumberOfObjects()%1000==0) std::cout << ' ';
                if (particleHandler.getNumberOfObjects()%10000==0) std::cout << addVolume << '\n';
            } else {
                fillHeight += 0.01*meanRadius; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());
    }

    //remove particles according to (max.) outflow rate
    void actionsAfterTimeStep() override
    {
        //store time when rolling should start
        static Mdouble rollerTime = 0;
        //make decisions based on compression stage
        static unsigned compressionStage_ = 0;

        if (compressionStage_==0)
        {
            // Functions after this statement only get executed every 100th time step (if counter==100)
            static unsigned counter = 0;
            if (++counter != 100) return;
            else counter = 0;

            if (getKineticEnergy() < 5e-3 * getElasticEnergy())
            {
                //change to compression state
                compressionStage_ = 1;

                //set roller time
                rollerTime = getTime() + 2.0*getZMax() / rollerVelocity_;

                //set final time such that the roller completes the domain
                setTimeMax(rollerTime + (getXMax() - getXMin() - rollerRadius_) / rollerVelocity_);

                logger(INFO, "Starting compression at t=%, rollerTime=%, tMax=%", getTime(), rollerTime, getTimeMax());

                //fix basal layer of particles
                for (const auto p : particleHandler)
                {
                    if (p->getPosition().Z < 0.0)
                        p->fixParticle();
                }

                //create roller
                AxisymmetricIntersectionOfWalls roller;
                roller.setSpecies(speciesHandler.getObject(1));
                roller.setPosition({rollerRadius_, 0, rollerRadius_ + 2.0*getZMax() + compaction_});
                roller.setAxis({0, 1, 0});//setAxis
                roller.addObject({-1, 0, 0}, {rollerRadius_, 0, 0});
                roller.setVelocity({0,0,-rollerVelocity_});//first move the roller down
                wallHandler.copyAndAddObject(roller);

                //add blade
                if (rollerType_ == RollerType::FORWARD_ROTATING_ROLLER_WITH_BLADE)
                {
                    setTimeMax(getTimeMax() - rollerVelocity_ * (bladeRollerDistance_ + bladeThickness_));
                    IntersectionOfWalls blade;
                    blade.setSpecies(speciesHandler.getObject(1));
                    blade.addObject({-1, 0, 0},
                                    {rollerRadius_ + rollerRadius_ + bladeRollerDistance_ + bladeThickness_, 0, 0});
                    blade.addObject({1, 0, 0}, {rollerRadius_ + rollerRadius_ + bladeRollerDistance_, 0, 0});
                    blade.addObject({0, 0, 1}, {0, 0, 2.0*getZMax() + compaction_});
                    blade.setPosition({0,0,0});//first move the roller down
                    blade.setVelocity({0,0,-rollerVelocity_});//first move the roller down
                    wallHandler.copyAndAddObject(blade);
                    logger(INFO, "z=% %",wallHandler.getObject(1)->getPosition().Z-rollerRadius_,2.0*getZMax() + compaction_);
                }
            }
        } else {
            if (getTime() > rollerTime) {
                logger(INFO, "Start forward movement");
                // set velocity
                wallHandler.getObject(1)->setVelocity({rollerVelocity_,0,0});
                wallHandler.getLastObject()->setVelocity({rollerVelocity_,0,0});
                //set angular velocity
                if (rollerType_ == RollerType::COUNTER_ROTATING_ROLLER)
                {
                    wallHandler.getObject(1)->setAngularVelocity({0,-rollerVelocity_/rollerRadius_,0});
                } else {
                    wallHandler.getObject(1)->setAngularVelocity({0,rollerVelocity_ /rollerRadius_,0});
                }
                logger(INFO, "a=%",wallHandler.getObject(1)->getAngularVelocity());
                logger(INFO, "z=%",wallHandler.getObject(1)->getPosition().Z-rollerRadius_);
                logger(INFO, "z=% %",compaction_,getZMax());
                // make sure this is set only once
                rollerTime = constants::inf;
            }
        }
    }

    // display time and eneRatio every time the files are printed
    void printTime() const override
    {
        std::cout << "t " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << " EneRatio " << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()/getElasticEnergy()
                  << std::endl;
    }

    void addRoller (RollerType rollerType, Mdouble rollerRadius, Mdouble rollerVelocity, Mdouble compaction, Mdouble bladeThickness=50e-6, Mdouble bladeRollerDistance=50e-6)
    {
        rollerType_ = rollerType;
        rollerRadius_ = rollerRadius;
        rollerVelocity_ = rollerVelocity;
        compaction_ = compaction;
        bladeThickness_ = bladeThickness;
        bladeRollerDistance_ = bladeRollerDistance;
    }

private:

    Mdouble relativeBasalLayerThickness_;
    Mdouble rollerRadius_;
    Mdouble rollerVelocity_;
    Mdouble compaction_;
    Mdouble bladeThickness_;
    Mdouble bladeRollerDistance_;
    RollerType rollerType_ = RollerType::NO_ROLLER;
};



int main(int argc UNUSED, char *argv[] UNUSED)
{
    //define domain size
    Mdouble domainLength = 1.0e-3;
    Mdouble domainWidth = 0.25e-3;
    Mdouble domainDepth = 0.1e-3;

    //define properties of particle-particle contacts
    Mdouble density = 1000;
    Mdouble collisionTime = 3e-4;
    Mdouble restitution = 0.2;
    Mdouble friction = 0.5;

    //put above properties into a contact law (do not modify)
    LinearViscoelasticSlidingFrictionSpecies particleSpecies;
    Mdouble mass = density*4./3.*pi*cubic(25e-6);
    particleSpecies.setDensity(density);
    particleSpecies.setCollisionTimeAndRestitutionCoefficient(collisionTime, restitution, mass);
    particleSpecies.setSlidingFrictionCoefficient(friction);
    particleSpecies.setSlidingStiffness(2.0/7.0*particleSpecies.getStiffness());
    particleSpecies.setSlidingDissipation(2.0/7.0*particleSpecies.getDissipation());
    logger(INFO,"Gravitational compression %\% of radius per layer",100*mass*9.8/particleSpecies.getStiffness()/25e-6);

    //create a contact law for particle-wall interactions
    LinearViscoelasticSlidingFrictionMixedSpecies wallSpecies = particleSpecies;

    //Create a solver and run the commands in the constructor PowderBed::PowderBed
    PowderBed pb (domainLength, domainWidth, domainDepth, particleSpecies, wallSpecies);
    pb.setTimeStep(0.05*collisionTime);
    pb.setSaveCount(5.*collisionTime/pb.getTimeStep()); //save every 20 collisions (good for visual output)
    pb.setTimeMax(0.1);

    //Create particles of a certain distribution
    Mdouble meanRadius = 25e-6;
    Mdouble stdRadius = 0.1 * meanRadius;
    pb.setGaussianDistribution(meanRadius,stdRadius);

    //Create particles of a certain distribution
    Mdouble rollerRadius = 250e-6;
    Mdouble rollerVelocity = 0.01; // m/s
    Mdouble compaction = 0.8*domainDepth; // m/s
    pb.addRoller(RollerType::COUNTER_ROTATING_ROLLER,rollerRadius,rollerVelocity,compaction);

    //run initial conditions
    pb.setParticlesWriteVTK(true);
    pb.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    pb.solve();

    return 0;
}
