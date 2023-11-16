//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#include <Mercury3D.h>
#include "Boundaries/CubeInsertionBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Particles/LiquidFilmParticle.h"
#include "Species/LinearViscoelasticFrictionLiquidMigrationWilletSpecies.h"

/**
 * Simulates a very thin rotating drum, to test liquid bridge migration.
 * You can change some of the default parameters via the command line interface. E.g., 
 *   ./RotatingDrumWet -bond 5 -distributeLiquidInhomogeneously
 * sets surface tension such that bond number is 5 and changes to an initially
 * inhomogeneous liquid distribution.
 */
class RotatingDrumWet : public Mercury3D {

private:
    // to allow read-in from command line
    int argc;
    char** argv;
    // store these values as class members
    // because they are set in setupInitialConditions, needed in actionsAfterTimeStep
    bool distributeLiquidInhomogeneously;
    double liquidVolume;
    double liquidBridgeVolumeMax;

public:
    RotatingDrumWet(int argc, char *argv[]) : argc(argc), argv(argv) {}

    void setupInitialConditions() override {
        // dimensional parameters
        double meanParticleRadius = 2e-3; // mean particle radius (m)
        double gravity = 9.8; // gravitational acceleration (m/s^2)
        double density = 1500; // particle density (kg/m^3)

        // dimensionless parameters
        double froude = helpers::readFromCommandLine(argc,argv,"-froude",0.01); // Froude number Fr = w^2R/g (to set drum angular velocity)
        logger.assert_always(froude>0,"Froude number needs to be positive");
        double bond = helpers::readFromCommandLine(argc,argv,"-bond",5.0); // bond number Bo = f_adh/f_grav (to set surface tension)
        logger.assert_always(bond>=0,"Bond number cannot be negative");
        double fillRatio = helpers::readFromCommandLine(argc,argv,"-fillRatio",0.4); // fraction of drum volume to be filled with particles (assumes a bulk solid fraction of 0.55)
        logger.assert_always(fillRatio>0,"Fill ratio needs to be positive");
        double slidingFriction = helpers::readFromCommandLine(argc,argv,"-slidingFriction",0.5); // sliding friction coefficient (-)
        logger.assert_always(slidingFriction>=0,"Sliding friction cannot be negative");
        double rollingFriction = helpers::readFromCommandLine(argc,argv,"-rollingFriction",0.5); // rolling friction coefficient (-)
        logger.assert_always(rollingFriction>=0,"Rolling friction cannot be negative");
        double restitution = helpers::readFromCommandLine(argc,argv,"-restitution",0.5); // restitution coefficient
        logger.assert_always(restitution>0 && restitution<=1,"Restitution has to be 0<r<=1");
        distributeLiquidInhomogeneously = helpers::readFromCommandLine(argc,argv,"-distributeLiquidInhomogeneously"); // whether liquid is initially distributed homogeneously
        double liquidSolidRatio = helpers::readFromCommandLine(argc,argv,"-liquidSolidRatio",0.05);; // volume of liquid over volume of solid (to set initial liquidFilmVolumes)
        logger.assert_always(liquidSolidRatio>=0,"liquid:solid volume ratio cannot be negative");
        double liquidBridgeVolumeMaxRelative = helpers::readFromCommandLine(argc,argv,"-liquidBridgeVolumeMaxRelative",0.15); // max volume of liquid bridge over volume of mean particle; used to set liquidBridgeVolumeMax; set to 1/6, assuming porosity of 4/9 and 4 times more contacts than particles.
        logger.assert_always(liquidBridgeVolumeMaxRelative>=0,"LiquidBridgeVolumeMax cannot be negative");

        // dependent parameters
        double drumRadius = helpers::readFromCommandLine(argc,argv,"-drumRadiusRelative",20.0)*meanParticleRadius; // drum radius (m)
        double drumWidth = 5*meanParticleRadius; // drum width (m)
        double angularVelocity = sqrt(froude*gravity/drumRadius); // rotation speed (rad/s)
        double timeMax = 10.0*constants::pi/angularVelocity; //simulation time (s), based on 5 rotations
        double stdParticleRadius = 0.1*meanParticleRadius; // standard deviation of particle radius (m)
        double solidVolume = constants::pi*mathsFunc::square(drumRadius)*drumWidth*fillRatio*0.55; // particle volume (m^3)
        liquidVolume = liquidSolidRatio * solidVolume; // liquid volume (m^3)
        logger(INFO,"drumRadius %, drumWidth %, angularVelocity %, timeMax %", drumRadius, drumWidth, angularVelocity, timeMax);

        // set name of output files
        setName("RotatingDrumWet");

        // turn on paraview output
        setParticlesWriteVTK(true);
        wallHandler.setWriteVTK(FileType::ONE_FILE);
        interactionHandler.setWriteVTK(true);

        // turn on gravity
        setGravity(Vec3D(0,0,-gravity));

        // set domain for visualisation
        setMin(Vec3D(-drumRadius, -0.5*drumWidth, -drumRadius));
        setMax(Vec3D(drumRadius, 0.5*drumWidth, drumRadius));

        // material and contact properties
        auto speciesParticle = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionLiquidMigrationWilletSpecies());
        speciesParticle->setDensity(density);
        double mass = speciesParticle->getMassFromRadius(meanParticleRadius);
        // collision time is set significantly smaller than gravitational time scale (factor of 0.1 is chosen to keep overlaps to 1% of particle radius)
        double collisionTime = 0.1*sqrt(meanParticleRadius/gravity);
        speciesParticle->setCollisionTimeAndRestitutionCoefficient(collisionTime, restitution, mass);
        speciesParticle->setSlidingFrictionCoefficient(slidingFriction);
        speciesParticle->setSlidingStiffness(2. / 7. * speciesParticle->getStiffness());
        speciesParticle->setSlidingDissipation(2. / 7. * speciesParticle->getDissipation());
        speciesParticle->setRollingFrictionCoefficient(rollingFriction);
        speciesParticle->setRollingStiffness(2. / 7. * speciesParticle->getStiffness());
        speciesParticle->setRollingDissipation(2. / 7. * speciesParticle->getDissipation());
        // Bo = f_adh/f_grav, f_adh = -2*pi*effectiveRadius*surfaceTension, f_grav = m*g
        speciesParticle->setSurfaceTension(bond * mass * gravity / (2.0 * constants::pi * meanParticleRadius));
        //s->setLiquidBridgeVolumeMin(0.01*liquidBridgeVolumeMaxRelative*mass/density);
        liquidBridgeVolumeMax = liquidBridgeVolumeMaxRelative * mass / density;
        speciesParticle->setLiquidBridgeVolumeMax(liquidBridgeVolumeMax);
        //s->setDistributionCoefficient(0.0);
        logger(INFO, "collisionTime %, surfaceTension %", speciesParticle->getCollisionTime(mass), speciesParticle->getSurfaceTension());

        // wall species is equal to particle species, except no liquid bridges can form
        auto speciesWall = speciesHandler.copyAndAddObject(*speciesParticle);
        auto speciesParticleWall = speciesHandler.getMixedObject(speciesParticle, speciesWall);
        speciesParticleWall->setLiquidBridgeVolumeMax(0);


        // particle size distribution
        PSD psd = PSD::getDistributionNormal(meanParticleRadius, stdParticleRadius, 50);

        // set final time, time step, and output frequency
        setTimeMax(timeMax);
        setTimeStep(0.05*collisionTime); //! this timestep is a bit too high
        setSaveCount(2.0*constants::pi/angularVelocity/getTimeStep()/100); // produce 100 output files per rotation

        // add cylindrical wall
        AxisymmetricIntersectionOfWalls cylinder;
        cylinder.setSpecies(speciesWall); // set material and contact properties
        cylinder.setAxis(Vec3D(0.0,1.0,0.0)); // axis of rotation
        cylinder.addObject(Vec3D(1.0,0.0,0.0),Vec3D(drumRadius,0.0,0.0));
        cylinder.setAngularVelocity(Vec3D(0.0,angularVelocity,0.0));
        wallHandler.copyAndAddObject(cylinder);

        // add periodic sides (so no effect of side walls)
        PeriodicBoundary boundary;
        boundary.set(Vec3D(0.,1.,0.),getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(boundary);

        // compute necessary particle number
        int particleNumber = solidVolume/(mass/density); // number of particles in drum (-)
        logger(INFO, "particleNumber %", particleNumber);
        double liquidBridgeVolumeMean = liquidVolume/particleNumber;
        logger(INFO, "liquidFilmVolumeMin %, liquidBridgeVolumeMean %, liquidBridgeVolumeMax %, liquidBridgeRadiusMax %", 
               speciesParticle->getLiquidBridgeVolumeMin(), 
               liquidBridgeVolumeMean, 
               speciesParticle->getLiquidBridgeVolumeMax(), 
               cbrt(speciesParticle->getLiquidBridgeVolumeMax())
               );

        // type of particles to insert
        LiquidFilmParticle p;
        p.setSpecies(speciesParticle);
        if (!distributeLiquidInhomogeneously) {
            p.setLiquidVolume(liquidVolume/particleNumber);
        }

        // insert dry particles
        CubeInsertionBoundary insertion;
        insertion.set(&p, 1000, getMin(), getMax(), Vec3D(0, 0, 0), Vec3D(0, 0, 0));
        insertion.setPSD(psd);
        insertion.setInitialVolume(solidVolume);
        insertion.insertParticles(this); // insert first batch of liquid
        boundaryHandler.copyAndAddObject(insertion);

        // remove old output files
        removeOldFiles();
    }

    void printTime() const override
    {
        Vec3D com = particleHandler.getCentreOfMass();
        double angleOfRepose = atan2(-com.X,-com.Z)*180.0/constants::pi;
        logger(INFO,"t %\tN %\tNc %\tNl %\to %\tAoR %",
               getTime()/getTimeMax(),
               particleHandler.getNumberOfObjects(),
               interactionHandler.getNumberOfObjects(),
               interactionHandler.getNumberOfLiquidBridges(),
               interactionHandler.getMeanOverlap()/particleHandler.getMeanRadius(),
               angleOfRepose);
    }

    // liquid is inserted after the third savecount, in clumps proportional to the liquidBridgeVolumeMax, so it has to re-distribute in the material
    void actionsAfterTimeStep() {
        if (distributeLiquidInhomogeneously && getNumberOfTimeSteps()==3*dataFile.getSaveCount()) {
            logger(INFO,"Distributing liquid");
            double liquidDistributed = 0;
            LiquidFilmParticle* p;
            do {
                int id = rand() % particleHandler.getSize();
                p = static_cast<LiquidFilmParticle*>(particleHandler.getObject(id));
                p->setLiquidVolume(liquidBridgeVolumeMax);
                liquidDistributed += liquidBridgeVolumeMax;
            } while (liquidDistributed<liquidVolume);
            p->addLiquidVolume(liquidDistributed-liquidVolume);
        }
    }
};

//! Defines a rotating drum simulation
int main(int argc, char *argv[]) {
    RotatingDrumWet drum(argc, argv);
    drum.solve();
    return 0;
}

