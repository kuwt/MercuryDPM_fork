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
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include "Boundaries/CubeInsertionBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"

class Drum : public Mercury3D {
public:
    void setupInitialConditions() override {
        // parameters
        double drumRadius = 40e-3; // drum radius (m)
        double drumWidth = 40e-3; // drum width (m)
        double angularVelocity = 2.0*constants::pi; // rotation speed (rad/s)
        double timeMax = 1; //simulation time (s)
        double meanParticleRadius = 2e-3; // mean particle radius (m)
        double stdParticleRadius = 0.2e-3; // standard deviation of particle radius (m)
        double stiffness = 64; // stiffness of particle contacts (N/m
        double dissipation = 10e-3; // dissipation of particle contacts (Ns/m)
        double density = 1500; // particle density (kg/m^3)
        double slidingFriction = 0.5; // sliding friction coefficient (-)
        double rollingFriction = 0.5; // rolling friction coefficient (-)
        int particleNumber = 1000; // number of particles in drum (-)

        // set name of output files
        setName("RotatingDrum");

        // remove old output files
        removeOldFiles();

        // turn on paraview output
        setParticlesWriteVTK(true);
        wallHandler.setWriteVTK(true);

        // turn on gravity
        setGravity(Vec3D(0,0,-9.8));

        // set domain for visualisation
        setMin(Vec3D(-drumRadius, -0.5*drumWidth, -drumRadius));
        setMax(Vec3D(drumRadius, 0.5*drumWidth, drumRadius));

        // material and contact properties
        LinearViscoelasticFrictionSpecies species;
        species.setDensity(density);
        species.setStiffness(stiffness);
        species.setDissipation(dissipation);
        species.setSlidingFrictionCoefficient(slidingFriction);
        species.setSlidingStiffness(2./7.*stiffness);
        species.setSlidingDissipation(2./7.*dissipation);
        species.setRollingFrictionCoefficient(rollingFriction);
        species.setRollingStiffness(2./7.*stiffness);
        species.setRollingDissipation(2./7.*dissipation);
        auto s = speciesHandler.copyAndAddObject(species);

        //
        double collisionTime = s->getCollisionTime(s->getMassFromRadius(meanParticleRadius));
        double restitution = s->getRestitutionCoefficient(s->getMassFromRadius(meanParticleRadius));
        logger(INFO,"collisionTime %, restitution coefficient %", collisionTime, restitution);

        // particle size distribution
        PSD psd = PSD::getDistributionNormal(meanParticleRadius, stdParticleRadius, 50);

        // set final time, time step, and output frequency
        setTimeMax(timeMax);
        setTimeStep(0.02*collisionTime);
        setSaveCount(getTimeMax()/getTimeStep()/100);

        // add cylindrical wall
        AxisymmetricIntersectionOfWalls cylinder;
        cylinder.setSpecies(s); // set material and contact properties
        cylinder.setAxis(Vec3D(0.0,1.0,0.0)); // axis of rotation
        cylinder.addObject(Vec3D(1.0,0.0,0.0),Vec3D(drumRadius,0.0,0.0));
        cylinder.setAngularVelocity(Vec3D(0.0,angularVelocity,0.0));
        wallHandler.copyAndAddObject(cylinder);

        //add side walls
        InfiniteWall side;
        side.setSpecies(s); // set material and contact properties
        side.setAngularVelocity(Vec3D(0.0,angularVelocity,0.0));
        side.set(Vec3D(0.,-1.,0.),Vec3D(0.0,-0.5*drumWidth,0.0));
        wallHandler.copyAndAddObject(side);
        side.set(Vec3D(0.,1.,0.),Vec3D(0.0,0.5*drumWidth,0.0));
        wallHandler.copyAndAddObject(side);

        //add type of particles to insert
        SphericalParticle p;
        p.setSpecies(s);

        // insert particles
        CubeInsertionBoundary insertion;
        insertion.set(&p, 1000, getMin(), getMax(), Vec3D(0, 0, 0), Vec3D(0, 0, 0));
        insertion.setPSD(psd);
        insertion.setInitialVolume(particleNumber*s->getVolumeFromRadius(meanParticleRadius));
        boundaryHandler.copyAndAddObject(insertion);
    }

    void printTime() const override
    {
        logger(INFO,"t %\tN %\to %",
               getTime()/getTimeMax(),
               particleHandler.getNumberOfObjects(),
               interactionHandler.getMeanOverlap()/particleHandler.getMeanRadius());
    }
};

//! Defines a rotating drum simulation
int main() {
    Drum drum;
    drum.solve();
    return 0;
}

