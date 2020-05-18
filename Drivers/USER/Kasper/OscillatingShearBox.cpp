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

#include <iostream>

#include "Mercury3D.h"
#include "Boundaries/ShearBoxBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "StatisticsVector.h"
#include "Walls/InfiniteWall.h"

/**
 * Defines a shear box with oscillating boundary conditions.
 *
 * Specifically, for the periodic boundary in x, all the velocity of all ghost particles is increased by A*sin(omega*t)*(z/zMean-1). Thus, particles get accellerated periodically to the left and right.
 */
class OscillatingShearBox : public Mercury3D
{
public:

    //uses the user-set parameters to define a Lees-Edwards shear box with oscillating boundary conditions.
    void setupInitialConditions() override
    {
        setMin({0,0,0});
        setMax({boxWidth,boxLength,boxHeight});
		
        setName("OscillatingShearBox");
        setSaveCount(50);
        setXBallsAdditionalArguments(" -v0 -solidf");
        setGravity({0,0,-1});

        logger(INFO,"Adding species");
        auto species =speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        species->setDensity(particleMass/species->getVolumeFromRadius(particleRadius));
        species->setCollisionTimeAndRestitutionCoefficient(collisionTime,restitutionCoefficient,particleMass);
        species->setSlidingStiffness(2.0/7.0*species->getStiffness());
        species->setSlidingDissipation(2.0/7.0*species->getDissipation());
        species->setSlidingFrictionCoefficient(slidingFrictionCoefficient);
        setTimeStep(timeStep);

        logger(INFO,"Adding Lees Edwards bc in x and periodic walls in y");
        Mdouble angularFrequency = 2.0*constants::pi*oscillationFrequency;
        auto shearBoxBoundary = boundaryHandler.copyAndAddObject(ShearBoxBoundary());
        shearBoxBoundary->set(
         /*velocity*/[=] (double time, double z) {
             return oscillationAmplitude*std::sin(angularFrequency*time)*(z/4.0-1);
         }, getXMin(),getXMax(),getYMin(),getYMax());
        wallHandler.copyAndAddObject(InfiniteWall(
         Vec3D(0.0,0.0,-1.0), Vec3D(0.0,0.0,getZMin()), species
        ));
        wallHandler.copyAndAddObject(InfiniteWall(
         Vec3D(0.0,0.0,1.0), Vec3D(0.0,0.0,getZMax()), species
        ));


        unsigned N=bulkVolumeFraction*boxHeight*boxLength*boxWidth/species->getVolumeFromRadius(particleRadius);
        logger(INFO,"Adding % particles", N);
        SphericalParticle p;
        p.setSpecies(species);
        p.setRadius(particleRadius);
        Vec3D position;
        unsigned i=0;
        for (unsigned i=0; i<N; i++){
            position.X = random.getRandomNumber(getXMin(),getXMax());
            position.Y = random.getRandomNumber(getYMin(),getYMax());
            position.Z = random.getRandomNumber(getZMin()+particleRadius,getZMax()-particleRadius);
            p.setPosition(position);
            particleHandler.copyAndAddObject(p);
        }
    }

    //overrides printTime such that console output shows the state of relaxation (Ekin/Eela)
    void printTime() const override
    {
        std::cout 
            << "t=" << std::setw(12) << getTime() 
            << " tmax=" << std::setw(12) << getTimeMax() 
            << " eneRatio " << std::setw(12) << getKineticEnergy()/getElasticEnergy()
            << std::endl;
    }

public:
    //parameters set by the user
    Mdouble particleRadius;
    Mdouble collisionTime;
    Mdouble particleMass;
    Mdouble timeStep;
    Mdouble restitutionCoefficient;
    Mdouble slidingFrictionCoefficient;
    Mdouble bulkVolumeFraction;
    Mdouble boxLength;
    Mdouble boxWidth;
    Mdouble boxHeight;
    Mdouble oscillationFrequency;
    Mdouble oscillationAmplitude;
};


int main(int argc UNUSED, char *argv[] UNUSED)
{
    OscillatingShearBox le;
    le.particleRadius = 0.5;
    le.collisionTime = 0.01;
    le.particleMass = 1.0;
    le.timeStep = 0.02*le.collisionTime;
    le.restitutionCoefficient = 0.1;
    le.slidingFrictionCoefficient = 0.1;
    le.bulkVolumeFraction = 0.6;
    le.boxLength = 8;
    le.boxWidth = 10;
    le.boxHeight = 8;
    le.oscillationFrequency = 1.0;
    le.oscillationAmplitude = 4.0;
    le.setTimeMax(2.0*le.oscillationFrequency);
    le.solve();
}
