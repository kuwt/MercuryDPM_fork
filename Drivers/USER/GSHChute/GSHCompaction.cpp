//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include <Boundaries/PeriodicBoundary.h>
#include "Mercury3D.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
using constants::pi;
using mathsFunc::cubic;

/*
 * Simulates a box of particles compressed step-by-step to high volume fraction
 */
class GSHCompaction : public Mercury3D{

public:
    
    void setupInitialConditions() override {
        setName("GSHCompaction");
        setTimeStep(1e-4);
        setTimeMax((initialDomainWidth-finalDomainWidth)/compactionVelocity);
        setMin(Vec3D(0,0,0));
        setMax(initialDomainWidth*Vec3D(1,1,1));
        setSaveCount(getTimeMax()/getTimeStep()/100);
        fStatFile.writeFirstAndLastTimeStep();
        restartFile.setFileType(FileType::MULTIPLE_FILES);
        setXBallsAdditionalArguments("-v0 -solidf");
        helpers::writeToFile(getName()+".gnu",
                "set xlabel 'volumeFraction'; set ylabel 'pressure'; plot 'GSHCompaction.out' using 1:3");

        PeriodicBoundary boundary;
        periodicBoundary[0] = boundaryHandler.copyAndAddObject(boundary);
        periodicBoundary[1] = boundaryHandler.copyAndAddObject(boundary);
        periodicBoundary[2] = boundaryHandler.copyAndAddObject(boundary);
        setPeriodicBoundaries();

        LinearViscoelasticSlidingFrictionSpecies species;
        species.setDensity(6.0/pi);
        species.setStiffness(2e5);
        species.setDissipation(25);
        species.setSlidingFrictionCoefficient(0.5);
        species.setSlidingDissipation(2./7.*species.getDissipation());
        species.setSlidingStiffness(2./7.*species.getStiffness());
        auto s = speciesHandler.copyAndAddObject(species);

        CubeInsertionBoundary insertionBoundary;
        insertionBoundary.set(BaseParticle(s),NEVER,getMin(),getMax(),Vec3D(0,0,0),Vec3D(0,0,0),0.475,0.525);
        insertionBoundary.setInitialVolume((double)n*pi/6.0);
        insertionBoundary.insertParticles(this);
    }

    void actionsBeforeTimeStep() override {
        setPeriodicBoundaries();
    }

    void computeInternalForces(BaseParticle* particle) override {
        particle->addForce(-25*particle->getVelocity());
        DPMBase::computeInternalForces(particle);
    }

    void printTime() const override {
        static std::ofstream outFile(getName()+".out");
        Mdouble solidFraction = cubic(finalDomainWidth/getXMax());
        Mdouble kineticPressure = getKineticStress().trace()/3.0;
        Mdouble staticPressure = getStaticStress().trace()/3.0;
        Mdouble pressure = kineticPressure + staticPressure;
        logger(INFO,"solidFraction % Pressure % Kinetics %", solidFraction, pressure, kineticPressure/staticPressure);
        outFile << solidFraction << ' ' << pressure << ' ' << staticPressure << std::endl;
    }

    void setPeriodicBoundaries() {

        Mdouble shrinkFactor = (initialDomainWidth-compactionVelocity*getTime())/getXMax();

        setMax(shrinkFactor*getMax());

        for (const auto particle : particleHandler) {
            particle->setPosition(shrinkFactor*particle->getPosition());
        }

        periodicBoundary[0]->set(Vec3D(1,0,0),getXMin(),getXMax());
        periodicBoundary[1]->set(Vec3D(0,1,0),getYMin(),getYMax());
        periodicBoundary[2]->set(Vec3D(0,0,1),getZMin(),getZMax());
    }

    //number of particles
    unsigned n = 1000;
    //compaction velocity
    Mdouble compactionVelocity = 0.5;

private:

    //stores pointers to the periodic boundaries
    std::array<PeriodicBoundary*,3> periodicBoundary;
    //stores initial and final domain width
    Mdouble initialDomainWidth = cbrt((double)n*1.8);
    Mdouble finalDomainWidth = cbrt((double)n*pi/6.0);
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    GSHCompaction compaction;
	compaction.solve();
}
