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

/// In this file 32^2 particles with the same velocity are placed in a bi-axial box. This makes them collide with the

class GSHRelax : public Mercury3D{
public:

    GSHRelax(int i) {
        readRestartFile("GSHCompaction.restart."+helpers::to_string(i));
        setTimeMax(inf);
        setFileType(FileType::NO_FILE);
        restartFile.setFileType(FileType::ONE_FILE);
        restartFile.setName("GSHRelax.restart."+i);
    }

    void computeInternalForces(BaseParticle* particle) override {
        particle->addForce(-25*particle->getVelocity());
        DPMBase::computeInternalForces(particle);
    }

    bool continueSolve() const override {
        Mdouble eneKin = getKineticEnergy();
        Mdouble eneEla = getElasticEnergy();
        static int counter = 10000;
        if (++counter>9999) {
            logger(INFO,"EneKin % EneEla %", eneKin, eneEla);
            counter=0;
        }
        return (eneKin>std::max(1e-7,1e-5*eneEla));
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    std::ofstream outFile("GSHRelax.out");

    for (int i=0; i<100; i++) {
        GSHRelax relax(i);
        relax.solve();
        Mdouble solidFraction = relax.particleHandler.getVolume()/relax.getTotalVolume();
        Mdouble kineticPressure = relax.getKineticStress().trace()/3.0;
        Mdouble staticPressure = relax.getStaticStress().trace()/3.0;
        outFile << solidFraction << ' ' << kineticPressure << ' ' << staticPressure << std::endl;
    }
}
