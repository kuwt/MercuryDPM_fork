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
#include <Species/LinearViscoelasticSpecies.h>
#include "Chute.h"

class DPM : public Mercury3D {

    void actionsAfterSolve() override {
        auto p = particleHandler.getLastObject();
        p->setPosition({2, 0, 0});
        hGridUpdateParticle(p);
        logger(INFO, "Particle p at x=% is at hGrid x-position %", p->getPosition().X, p->getHGridX());
        p->setPosition({2.0000001, 0, 0});
        hGridUpdateParticle(p);
        logger(INFO, "Particle p at x=% is at hGrid x-position %", p->getPosition().X, p->getHGridX());
    }

    void actionsBeforeTimeStep() override {
        auto p = particleHandler.getLastObject();
        static Mdouble disp0 = getHGridTotalCurrentMaxRelativeDisplacement();
        Mdouble disp = getHGridTotalCurrentMaxRelativeDisplacement();
        if (disp<disp0) {
            logger(INFO,"HGrid has been updated after t=% (x=%), displacement %",getTime()-getTimeStep(),p->getPosition().X-p->getVelocity().X*getTimeStep(),disp0);
        }
        disp0 = disp;
        static int x0 = p->getHGridX();
        int x = p->getHGridX();
        if (x0!=x) {
            logger(INFO,"Particle p at x=% has been updated to hGrid x-position % at t=%",p->getPosition().X-p->getVelocity().X*getTimeStep(), x, getTime()-getTimeStep());
            x0 = x;

            if (x==1) {
                helpers::check(p->getPosition().X-p->getVelocity().X*getTimeStep(),2.0005,1e-10,"x-position");
            } else if (x==2) {
                helpers::check(p->getPosition().X-p->getVelocity().X*getTimeStep(),4.0009,1e-10,"x-position");
            } else if (x==3) {
                helpers::check(p->getPosition().X-p->getVelocity().X*getTimeStep(),6.0013,1e-10,"x-position");
            } else if (x==4) {
                helpers::check(p->getPosition().X-p->getVelocity().X*getTimeStep(),8.0017,1e-10,"x-position");
            }
        }
    }
};


int main()
{

    // Problem parameters
    DPM problem;
    problem.setFileType(FileType::NO_FILE);
    problem.setHGridCellOverSizeRatio(2);
    problem.setHGridUpdateEachTimeStep(false);
    problem.setName("HGridUpdateUnitTest");   // data output file name
    problem.setTimeStep(1e-4);
    problem.setTimeMax(10);
    problem.setDomain({0,0,0},{1,1,1});

    // Particle species
    LinearViscoelasticSpecies species;              // initialise species
    species.setDensity(2000);                       // particle density
    species.setCollisionTimeAndRestitutionCoefficient( 5e-3, 0.8, 1);
    auto s = problem.speciesHandler.copyAndAddObject(species);   // assign species to problem species handler

    SphericalParticle particle;
    particle.setSpecies(s);
    particle.setRadius(0.5);
    particle.setVelocity({1,0,0});
    auto p = problem.particleHandler.copyAndAddObject(particle);

    //logger(INFO,"Cell %",)


    //solve
    problem.solve();
    problem.hGridInfo();
    return 0;
}
