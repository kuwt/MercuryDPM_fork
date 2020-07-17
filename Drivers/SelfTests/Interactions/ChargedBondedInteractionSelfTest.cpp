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
#include "Species/LinearViscoelasticFrictionChargedBondedSpecies.h"
#include "DPMBase.h"

/** ChargedBondedInteractionSelfTest
 * Two particle pairs are bonded together, then collide.
 * The output is tested with fstatistics.
 */
class ChargedBondedInteractionSelfTest : public DPMBase
{
public:
    void setupInitialConditions() override {
        //set stiffness, radius, density
        double radius = 0.5;
        double density = 6./constants::pi;
        double stiffness = 2e5;  //collision time 1/200
        double dissipation = 200; //for 50 we get restitution coefficient 0.88
        double timeMax = 0.1;
        //set bond strength, such that equilibrium overlap is a fractiono of the radius
        double equilibriumOverlap = 0.1*radius;
        double bondStrength = stiffness*equilibriumOverlap;
        double velocity = 2.0*equilibriumOverlap/(timeMax-0.03);

        //setting the species of particles
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionChargedBondedSpecies());
        //setting the material properties of the particles
        species->setDensity(density);
        species->setStiffness(stiffness);
        species->setBondForceMax(bondStrength);
        species->setBondDissipation(dissipation);

        //set domain size, time max and time step, save count, name
        setMin(-radius*Vec3D(4,1,1));
        setMax(+radius*Vec3D(4,1,1));
        setName("ChargedBondedInteractionSelfTest");
        setTimeStep(1e-5);
        setSaveCount(30);
        setTimeMax(timeMax);

        SphericalParticle p;
        p.setSpecies(species);
        p.setRadius(radius);
        
        //first particle pair
        p.setVelocity(Vec3D(velocity,0,0));
        p.setPosition(Vec3D(-3*radius,0,0));
        particleHandler.copyAndAddObject(p);
        p.setPosition(Vec3D(-1*radius-equilibriumOverlap,0,0));
        particleHandler.copyAndAddObject(p);
        
        //second particle pair
        p.setVelocity(Vec3D(-velocity,0,0));
        p.setPosition(Vec3D(3*radius,0,0));
        particleHandler.copyAndAddObject(p);
        p.setPosition(Vec3D(1*radius+equilibriumOverlap,0,0));
        particleHandler.copyAndAddObject(p);

        // bond particles
        dynamic_cast<ChargedBondedInteraction*>(
                interactionHandler.getInteraction(particleHandler.getObject(0),particleHandler.getObject(1),0))->bond();
        dynamic_cast<ChargedBondedInteraction*>(
                interactionHandler.getInteraction(particleHandler.getObject(2),particleHandler.getObject(3),0))->bond();
    }
};

int main(int argc, char* argv[])
{
    ChargedBondedInteractionSelfTest dpm;
    dpm.solve(argc, argv);

    helpers::writeToFile("ChargedBondedInteractionSelfTestEne.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'Energy [J]'\n"
                         "p 'ChargedBondedInteractionSelfTest.ene' u 1:3 w l t 'kinetic', '' u 1:5 w l t 'potential'\n"
    );
    logger(INFO, "Run gnuplot ChargedBondedInteractionSelfTestTestEne.gnu --persist to view energy statistics");

    helpers::writeToFile("ChargedBondedInteractionSelfTestForces.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'forces [J]'\n"
                         "p 'ChargedBondedInteractionSelfTest.fstat' u 1:7\n"
//                         "p 'ChargedBondedInteractionSelfTest.fstat' u 1:(($2==0)&($3==1)?$9:(1/0))\n"
    );
    logger(INFO, "Run gnuplot ChargedBondedInteractionSelfTestTestForces.gnu --persist to view force statistics");
}
