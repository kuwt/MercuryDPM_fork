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

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include "Species/Species.h"
#include "Species/LinearViscoelasticFrictionBondedSpecies.h"
using constants::pi;
using mathsFunc::cubic;

/**
 * In this file two bonded particle pairs are placed in a box, and are allowed to jump around under gravity.
 */
class TorsionFrictionSelfTest : public Mercury3D {

	void setupInitialConditions() override
	{
        setName("TorsionFrictionSelfTest");

        // make species such
        // - A particle of diameter 1 has a mass of 1
        // - Elastic force, no energy loss by dissipation
        // - Bond force such that two bonded particles are in equilibrium if delta=0.1, f_n=1
        // - Torsion friction, such that the max torque is mut*f_n*r=0.1
        LinearViscoelasticFrictionBondedSpecies species;
        species.setDensity(6.0/pi);
        //species.setCollisionTimeAndRestitutionCoefficient(1,1,1);
        species.setStiffness(10);
        species.setTorsionStiffness(4);
        species.setTorsionFrictionCoefficient(.2);
        species.setBondForceMax(species.getStiffness()*0.1);
        // Use s to set the species of particles and walls
        auto s = speciesHandler.copyAndAddObject(species);

        //set time-stepping properties properties to simulate 10 collision times
        setSaveCount(1);
        setTimeStep(0.02*s->getCollisionTime(1.0));
        setTimeMax(4.0*s->getCollisionTime(1.0));

        //set system dimension
		setMax(Vec3D(1,1,1));
		setMin(-getMax());

		//define two particles at equilibrium overlap twisting around each other
		SphericalParticle particle;
        particle.setSpecies(s);
		particle.setPosition(Vec3D(0,0,+0.45));
		particle.setRadius(0.5);
		auto p0 = particleHandler.copyAndAddObject(particle);
		particle.setPosition(Vec3D(0,0,-0.45));
		particle.setAngularVelocity(Vec3D(0,0,1));
        auto p1 = particleHandler.copyAndAddObject(particle);

        //bond the two particles particles
		dynamic_cast<BondedInteraction*>(interactionHandler.getInteraction(p0,p1,0))->bond();
	}

    void actionsAfterTimeStep() override
    {
        Mdouble torque = dynamic_cast<BondedInteraction*>(interactionHandler.getLastObject())->getTorque().Z;
        static std::ofstream ofile ("TorsionFrictionSelfTest.torque");
        ofile << getTime() << ' ' << torque << '\n';
    }

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	//create an instance of the DPM class
    TorsionFrictionSelfTest dpm;
    dpm.solve();

    helpers::writeToFile("TorsionFrictionSelfTest.gnu",
            "set xlabel 'time'\n"
            "set ylabel 'torque'\n"
            "p 'TorsionFrictionSelfTest.torque' u 1:2 t 'TorsionFrictionSelfTest'\n"
            );

    helpers::writeToFile("TorsionFrictionSelfTestEne.gnu",
                         "set xlabel 'time'\n"
                         "set ylabel 'Potential + Kinetic Energy'\n"
                         "p 'TorsionFrictionSelfTest.ene' u 1:($4+$5+$6) t ''\n"
    );

    return 0;
}
