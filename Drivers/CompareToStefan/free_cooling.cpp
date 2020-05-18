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

#include<iostream>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include "Mercury2D.h"
class free_cooling : public Mercury2D{
public:

	void actionsBeforeTimeStep(){};
		
	void setupInitialConditions()
	{
        ///\todo TW check that all user codes are compiling AND running; but many of them need input files.
        if (readDataFile("c3d.ini",7) || readDataFile("../c3d.ini",7)) {
			//~ for (vector<CParticle>::iterator it = Particles.begin(); it!=Particles.end(); ++it) {
				//~ it->Position.Y =  it->Position.Z;
				//~ it->Position.Z = 0.0;
				//~ it->Velocity.Y =  it->Velocity.Z;
				//~ it->Velocity.Z = 0.0;
			//~ }
			//setStiffnessAndRestitutionCoefficient(1e6, .9, Particles[0].get_mass());
		    write(std::cout,true);
            auto species = dynamic_cast<LinearViscoelasticFrictionSpecies*>(speciesHandler.getObject(0));
            if (species!= nullptr)
            {
                std::cout << "Mass: " << particleHandler.getObject(0)->getMass() << std::endl;
                std::cout << "Collision time: " << species->getCollisionTime(particleHandler.getObject(0)->getMass()) << std::endl;
                std::cout << "Restitution coefficient: " << species->getRestitutionCoefficient(particleHandler.getObject(0)->getMass()) << std::endl;
            }
		} else {
		  std::cerr << "Input data not found exiting " << std::endl;
			exit(-1);
		}
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
 	free_cooling problem;
 	problem.setName("c3d");
	problem.readRestartFile();//load_restart_data();
	problem.setSaveCount(1);
	problem.setTimeMax(4e-3); 	
	problem.solve();
	
	problem.write(std::cout,false);
	problem.writeRestartFile();
}
