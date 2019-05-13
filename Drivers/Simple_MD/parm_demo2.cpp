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

#include "DPMBase.h"
#include <iostream>
#include <vector>

class two_particle_elastic_collision : public DPMBase {

	void setupInitialConditions()
	{
		
		//Check if the run has been done before
		if (FileExists(data_filename.str()))
		{
			//If it has move on to teh next run immedently
			std::cout << "This run has been done " << std::endl;
			launchNewRun("parm_demo",true);
			exit(0);
		}
		
		SphericalParticle P0,P1;
		
		particleHandler.copyAndAddObject(P0);
		particleHandler.copyAndAddObject(P1);
		
		//Set up the stuff that is the same for all runs
		particleHandler.getObject(0)->setPosition(Vec3D(0.006,0.005,0.0));
		particleHandler.getObject(1)->setPosition(Vec3D(0.004,0.005,0.0));
	
		particleHandler.getObject(0)->setVelocity(Vec3D(-0.1,0.0,0.0));
		particleHandler.getObject(1)->setVelocity(Vec3D( 0.1,0.0,0.0));
	
		//Set up a 2 by 3 study
		vector<int> study_num=get2DParametersFromRunNumber(2,3);
		
		//We have three different size particles for both left and right particles
		particleHandler.getObject(0)->setRadius(0.0001*study_num[1]);
		particleHandler.getObject(1)->setRadius(0.0001*study_num[2]);
		
		//If study 0 is complete quit
		if (study_num[0] > 0)
			{
				std::cout << "Whole study has started" << std::endl;
				exit(0);
			}
		else
		//If the study is not complete save the data to disk and move on
			{
				writeRestartFile();
				launchNewRun("parm_demo");
			}
		
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	

	///Start off by solving the default problem
	two_particle_elastic_collision problem;
	
	///Autonumber turns on file numbering
	problem.autoNumber();
	problem.set_name("parm_demo"); //elastic_collision
	problem.setTimeStep(1e-5);
	problem.setTimeMax(1.0);
	problem.setSaveCount(10);
	
	problem.solve();
	
	
	
	
}
