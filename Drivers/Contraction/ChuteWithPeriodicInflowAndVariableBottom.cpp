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

#include "ChuteWithPeriodicInflow.h"

class ChuteWithPeriodicInflowAndVariableBottom : public ChuteWithPeriodicInflow
{
public:
  ChuteWithPeriodicInflowAndVariableBottom(std::string restart_file) : ChuteWithPeriodicInflow(restart_file) {
		// creates bottom outside periodic chute of species 1
		int numRepetitions = 3; //total length 30x20=600
		double lengthPeriodicChute = getXMax();
		//set_Nmax(get_N()*(numRepetitions+2.));
		int N=particleHandler.getNumberOfObjects();//get_N();
		int j;
		for (j=1; j<std::max(2.,numRepetitions*.15); j++) {
			for (int i=0; i<N; i++) {
				if (inflowParticle_.isFixed()) {
				  	  particleHandler.addObject(particleHandler.getObject(i));//Particles.push_back(Particles[i]);
					  particleHandler.getLastObject()->setSpecies(speciesHandler.getObject(1));
					  particleHandler.getLastObject()->move(Vec3D(lengthPeriodicChute*j,0.,0.));//Particles.back().Position.X+=lengthPeriodicChute*j;

				  //Particles.push_back(Particles[i]);
				  //Particles.back().indSpecies=1;
				  //Particles.back().Position.X+=lengthPeriodicChute*j;
				}
			}
		}

		Chute SlowBottom;
		SlowBottom.setName("../ChuteWithPeriodicInflow/H20A22L1M0.5B0.5");
		SlowBottom.readRestartFile();
		for (; j<numRepetitions; j++) {
			for (int i=0; i<SlowBottom.particleHandler.getNumberOfObjects(); i++) {
			  if (SlowBottom.particleHandler.getObject(i)->isFixed()) {
				  					  
			                  particleHandler.addObject(particleHandler.getObject(i));//Particles.push_back(Particles[i]);
					  particleHandler.getLastObject()->setSpecies(speciesHandler.getObject(1));
					  particleHandler.getLastObject()->move(Vec3D(lengthPeriodicChute*j,0.,0.));//Particles.back().Position.X+=lengthPeriodicChute*j;

					  //Particles.push_back(SlowBottom.getObjects()[i]);
					  //Particles.back().indSpecies=1;
					  //Particles.back().Position.X+=lengthPeriodicChute*j;
				}
			}
		}
		setXMax(numRepetitions*getXMax());
		//set_HGRID_num_buckets_to_power(get_Nmax());
		//Walls.resize(0);
		setName("ChuteWithPeriodicInflowAndVariableBottom");
	}
};



int main(int, char**)
{
	ChuteWithPeriodicInflowAndVariableBottom problem("../ChuteWithPeriodicInflow/H20A22L0.5M0.5B0.5");
	problem.setTimeMax(1e4);
	problem.setSaveCount(1000);
	problem.setXBallsAdditionalArguments("-v0 -solidf");
	problem.write(std::cout,false);
	problem.dataFile.setFileType(FileType::ONE_FILE);
	problem.restartFile.setFileType(FileType::ONE_FILE);
	problem.fStatFile.setFileType(FileType::ONE_FILE);
	problem.eneFile.setFileType(FileType::ONE_FILE);
	//problem.restartFile.setFileType(FileType::ONE_FILE);
	//problem.dataFile.setFileType(FileType::ONE_FILE);
	//problem.fStatFile.setFileType(FileType::NO_FILE);
	problem.solve();
}
