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

#include<iostream>
#include "scr/DPMBase.h"
#include "scr/Mercury3D.h"
#include "scr/Chute.h"
using namespace std;

class ChutePeriodic : public Chute{
public:

	void actionsBeforeTimeStep(){};
		
	void setupInitialConditions(){
		//load initial particles
		readDataFile("c3d.ini") || readDataFile("../c3d.ini");
		//fix Particles with zero velocity
		int count=0;
		for (vector<CParticle>::iterator it = Particles.begin(); it!=Particles.end(); ++it) {
			if (it->Velocity.X==0.0 && it->Velocity.Y==0.0 && it->Velocity.Z==0.0) {
				it->fixParticle();
				count++;
			}
		}
		//set timestep to an eigth
		setTimeStep(getTimeStep()/8);
		//save every .2 time units
		setSaveCount(floor(0.25/getTimeStep()));
		//set number of buckets to power of two
		int NUM_BUCKETS = 1024;
		while (NUM_BUCKETS<2*get_N()) {
			NUM_BUCKETS *= 2;
		}
		set_HGRID_num_buckets(NUM_BUCKETS);
		//output
		cout << "fixed particles: " << count << endl;
		cout << "dt=" << getTimeStep() << ", savecount=" << get_savecount()  << ", NUM_BUCKETS=" << NUM_BUCKETS << endl; 
	  write(std::cout,false);
	};

};

int main(int argc, char *argv[])
{
	if (argc<2) {cerr << "Arguments needed!" << endl; exit(-1);}

 	ChutePeriodic problem;
 	problem.setName(argv[1]);
	problem.load_restart_data();
	problem.solve();
	problem.write(std::cout,false);
	problem.writeRestartFile();
}
